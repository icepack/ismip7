r"""Regridding from Firedrake unstructured mesh to ISMIP7 regular grid."""

import numpy as np

ISMIP_VARIABLES = {
    "lithk":        {"long_name": "Land ice thickness", "units": "m"},
    "orog":         {"long_name": "Surface altitude", "units": "m"},
    "base":         {"long_name": "Base altitude", "units": "m"},
    "topg":         {"long_name": "Bedrock altitude", "units": "m"},
    "acabf":        {"long_name": "Surface mass balance flux", "units": "kg m-2 s-1"},
    "libmassbfgr":  {"long_name": "Basal mass balance flux (grounded)", "units": "kg m-2 s-1"},
    "libmassbffl":  {"long_name": "Basal mass balance flux (floating)", "units": "kg m-2 s-1"},
    "dlithkdt":     {"long_name": "Tendency of land ice thickness", "units": "m s-1"},
    "xvelsurf":     {"long_name": "Surface x-velocity", "units": "m s-1"},
    "yvelsurf":     {"long_name": "Surface y-velocity", "units": "m s-1"},
    "xvelbase":     {"long_name": "Basal x-velocity", "units": "m s-1"},
    "yvelbase":     {"long_name": "Basal y-velocity", "units": "m s-1"},
    "xvelmean":     {"long_name": "Vertical mean x-velocity", "units": "m s-1"},
    "yvelmean":     {"long_name": "Vertical mean y-velocity", "units": "m s-1"},
    "strbasemag":   {"long_name": "Basal drag", "units": "Pa"},
    "licalvf":      {"long_name": "Calving flux", "units": "kg m-2 s-1"},
    "sftgif":       {"long_name": "Land ice area fraction", "units": "1"},
    "sftgrf":       {"long_name": "Grounded ice area fraction", "units": "1"},
    "sftflf":       {"long_name": "Floating ice area fraction", "units": "1"},
}

ISMIP7_NX = 761
ISMIP7_NY = 761
ISMIP7_X0 = -3_040_000.0
ISMIP7_X1 = 3_040_000.0
ISMIP7_Y0 = -3_040_000.0
ISMIP7_Y1 = 3_040_000.0
ISMIP7_DX = 8000.0


def ismip7_grid_coords():
    r"""Return ISMIP7 AIS grid cell center coordinates."""
    x = np.linspace(ISMIP7_X0, ISMIP7_X1, ISMIP7_NX)
    y = np.linspace(ISMIP7_Y0, ISMIP7_Y1, ISMIP7_NY)
    return x, y


class ISMIP7Regridder:
    r"""Regrid Firedrake fields to the ISMIP7 AIS standard 8km grid."""

    def __init__(self, source_mesh):
        import firedrake as fd
        from firedrake.petsc import PETSc

        PETSc.Sys.Print("Building ISMIP7 regridder (VertexOnlyMesh)...")

        self._mesh = source_mesh
        self.x, self.y = ismip7_grid_coords()

        coords = source_mesh.coordinates.dat.data_ro
        xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
        ymin, ymax = coords[:, 1].min(), coords[:, 1].max()

        pad = ISMIP7_DX
        self._ii = np.where(
            (self.x >= xmin - pad) & (self.x <= xmax + pad)
        )[0]
        self._jj = np.where(
            (self.y >= ymin - pad) & (self.y <= ymax + pad)
        )[0]

        xx, yy = np.meshgrid(self.x[self._ii], self.y[self._jj])
        pts = np.column_stack([xx.ravel(), yy.ravel()])

        PETSc.Sys.Print(
            f"  Grid: {ISMIP7_NX}x{ISMIP7_NY}, "
            f"querying {len(self._ii)}x{len(self._jj)} = {len(pts)} points"
        )

        self._vom = fd.VertexOnlyMesh(
            source_mesh, pts, missing_points_behaviour="warn"
        )
        self._Q_vom = fd.FunctionSpace(self._vom, "DG", 0)
        self._Q_input = fd.FunctionSpace(self._vom.input_ordering, "DG", 0)

        n_found = self._vom.num_vertices()
        PETSc.Sys.Print(f"  Points inside mesh: {n_found} / {len(pts)}")

        self._n_query = len(pts)
        self._query_shape = (len(self._jj), len(self._ii))

    def regrid(self, field, fill_value=np.nan):
        r"""Regrid a scalar CG1 function to the ISMIP7 grid."""
        import firedrake as fd
        from firedrake.__future__ import interpolate

        f_vom = fd.assemble(interpolate(field, self._Q_vom))
        f_input = fd.Function(self._Q_input)
        f_input.dat.data_wo[:] = fill_value
        f_input.interpolate(f_vom)
        crop = np.array(f_input.dat.data_ro)

        grid = np.full((ISMIP7_NY, ISMIP7_NX), fill_value)
        grid[np.ix_(self._jj, self._ii)] = crop.reshape(self._query_shape)

        return grid

    def regrid_vector(self, field, fill_value=np.nan):
        r"""Regrid a vector CG1 function to the ISMIP7 grid."""
        import firedrake as fd

        data = field.dat.data_ro
        if data.ndim == 1:
            return self.regrid(field, fill_value=fill_value)

        dim = data.shape[1]
        Q_src = fd.FunctionSpace(self._mesh, "CG", 1)
        components = []
        for d in range(dim):
            f_comp = fd.Function(Q_src)
            f_comp.dat.data[:] = data[:, d]
            components.append(self.regrid(f_comp, fill_value=fill_value))

        return np.stack(components, axis=-1)

    def to_dataset(self, fields_dict, time=None, attrs=None):
        r"""Regrid fields and return as an xarray Dataset with ISMIP conventions."""
        import xarray as xr

        data_vars = {}
        for name, field in fields_dict.items():
            if isinstance(name, tuple):
                grid = self.regrid_vector(field)
                for d, comp_name in enumerate(name):
                    var_attrs = ISMIP_VARIABLES.get(comp_name, {})
                    data_vars[comp_name] = xr.Variable(
                        ["y", "x"], grid[:, :, d], attrs=var_attrs,
                    )
            else:
                data = field.dat.data_ro
                if data.ndim > 1 and data.shape[1] > 1:
                    grid = self.regrid_vector(field)
                    for d, suffix in enumerate(["x", "y"]):
                        comp_name = f"{suffix}vel{name}" if name in (
                            "surf", "base", "mean"
                        ) else f"{name}_{suffix}"
                        var_attrs = ISMIP_VARIABLES.get(comp_name, {})
                        data_vars[comp_name] = xr.Variable(
                            ["y", "x"], grid[:, :, d], attrs=var_attrs,
                        )
                else:
                    var_attrs = ISMIP_VARIABLES.get(name, {})
                    data_vars[name] = xr.Variable(
                        ["y", "x"], self.regrid(field), attrs=var_attrs,
                    )

        coords = {"x": self.x, "y": self.y}
        if time is not None:
            coords["time"] = time

        ds_attrs = {
            "Grid": "ISMIP7 AIS standard 8km",
            "proj4": "epsg:3031",
            "Conventions": "CF-1.10",
            "grid_resolution_m": ISMIP7_DX,
            "regridding_method": "Firedrake VertexOnlyMesh interpolation",
        }
        if attrs:
            ds_attrs.update(attrs)

        return xr.Dataset(data_vars, coords=coords, attrs=ds_attrs)

    def save_netcdf(self, fields_dict, filename, time=None, attrs=None):
        r"""Regrid fields and save to NetCDF."""
        ds = self.to_dataset(fields_dict, time=time, attrs=attrs)
        ds.to_netcdf(filename)
        ds.close()
