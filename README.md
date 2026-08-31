# ISMIP7

Ice sheet simulations for [ISMIP7](https://www.ismip6.org/ismip7/) using
[icepack2](https://github.com/icepack/icepack2) with Firedrake and tlm_adjoint.

## Structure

```
ismip7/
├── icepack2_tools/         # Reusable utilities (mesh, forcing, regrid, eikonal, grounding zone)
├── antarctica/             # Antarctic continent simulations
│   ├── scripts/            # Pipeline scripts
│   ├── data/               # Downloaded datasets (see data/README.md)
│   ├── mesh/               # Generated meshes
│   ├── results/            # Inversion and sensitivity results
│   └── figs/               # Figures
```

## Dependencies

- [Firedrake](https://firedrakeproject.org)
- [icepack2](https://github.com/icepack/icepack2)
- [tlm_adjoint](https://github.com/jrmaddison/tlm_adjoint)
- gmsh, rasterio, icepack, geopandas, earthaccess
