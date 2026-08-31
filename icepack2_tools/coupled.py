r"""Coupled ice sheet -- subshelf plume model for ISMIP7 simulations.

Couples the icepack2 dual-form SSA (diagnostic velocity + DG0 thickness
transport) with the plume depth-integrated subshelf ocean model to compute
physically consistent melt rates driven by cavity circulation.

The coupling is asynchronous: within each ice sheet time step (dt_ice ~ 1 yr),
the plume model is subcycled at its own timestep (dt_plume ~ 10-30 min) to
reach a quasi-steady melt distribution for the current cavity geometry. The
melt rate field is then passed to the ice thickness evolution equation.

Coupling interface:
    Ice -> Plume:  ice_draft, bed_elevation, floating_mask
    Plume -> Ice:  melt_rate (m/yr ice equivalent)
    Ocean -> Plume: ambient T, S at ice draft depth (from ISMIP7 ocean data)
"""

import numpy as np
from firedrake import Constant, Function, max_value, assemble, dx


def compute_ice_draft(h, s, rho_I=917.0, rho_W=1024.0):
    r"""Compute ice shelf draft from thickness and surface elevation."""
    return s - h


def compute_floating_mask(h, b, s, rho_I=917.0, rho_W=1024.0):
    r"""Return 1 where ice is floating (haf <= 0), 0 where grounded."""
    s_float = b + (rho_W / rho_I) * max_value(-b, Constant(0.0))
    haf = s - s_float
    from firedrake import conditional, le
    return conditional(le(haf, 0.0), 1.0, 0.0)


class CoupledModel:
    r"""Coupled ice sheet - plume model.

    Wraps the ice sheet simulation context and a PlumeModel instance,
    handling the data transfer and subcycling between them.
    """

    def __init__(self, ice_ctx, plume_model, ocean_forcing=None,
                 plume_dt=600.0, plume_spinup_time=86400.0):
        """
        Parameters
        ----------
        ice_ctx : dict
            Ice sheet model context from simulation.setup_model().
        plume_model : plume.PlumeModel
            Initialized plume model on the same mesh or a submesh.
        ocean_forcing : ISMIP7Ocean, optional
            Ocean forcing reader for ambient T/S.
        plume_dt : float
            Plume model timestep in seconds (default 600 = 10 min).
        plume_spinup_time : float
            Total plume integration time per ice step in seconds
            (default 86400 = 1 day). Should be long enough for the
            plume to reach quasi-steady state.
        """
        self.ice_ctx = ice_ctx
        self.plume = plume_model
        self.ocean = ocean_forcing
        self.plume_dt = plume_dt
        self.plume_spinup_time = plume_spinup_time
        self._sec_per_year = 31556926.0

    def update_plume_geometry(self):
        r"""Transfer current ice geometry to the plume model."""
        h = self.ice_ctx["h"]
        s = self.ice_ctx["s"]
        draft = compute_ice_draft(h, s)
        self.plume.ice_draft.interpolate(draft)

        if self.plume.bed_elevation is not None:
            b = self.ice_ctx["b"]
            self.plume.bed_elevation.interpolate(b)

        if self.plume.wet_mask is not None:
            h_val = self.ice_ctx["h"]
            b_val = self.ice_ctx["b"]
            s_val = self.ice_ctx["s"]
            floating = compute_floating_mask(h_val, b_val, s_val)
            self.plume.wet_mask.interpolate(floating)

    def update_ambient_from_ocean(self, year):
        r"""Update plume ambient T/S from ISMIP7 ocean forcing data."""
        if self.ocean is None:
            return

        plume_mesh = self.plume.ice_draft.function_space().mesh()
        mesh_x = plume_mesh.coordinates.dat.data_ro[:, 0]
        mesh_y = plume_mesh.coordinates.dat.data_ro[:, 1]
        draft = self.plume.ice_draft.dat.data_ro

        # For the plume model we need absolute T and S, not just TF.
        # Use the ISMIP7 thetao and so fields directly if available, via
        # the same cached in-memory interpolators as get_thermal_forcing.
        # Fallback: approximate from TF assuming T_f ~ -1.9 C at surface.
        if hasattr(self.ocean, "_lookup"):
            t_amb = self.ocean._lookup(
                "thetao", year, mesh_x, mesh_y, draft=draft, nan_fill=-1.9
            )
            s_amb = self.ocean._lookup(
                "so", year, mesh_x, mesh_y, draft=draft, nan_fill=34.5
            )
            if t_amb is not None and s_amb is not None:
                self.plume.ambient_temperature.dat.data[:] = t_amb
                self.plume.ambient_salinity.dat.data[:] = s_amb
                return

        # Fallback: estimate T from TF
        tf = self.ocean.get_thermal_forcing(year, mesh_x, mesh_y, draft=draft)
        self.plume.ambient_temperature.dat.data[:] = tf - 1.9
        self.plume.ambient_salinity.dat.data[:] = 34.5

    def subcycle_plume(self):
        r"""Run the plume model to quasi-steady state."""
        n_steps = int(self.plume_spinup_time / self.plume_dt)
        for _ in range(n_steps):
            self.plume.step(self.plume_dt)

    def get_melt_rate(self):
        r"""Get the plume melt rate in m/yr ice equivalent on the ice mesh."""
        return self.plume.melt_rate_ice_equivalent()

    def forcing_callback(self, ctx, t_yr):
        r"""Forcing callback compatible with simulation.run_simulation().

        Called each ice time step. Updates plume geometry, ambient
        conditions, subcycles the plume, and sets the ocean_melt field.
        """
        self.update_plume_geometry()
        self.update_ambient_from_ocean(t_yr)
        self.subcycle_plume()

        melt_ice = self.get_melt_rate()
        # Convert from m/s to m/yr
        melt_myr = Function(ctx["Q"])
        melt_myr.dat.data[:] = melt_ice.dat.data_ro * self._sec_per_year
        ctx["ocean_melt"].dat.data[:] = melt_myr.dat.data_ro

        # Also apply atmosphere forcing if available
        if hasattr(self, 'atm') and self.atm is not None:
            from icepack2_tools.forcing import smb_kgm2s_to_myr
            mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
            mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]
            smb = self.atm.get_smb(t_yr, mesh_x, mesh_y, anomaly=True)
            ctx["accum"].dat.data[:] = smb
