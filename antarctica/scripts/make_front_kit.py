#!/usr/bin/env python3
r"""Write the front-scoring subset of the ISMIP7 observations kit.

The kit is 9 GB; ``score_front.py`` reads five of its variables (the grid,
the Greene masks and their dates, the Mouginot basins), which compress to
9 MB because the masks are binary. The subset travels to a machine that
holds the forward runs but not the kit (``--kit`` on the scorer). Keep it
out of ``obs/mipkit/``: the forward's dH/dt reader takes any
``AntarcticaObsISMIP7*.nc`` there for the full kit, and the subset has no
dH/dt field.

    python antarctica/scripts/make_front_kit.py \
        ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc \
        ISMIP7/AIS/obs/front/AntarcticaObsISMIP7-v1.2-front.nc
"""
import sys

import netCDF4
import numpy as np

VARIABLES = ("x", "y", "greene_mask_time", "mouginot_basins", "icemask_greene")


def main(src, dst):
    s = netCDF4.Dataset(src)
    d = netCDF4.Dataset(dst, "w", format="NETCDF4")
    d.setncatts({k: s.getncattr(k) for k in s.ncattrs()})
    d.source_kit = src.split("/")[-1]
    d.note = "front-scoring subset of the ISMIP7 observations kit; see score_front.py"
    for dim in ("x", "y", "greene_mask_time"):
        d.createDimension(dim, len(s.dimensions[dim]))
    for name in VARIABLES:
        v = s[name]
        chunks = {"mouginot_basins": (1024, len(s.dimensions["x"])),
                  "icemask_greene": (1, 1024, len(s.dimensions["x"]))}.get(name)
        o = d.createVariable(name, v.dtype, v.dimensions, zlib=chunks is not None, complevel=4,
                             chunksizes=chunks)
        o.setncatts({k: v.getncattr(k) for k in v.ncattrs()})
        if name == "icemask_greene":
            for k in range(v.shape[0]):
                o[k, :, :] = v[k, :, :]
        else:
            o[:] = v[:]
    d.close()
    d = netCDF4.Dataset(dst)
    for name in VARIABLES:
        if name == "icemask_greene":
            for k in (0, s[name].shape[0] // 2, s[name].shape[0] - 1):
                assert np.array_equal(np.asarray(s[name][k]), np.asarray(d[name][k])), (name, k)
        else:
            assert np.array_equal(np.asarray(s[name][:]), np.asarray(d[name][:])), name
    print("wrote", dst)


if __name__ == "__main__":
    main(*sys.argv[1:3])
