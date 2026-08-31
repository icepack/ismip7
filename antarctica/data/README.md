# Data

Observational inputs downloaded by `python scripts/download_data.py`
(BedMachine, MEaSUREs velocity, RACMO SMB). See §1 of
[`antarctica/README.md`](../README.md) for the dataset table, accounts,
and where each product lands.

`scripts/gl_sensitivity.py` additionally reads the MEaSUREs ice-shelf
polygons ([NSIDC-0709](https://nsidc.org/data/nsidc-0709), not fetched by
`download_data.py`) from the path in its `SHAPEFILE` constant.
