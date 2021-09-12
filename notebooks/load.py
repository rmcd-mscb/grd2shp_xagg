"""Load data."""
# %%
import gridmet_cfsv2 as gm

# import grd2shp_xagg
# import geopandas as gpd
# import numpy as np
# import time
# import xagg as xa
# import pickle

gm_vars = [
    "air_temperature",
    "air_temperature",
    "precipitation_amount",
    "specific_humidity",
]
m = gm.Gridmet(type=3)

# %%
ds = m.tmin
ds.load()
day_dim = ds.dims["day"]
ds.air_temperature[0:15, 0, :, :] = ds.air_temperature[32:47, 0, :, :]
ds.air_temperature[0:15, 1, :, :] = ds.air_temperature[16:31, 1, :, :]
ds.air_temperature[16:31, 0, :, :] = ds.air_temperature[32:47, 0, :, :]

ds.air_temperature[16:31, day_dim - 1, :, :] = ds.air_temperature[
    0:15, day_dim - 1, :, :
]
ds.air_temperature[32:47, day_dim - 1, :, :] = ds.air_temperature[
    0:15, day_dim - 1, :, :
]
ds.air_temperature[32:47, day_dim - 2, :, :] = ds.air_temperature[
    16:31, day_dim - 2, :, :
]
ds.to_netcdf("../data/tmin_full.nc")
# %%
