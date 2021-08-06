"""A file."""
import gridmet_cfsv2 as gm
import pickle  # noqa: S403
import time
import xagg as xa
from xagg.export import prep_for_nc

import geopandas as gpd
import numpy as np
import xarray as xr

# import grd2shp_xagg


# %matplotlib inline


gm_vars = [
    "air_temperature",
    "air_temperature",
    "precipitation_amount",
    "specific_humidity",
]
m = gm.Gridmet(type=3)
# %%

""" To save time, this block of code exports the resulting netcdf file so that it
doesn't need to be re-downloaded. This block of code loads tmax into memory so the
matrix of ensemble (48 ensembles) vs day (30 days for each ensemble,
 but the ensembles are established from three seperate days of simulation
so there is a need to fill in missing days) """

# ds = m.tmax
# ds.load()
# day_dim = ds.dims['day']
# ds.air_temperature[0:15, 0, :, :] = ds.air_temperature[32:47, 0, :, :]
# ds.air_temperature[0:15, 1, :, :] = ds.air_temperature[16:31, 1, :, :]
# ds.air_temperature[16:31, 0, :, :] = ds.air_temperature[32:47, 0, :, :]

# ds.air_temperature[16:31, day_dim-1, :, :] = ds.air_temperature[0:15, day_dim-1, :, :]
# ds.air_temperature[32:47, day_dim-1, :, :] = ds.air_temperature[0:15, day_dim-1, :, :]
# ds.air_temperature[32:47, day_dim-2, :, :] = ds.air_temperature[16:31, day_dim-2, :, :]
# ds.to_netcdf('../data/tmax_full.nc')

# %%
ds2 = xr.open_dataset("../data/tmax_full.nc")
ds2.air_temperature.isel(time=0, day=0).plot()
# %%
with open("../data/weight_gfv1_1.txt", "rb") as file:
    agg = pickle.load(file)  # noqa: S301
# %%
start = time.perf_counter()
aggragated = xa.aggregate(ds2, agg)
end = time.perf_counter()
print(f"finished agg in {round(end-start, 2)} second(s)")
# %%

xrd = prep_for_nc(aggragated)
xrd
# %%
aggragated.to_netcdf("../data/tmax_full_mapped_out.nc")
# %%
"""the following code failes because it expects aggragated to be a 2D array and
   instead it's ndim eqauls three because there are 48 ensembles in this data set
   so as below if we want to save to a dataframe we sould neet to select one of the
   ensembles"""
mapped = aggragated.to_dataframe()
mapped

# %%

# %%
start = time.perf_counter()
aggragated = xa.aggregate(ds2.isel(time=0), agg)
end = time.perf_counter()
print(f"finished agg in {round(end-start, 2)} second(s)")
# %%
mapped = aggragated.to_dataframe()
mapped
# %%
gdf = gpd.read_file("../data/GFv1.1_v2e_geographic.shp")
gdf["temp1"] = mapped.air_temperature0.values
gdf.plot(column="temp1")
# %%
print(type(aggragated.agg.air_temperature))
print(aggragated.agg.air_temperature.size)
print(type(aggragated.agg.air_temperature.values))
print(type(aggragated.agg.air_temperature.values[0]))
print(type(aggragated.agg.air_temperature.values[0][0]))
print(type(aggragated.agg.air_temperature.values[0][0][0]))

# %%
aggragated.agg.air_temperature[0]

# %%

# %%
dvar = 1
dgeom = 139807
dday = 32
test = np.zeros((dvar, dgeom, dday))
start = time.perf_counter()
for index, value in aggragated.agg.air_temperature.items():
    test[0, index, :] = value[0]
end = time.perf_counter()
print(f"finished array in {round(end-start, 2)} second(s)")
