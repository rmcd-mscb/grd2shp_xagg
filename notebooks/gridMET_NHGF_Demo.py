"""Module to test gridmet workflow."""
# %%
import pickle  # noqa S403
import time

import geopandas as gpd
import xagg as xa
import xarray as xr

import grd2shp_xagg

# import numpy as np

# %%
tmax_url = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmx_1979_CurrentYear_CONUS.nc"  # noqa B950
tmin_url = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmn_1979_CurrentYear_CONUS.nc"  # noqa B950
prcp_url = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_pr_1979_CurrentYear_CONUS.nc"  # noqa B950

# %%
ds_prcp = xr.open_dataset(prcp_url + "#fillmismatch")
ds_tmin = xr.open_dataset(tmin_url + "#fillmismatch")
ds_tmax = xr.open_dataset(tmax_url + "#fillmismatch")
# %%
del_gdf = gpd.read_file("/mnt/c/Users/rmcd/git/grd2shp_xagg/data/nhru_01.shp")
ds_wght = xr.open_dataset(
    "/mnt/c/Users/rmcd/git/grd2shp_xagg/data/agg_met_tmmx_1979_CurrentYear_CONUS.nc"
)
# %%
del_gdf
# %%
# uncomment below if first time through notebook to generate weights

start = time.perf_counter()
weightmap = xa.pixel_overlaps(ds_wght, del_gdf)
end = time.perf_counter()
print(f"finished agg in {round(end-start, 2)} second(s)")

# save weights for future use
with open("../data/nhru_01_weights.pickle", "wb") as file:
    pickle.dump(weightmap, file)

# %%
# setup grd2shp_xagg

start_date = "1979-01-01"
end_date = "1979-01-07"
lon_min = -74
lon_max = -67
lat_min = 40
lat_max = 49

# subset gridMET data
ds_prcp_1 = ds_prcp.sel(
    day=slice(start_date, end_date),
    lon=slice(lon_min, lon_max),
    lat=slice(lat_max, lat_min),
)
ds_tmin_1 = ds_tmin.sel(
    day=slice(start_date, end_date),
    lon=slice(lon_min, lon_max),
    lat=slice(lat_max, lat_min),
)
ds_tmax_1 = ds_tmax.sel(
    day=slice(start_date, end_date),
    lon=slice(lon_min, lon_max),
    lat=slice(lat_max, lat_min),
)

ds_tmax_1.load()
ds_tmin_1.load()
ds_prcp_1.load()

g2s = grd2shp_xagg.Grd2ShpXagg()
g2s.initialize(
    grd=[ds_tmax_1, ds_tmin_1, ds_prcp_1],
    shp=del_gdf,
    wght_file="../data/nhru_01_weights.pickle",
    time_var="day",
    lat_var="lat",
    lon_var="lon",
    var=[
        "daily_maximum_temperature",
        "daily_minimum_temperature",
        "precipitation_amount",
    ],
    var_output=["tmax", "tmin", "prcp"],
    ctype=0,
)

# %%
g2s.run_weights()
# %%
val = g2s.mapped_data("prcp")
val
# %%
g2s.write_gm_file(opath="../data", prefix="test1_")
# %%
ds_tmax_1
# %%
