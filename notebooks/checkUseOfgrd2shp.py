"""A file."""
# %%
import pickle  # noqa: S403

import geopandas as gpd
import xarray as xr

import grd2shp_xagg


# %%
# %matplotlib inline


ds2 = xr.open_dataset("../data/tmax_full.nc")

with open("../data/weight_gfv1_1.txt", "rb") as file:
    agg = pickle.load(file)  # noqa: S301

gdf = gpd.read_file("../data/GFv1.1_v2e_geographic.shp")

# %%
gm_vars = ["air_temperature"]

g2s = grd2shp_xagg.Grd2ShpXagg()

g2s.initialize(
    grd=[ds2],
    shp=gdf,
    wght_file="../data/weight_gfv1_1.txt",
    time_var="day",
    lat_var="lat",
    lon_var="lon",
    var=gm_vars,
    var_output=["tmax"],
    ctype=0,
)
# %%
g2s.run_weights()
# %%
val = g2s.mapped_data("tmax")

# %%
g2s.write_gmcfsv2_file(
    opath="../data", prefix="test_", elev_file="../data/package.gpkg"
)
# %%
val
# %%
val2 = val.isel(time=0)
val2
# %%
val.to_netcdf("../data/tmax_out_test.nc")
# %%
