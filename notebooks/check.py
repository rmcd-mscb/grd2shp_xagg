"""check."""
import gridmet_cfsv2 as gm
import time

import geopandas as gpd

import grd2shp_xagg

# import pickle  # noqa: S403
# import xarray as xr

# import xagg


def fill_ds(ds, var):
    """Fill missing days in ensembles."""
    day_dim = ds.dims["day"]
    ds[var][0:15, 0, :, :] = ds[var][32:47, 0, :, :]
    ds[var][0:15, 1, :, :] = ds[var][16:31, 1, :, :]
    ds[var][16:31, 0, :, :] = ds[var][32:47, 0, :, :]

    ds[var][16:31, day_dim - 1, :, :] = ds[var][0:15, day_dim - 1, :, :]
    ds[var][32:47, day_dim - 1, :, :] = ds[var][0:15, day_dim - 1, :, :]
    ds[var][32:47, day_dim - 2, :, :] = ds[var][16:31, day_dim - 2, :, :]


if __name__ == "__main__":

    gm_vars = [
        "air_temperature",
        "air_temperature",
        "precipitation_amount",
        "specific_humidity",
    ]

    # %matplotlib inline
    m = gm.Gridmet(type=3)

    ds1 = m.tmax
    ds1.load()
    fill_ds(ds1, gm_vars[0])
    ds1.to_netcdf("../data/tmax_full3.nc")

    ds2 = m.tmin
    ds2.load()
    fill_ds(ds2, gm_vars[1])
    ds2.to_netcdf("../data/tmin_full3.nc")

    ds3 = m.prcp
    ds3.load()
    fill_ds(ds3, gm_vars[2])
    ds3.to_netcdf("../data/prcp_full3.nc")

    ds4 = m.specific_humidity
    ds4.load()
    fill_ds(ds4, gm_vars[3])
    ds4.to_netcdf("../data/shum_full3.nc")

    # ds2 = xr.open_dataset("./data/tmax_full3.nc")
    # ds1 = xr.open_dataset("./data/tmin_full3.nc")
    # gdf = gpd.read_file("./data/nhru_01.shp")
    gdf = gpd.read_file("./data/GFv1.1_v2e_geographic.shp")
    # weightmap = xagg.pixel_overlaps(ds2, gdf)
    # with open("../data/weight_nhru1.txt", 'wb') as file:
    #     pickle.dump(weightmap, file)
    # with open("./data/weight_nhru1.txt", "rb") as file:
    #     agg = pickle.load(file)  # noqa: S301

    # gm_vars = ["air_temperature", "air_temperature"]
    # g2s = grd2shp_xagg.Grd2ShpXagg()
    # g2s.initialize(
    #     grd=[ds2],
    #     shp=gdf,
    #     wght_file="./data/weight_nhru1.txt",
    #     ens_var="time",
    #     time_var="day",
    #     lat_var="lat",
    #     lon_var="lon",
    #     var=gm_vars,
    #     var_output=["tmax"],
    #     ctype=0,
    #     threaded=False
    # )

    # g2s.run_weights()
    # start = time.perf_counter()
    # g2s.write_gmcfsv2_file(opath='./data', prefix='test_',
    #                        elev_file='./data/package.gpkg')
    # end = time.perf_counter()
    # print(f'write in {round(end-start, 2)} second(s)')

    g2sa = grd2shp_xagg.Grd2ShpXagg()
    g2sa.initialize(
        grd=[ds1, ds2, ds3, ds4],
        shp=gdf,
        # wght_file="./data/weight_nhru1.txt",
        wght_file="./data/weight_gfv1_1.txt",
        ens_var="time",
        time_var="day",
        lat_var="lat",
        lon_var="lon",
        var=gm_vars,
        var_output=["tmin", "tmax", "prcp", "shum"],
        ctype=0,
        threaded=True,
    )
    start = time.perf_counter()
    g2sa.run_weights_pp()
    end = time.perf_counter()
    print(f"agg in {round(end-start, 2)} second(s)")

    start = time.perf_counter()
    g2sa.write_gmcfsv2_file(
        opath="./data", prefix="testa_", elev_file="./data/package.gpkg"
    )
    end = time.perf_counter()
    print(f"write in {round(end-start, 2)} second(s)")
