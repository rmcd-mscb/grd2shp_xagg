"""check."""
import pickle  # noqa: S403
import time

import geopandas as gpd
import xarray as xr

import grd2shp_xagg

# import xagg

if __name__ == "__main__":

    # %matplotlib inline
    ds2 = xr.open_dataset("./data/tmax_full.nc")
    gdf = gpd.read_file("./data/nhru_01.shp")
    # weightmap = xagg.pixel_overlaps(ds2, gdf)
    # with open("../data/weight_nhru1.txt", 'wb') as file:
    #     pickle.dump(weightmap, file)
    with open("./data/weight_nhru1.txt", "rb") as file:
        agg = pickle.load(file)  # noqa: S301

    gm_vars = ["air_temperature"]
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
        grd=[ds2],
        shp=gdf,
        wght_file="./data/weight_nhru1.txt",
        ens_var="time",
        time_var="day",
        lat_var="lat",
        lon_var="lon",
        var=gm_vars,
        var_output=["tmax"],
        ctype=0,
        threaded=True,
    )

    g2sa.run_weights()

    start = time.perf_counter()
    g2sa.write_gmcfsv2_file(
        opath="./data", prefix="testa_", elev_file="./data/package.gpkg"
    )
    end = time.perf_counter()
    print(f"write in {round(end-start, 2)} second(s)")
