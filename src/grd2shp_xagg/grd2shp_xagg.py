"""Module to interpolate gridded climate forcings to geometry."""
import datetime
import pickle  # noqa: S403
import sys
import xagg as xa
from pathlib import Path
from xagg.export import prep_for_nc

import geopandas as gpd
import metpy.calc as mpcalc
import netCDF4
import numpy as np
import xarray as xr
from metpy.units import units

# import pandas as pd


# prsr = 101.3 * (((293.0-0.0065*Hru_elev_meters(i))/293.0)**5.26)
def std_pres(elev):
    """Calculate standard pressure.

    Args:
        elev (float): Elevation above sea-level in meters

    Returns:
        float: Standard Pressure (kPa)
    """
    return 101.325 * (((293.0 - 0.0065 * elev) / 293.0) ** 5.26)


def getaverage(data, wghts):
    """Get weighted avaerage.

    Args:
        data (float): array of values
        wghts (float): weights for each value in data

    Returns:
        float: weighted average
    """
    try:
        v_ave = np.average(data, weights=wghts)
    except ZeroDivisionError:
        v_ave = netCDF4.default_fillvals["f8"]
    return v_ave


def np_get_wval(ndata, wghts, hru_id=0, verbose=False):
    """Returns weighted average of ndata with weights.

    Args:
        ndata (float): The subset of values associated with the gridmet id's
                       that are mapped to hru_id
        wghts (float): Interpolation weights, user provided
        hru_id (int, optional): geometry id - can be used in debugging]. Defaults to 0.
        verbose(bool, optional): If True print geometry id of masked values

    Returns:
        [float]: [The weighted average of ndata based on weights wghts]
    """
    mdata = np.ma.masked_array(ndata, np.isnan(ndata))
    tmp = np.ma.average(mdata, weights=wghts)

    if tmp is np.ma.masked:
        if verbose:
            print(f"returning masked value: {hru_id}", ndata)
        return netCDF4.default_fillvals["f8"]

    else:
        return tmp


class Grd2ShpXagg:
    """Class to map or interpolate gridded data (netCDF, Grib).

    Data (focused on climate for now) onto geometry (Polygon) using
    area-weighted averages or rasterstats zonal statistics.
    """

    def __init__(self):
        """Init."""
        self.shp = None
        self.grd = None
        self.calctype = {0: "area-weighted average", 1: "zonal average"}
        self.type = None
        self.wieghts = None
        self.wght_id = None
        self.gdf = None
        self.gdf1 = None
        self.numgeom = 0
        self.numtimesteps = None
        self.unique_geom_ids = None
        self.time_var = None
        self.lat_var = None
        self.lon_var = None
        self.var = None
        self.var_output = None
        self._np_var = None
        self.str_start = None
        self.numvars = None
        self._start_date = None
        self._end_date = None
        self.current_time = None
        self.current_time_index = None
        self.mapped_vars = []

    def initialize(
        self,
        grd,
        shp,
        wght_file,
        time_var,
        lat_var,
        lon_var,
        var,
        var_output,
        opath,
        fileprefix,
        ctype,
    ):
        """Initialize.

        Args:
            grd ([type]): [description]
            shp ([type]): [description]
            wght_file ([type]): [description]
            time_var ([type]): [description]
            lat_var ([type]): [description]
            lon_var ([type]): [description]
            var ([type]): [description]
            var_output ([type]): [description]
            opath ([type]): [description]
            fileprefix ([type]): [description]
            ctype ([type]): [description]

        Raises:
            ValueError: [description]
            IOError: [description]
            IOError: [description]
        """
        self.valid_grd(grd)
        self.valid_shp(shp)

        self.shp = shp
        stringdict = {
            "time_var": time_var,
            "lat_var": lat_var,
            "lon_var": lon_var,
            "fileprefix": fileprefix,
        }

        for key, value in stringdict.items():
            if not isinstance(value, str):
                raise ValueError(f"arguement: {key}:{value} must be a string")
        self.time_var = time_var
        self.lat_var = lat_var
        self.lon_var = lon_var
        self.fileprefix = fileprefix

        self.valid_var(var=var, grd=grd)
        self.valid_var_ouput(var_output=var_output, var=var)

        self.opath = Path(opath)
        if self.opath.exists():
            print("output path exists", flush=True)
        else:
            sys.exit(f"Output Path does not exist: {self.opath} - EXITING")

        # assert(calctype in self.calctype)
        self.valid_calctype(ctype=ctype)

        try:
            with open(wght_file, "rb") as file:
                self.wieghts = pickle.load(file)  # noqa: S301
        except IOError as ie:
            raise IOError(f"Weight File error: {ie}")

        try:
            self.gdf = shp
            self.gdf.reset_index(drop=True, inplace=True)
        except IOError as ie:
            raise IOError(f"Geometry File error: {ie}")

        # grab some helpful vars. Assumption is dims are same for all vars!
        self.numvars = len(self.var)
        self.numtimesteps = self.grd[0].dims[self.time_var]
        self.str_start = np.datetime_as_string(self.grd[0][self.time_var][0], unit="D")
        self._start_date = self.grd[0][self.time_var][0]
        self._end_date = self.grd[0][self.time_var][self.numtimesteps - 1]

        self._np_var = np.zeros((self.numvars, self.numgeom, self.numtimesteps))
        print(f"numtimesteps: {self.numtimesteps} and Start date: {self.str_start}")

    def valid_grd(self, grd):
        """Test for valid grd.

        Args:
            grd (list): List of xarray grids

        Raises:
            ValueError: must be list of xarray.Dataset
            ValueError: members of list must be xarray.Dataset
        """
        if not isinstance(grd, list):
            raise ValueError(f"grd: {type(grd)} must be a list[xarray.Dataset]")
        else:
            for item in grd:
                if not isinstance(item, xr.Dataset):
                    raise ValueError(
                        f"Item {type(item)} - arg: grd must be an xarray.Dataset"
                    )
        self.grd = grd

    def valid_shp(self, shp):
        """Test for valid shp parameter.

        Must be a geopandas dataframe.

        Args:
            shp (Geopandas.GeoDataFrame): Must be geopandas.GeoDataFrame

        Raises:
            ValueError: Returns ValueError if not
        """
        if not isinstance(shp, gpd.GeoDataFrame):
            raise ValueError(f"shp: {shp} must be a Geopandas Datafram")

    def valid_var(self, var, grd):
        """Test for valid var.

        Args:
            var (list): List of strings
            grd (list): List of xarray.Dataset

        Raises:
            ValueError: var must be list
            ValueError: var must have same len as grd
            ValueError: var must be list of strings
        """
        if not isinstance(var, list):
            raise ValueError(
                f"Arguement var:{var} must be a list of strings.  Enclose in []"
            )
        if not len(grd) == len(var):
            raise ValueError(
                f"Length of grd: {len(grd)} must be equal to length of var: {len(var)}"
            )
        for item in var:
            if not isinstance(item, str):
                raise ValueError(f"Item in var: {item} must be a string")
        self.var = var

    def valid_var_ouput(self, var_output, var):
        """Test for valid var_ouput.

        Args:
            var_output (list): List of names variable names used in output
            var (list): List of input variable names

        Raises:
            ValueError: var_ouput must be a list
            ValueError: len(var_ouput) must equal len(var)
            ValueError: var_output members must be strings
        """
        if not isinstance(var_output, list):
            raise ValueError(
                f"Arguement var:{var_output} must be a list of strings.  Enclose in []"
            )
        if not len(var_output) == len(var):
            raise ValueError(
                f"Length of var_output: {len(var_output)} "
                f"must be equal to length of var: {len(var)}"
            )
        for item in var_output:
            if not isinstance(item, str):
                raise ValueError(f"Item in var: {item} must be a string")
        self.var_output = var_output

    def valid_calctype(self, ctype):
        """Test for valid calctype.

        Args:
            ctype (int): One of values in {0: "area-weighted average", 1: "zonal average"}

        Raises:
            ValueError: Not in {0: "area-weighted average", 1: "zonal average"}
        """
        if ctype not in self.calctype:
            raise ValueError(
                f"calctype {ctype} must be one of Calculation Types {self.calctype}"
            )
        self.type = ctype

    def run_weights(self):
        """Run weights."""
        for index, tvar in enumerate(self.var):
            print(f"generating mapped vales for {tvar} ...")
            grid = self.grd[index]
            aggragated = xa.aggregate(grid, self.wieghts)
            xr_agg = prep_for_nc(aggragated)
            self.mapped_vars.append(xr_agg)
            print(f"finished mapped values for {tvar}")

    @property
    def start_date(self):
        """Return Start Date.

        Returns:
            datetime: Start date of gridded data.
        """
        return self._start_date

    @property
    def end_date(self):
        """Return End Date.

        Returns:
            datetime: End date of gridded data.
        """
        return self._end_date

    @property
    def current_date(self):
        """Return Current Date.

        Returns:
            datetime: Current date of gridded data.
        """
        return self.current_date

    @property
    def num_timesteps(self):
        """Return number of time-steps.

        Returns:
            int: Number of time-steps
        """
        return self.numtimesteps

    def mapped_data(self, var):
        """Return xarray of mapped data by var.

        Args:
            var (string): Should be value in self.var

        Returns:
            xarray: Mapped values
        """
        print(var)
        try:
            index = self.var_output.index(var)
        except ValueError as ve:
            print(f"val not in {self.var}", ve)
            return

        return self.mapped_vars[index]

    def write_file(  # noqa: C901
        self, elev_file, punits=0, datetag=None, filename=None, append=False
    ):
        """Write netcdf file of resulting interpolation.

        Args:
            elev_file ([type]): [description]
            punits (int, optional): [description]. Defaults to 0.
            datetag ([type], optional): [description]. Defaults to None.
            filename ([type], optional): [description]. Defaults to None.
            append (bool, optional): [description]. Defaults to False.

        """
        if datetag is None:
            datetag = str(datetime.datetime.now().strftime("%Y_%m_%d"))

        if not append:
            ncfile = self.create_ncf(datetag)
        else:
            ncfile = netCDF4.Dataset(
                self.opath
                / (
                    self.fileprefix
                    + "climate_"
                    + str(datetime.datetime.now().strftime("%Y%m%d"))
                    + ".nc"
                ),
                mode="a",
                format="NETCDF_CLASSIC",
            )

        for index, tvar in enumerate(self.var_output):
            vartype = self.grd[index][self.var[index]].dtype
            ncvar = ncfile.createVariable(tvar, vartype, ("geomid", "time"))
            ncvar.fill_value = netCDF4.default_fillvals["f8"]
            ncvar.long_name = self.grd[index][self.var[index]].long_name
            ncvar.standard_name = self.grd[index][self.var[index]].standard_name
            ncvar.description = self.grd[index][self.var[index]].description
            # ncvar.grid_mapping = 'crs'
            ncvar.units = self.grd[index][self.var[index]].units
            if tvar in ["tmax", "tmin"]:
                if punits == 1:
                    conv = units.degC
                    ncvar[:, :] = (
                        units.Quantity(
                            self._np_var[index, 0 : self.current_time_index, :],
                            ncvar.units,
                        )
                        .to(conv)
                        .magnitude
                    )
                    ncvar.units = conv.format_babel()
                else:
                    conv = units.degF
                    ncvar[:, :] = (
                        units.Quantity(
                            self._np_var[index, 0 : self.current_time_index, :],
                            ncvar.units,
                        )
                        .to(conv)
                        .magnitude
                    )
                    ncvar.units = conv.format_babel()
            elif tvar == "prcp":
                if punits == 1:
                    conv = units("mm")
                    ncvar[:, :] = (
                        units.Quantity(
                            self._np_var[index, 0 : self.current_time_index, :],
                            ncvar.units,
                        )
                        .to(conv)
                        .magnitude
                    )
                    ncvar.units = conv.units.format_babel()
                else:
                    conv = units("inch")
                    ncvar[:, :] = (
                        units.Quantity(
                            self._np_var[index, 0 : self.current_time_index, :],
                            ncvar.units,
                        )
                        .to(conv)
                        .magnitude
                    )
                    ncvar.units = conv.units.format_babel()
            else:
                ncvar[:, :] = self._np_var[index, 0 : self.current_time_index, :]
                ncvar.units = self.grd[index][self.var[index]].units

        elevf = gpd.read_file(elev_file, layer="hru_elev")
        elev = elevf["hru_elev"].values

        if all(x in self.var_output for x in ["tmax", "tmin", "shum"]):
            tmax_ind = self.var_output.index("tmax")
            tmin_ind = self.var_output.index("tmin")
            shum_ind = self.var_output.index("shum")

        print(f"tmaxind: {tmax_ind}, tminind: {tmin_ind}, shumind: {shum_ind}")

        rel_h = np.zeros((self.current_time_index, self.numgeom))
        for j in np.arange(np.int(self.numgeom)):
            pr = mpcalc.height_to_pressure_std(units.Quantity(elev[j], "m"))
            for i in np.arange(np.int(self.current_time_index)):
                tmax = units.Quantity(self._np_var[tmax_ind, i, j], units.kelvin)
                tmin = units.Quantity(self._np_var[tmin_ind, i, j], units.kelvin)
                spch = units.Quantity(self._np_var[shum_ind, i, j], "kg/kg")
                rhmax = mpcalc.relative_humidity_from_specific_humidity(pr, tmax, spch)
                rhmin = mpcalc.relative_humidity_from_specific_humidity(pr, tmin, spch)
                rel_h[i, j] = (rhmin.magnitude + rhmax.magnitude) / 2.0

        ncvar = ncfile.createVariable("humidity", rel_h.dtype, ("geomid", "time"))
        ncvar.units = "1"
        ncvar.fill_value = netCDF4.default_fillvals["f8"]
        ncvar[:, :] = rel_h[0 : self.current_time_index, :]

        ncfile.close()

    def create_ncf(self, datetag):
        """Create netCDF file.

        Args:
            datetag ([type]): [description]

        Returns:
            [type]: [description]
        """
        ncfile = netCDF4.Dataset(
            self.opath
            / (self.fileprefix + "climate_" + datetag.strftime("%Y_%m_%d") + ".nc"),
            mode="w",
            format="NETCDF4_CLASSIC",
        )

        def getxy(pt):
            return pt.x, pt.y

        centroidseries = self.gdf1.geometry.centroid.to_crs(epsg=4327)
        tlon, tlat = [list(t) for t in zip(*map(getxy, centroidseries))]

        # Global Attributes
        ncfile.Conventions = "CF-1.8"
        ncfile.featureType = "timeSeries"
        ncfile.history = ""

        sp_dim = len(self.gdf1.index)
        # Create dimensions

        ncfile.createDimension("geomid", size=sp_dim)  # hru_id
        ncfile.createDimension(
            "time", size=None
        )  # unlimited axis (can be appended to).

        # Create Variables
        time = ncfile.createVariable("time", "f4", ("time",))
        time.long_name = "time"
        time.standard_name = "time"
        time.units = "days since " + self.str_start
        time.calendar = "standard"
        time[:] = np.arange(0, self.current_time_index, dtype=np.float)

        hru = ncfile.createVariable("geomid", "i", ("geomid",))
        hru.cf_role = "timeseries_id"
        hru.long_name = "local model hru id"
        hru[:] = np.asarray(self.gdf1.index)

        lat = ncfile.createVariable("hru_lat", np.dtype(np.float32).char, ("geomid",))
        lat.long_name = "Latitude of HRU centroid"
        lat.units = "degrees_north"
        lat.standard_name = "hru_latitude"
        lat[:] = tlat

        lon = ncfile.createVariable("hru_lon", np.dtype(np.float32).char, ("geomid",))
        lon.long_name = "Longitude of HRU centroid"
        lon.units = "degrees_east"
        lon.standard_name = "hru_longitude"
        lon[:] = tlon

        return ncfile
