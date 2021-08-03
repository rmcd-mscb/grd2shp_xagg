import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
import sys
import netCDF4
import datetime
import metpy.calc as mpcalc
from metpy.units import units


# prsr = 101.3 * (((293.0-0.0065*Hru_elev_meters(i))/293.0)**5.26)
def std_pres(elev):
    return 101.325 * (((293.0-0.0065*elev)/293.0)**5.26)


def getaverage(data, wghts):
    try:
        v_ave = np.average(data, weights=wghts)
    except ZeroDivisionError:
        v_ave = netCDF4.default_fillvals['f8']
    return v_ave


def np_get_wval(ndata, wghts, hru_id=0, verbose=False):
    """
    [Returns weighted average of ndata with weights]

    Args:
        ndata ([float]): [The subset of values associated with the gridmet id's that are mapped to hru_id]
        wghts ([float]): [Interpolation weights, user provided]
        hru_id (int, optional): [geometry id - can be used in debugging]. Defaults to 0.
        verbose(bool, optional): [If True print geometry id of masked values]

    Returns:
        [float]: [The weighted average of ndata based on weights wghts]
    """
    mdata = np.ma.masked_array(ndata, np.isnan(ndata))
    tmp = np.ma.average(mdata, weights=wghts)

    if tmp is np.ma.masked:
        if verbose:
            print(f'returning masked value: {hru_id}', ndata)
        return netCDF4.default_fillvals['f8']

    else:
        return tmp


class Grd2Shp:
    """
    Class to map or interpolate gridded data (netCDF, Grib) data (focused on climate for now) onto
    geometry (Polygon) using area-weighted averages or rasterstats zonal statistics
    """

    def __init__(self):
        self.shp = None
        self.grd = None
        self.calctype = {
            0: "area-weighted average",
            1: "zonal average"
            }
        self.type = None
        self.wght_file = None
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

    def initialize(self,
                   grd,
                   shp,
                   wght_file,
                   time_var: str,
                   lat_var: str,
                   lon_var: str,
                   var,
                   var_output,
                   opath: str,
                   fileprefix: str = '',
                   calctype: int = 0,
                   ):

        if not isinstance(grd, list):
            raise ValueError(f'grd: {type(grd)} must be a list[xarray.Dataset]')
        else:
            for item in grd:
                if not isinstance(item, xr.Dataset):
                    raise ValueError(f'Item {type(item)} - arg: grd must be an xarray.Dataset')
        self.grd = grd

        if not isinstance(shp, gpd.GeoDataFrame):
            raise ValueError(
                f'shp: {shp} must be a Geopandas Datafram')
        self.shp = shp
        stringdict = {'time_var': time_var, 'lat_var': lat_var,
                      'lon_var': lon_var, 'fileprefix': fileprefix}
        for key, value in stringdict.items():
            if not isinstance(value, str):
                raise ValueError(f'arguement: {key}:{value} must be a string')
        self.time_var = time_var
        self.lat_var = lat_var
        self.lon_var = lon_var
        self.fileprefix = fileprefix

        if not isinstance(var, list):
            raise ValueError(f'Arguement var:{var} must be a list of strings.  Enclose in []')
        if not len(grd) == len(var):
            raise ValueError(f'Length of grd: {len(grd)} must be equal to length of var: {len(var)}')
        for item in var:
            if not isinstance(item, str):
                raise ValueError(f'Item in var: {item} must be a string')
        self.var = var

        if not isinstance(var_output, list):
            raise ValueError(f'Arguement var:{var_output} must be a list of strings.  Enclose in []')
        if not len(var_output) == len(var):
            raise ValueError(f'Length of var_output: {len(var_output)} must be equal to length of var: {len(var)}')
        for item in var_output:
            if not isinstance(item, str):
                raise ValueError(f'Item in var: {item} must be a string')
        self.var_output = var_output

        self.opath = Path(opath)
        if self.opath.exists():
            print('output path exists', flush=True)
        else:
            sys.exit(f'Output Path does not exist: {self.opath} - EXITING')

        try:
            calctype in self.calctype
        except ValueError as ve:
            print(f'calctype {calctype} must be one of Calculation Types {self.calctype}')
            raise ve
        self.type = calctype

        try:
            self.wght_file = pd.read_csv(wght_file)
        except IOError as ie:
            raise IOError(f'Weight File error: {ie}')

        try:
            self.gdf = shp
            self.gdf.reset_index(drop=True, inplace=True)
        except IOError as ie:
            raise IOError(f'Geometry File error: {ie}')

        # grab the geom_id from the weights file and use as identifier below
        self.wght_id = self.wght_file.columns[1]

        # this geodataframe merges all hru-ids and dissolves so the length of the index
        # equals the number of hrus
        self.gdf1 = self.gdf.sort_values(self.wght_id).dissolve(by=self.wght_id)
        self.numgeom = len(self.gdf1.index)
        self.geomindex = np.asarray(self.gdf1.index, dtype=np.int)

        # group by the weights_id for processing
        self.unique_geom_ids = self.wght_file.groupby(self.wght_id)

        # grab some helpful vars. Assumption is dims are same for all vars!
        self.numvars = len(self.var)
        self.numtimesteps = self.grd[0].dims[self.time_var]
        self.str_start = np.datetime_as_string(self.grd[0][self.time_var][0], unit='D')
        self._start_date = self.grd[0][self.time_var][0]
        self._end_date = self.grd[0][self.time_var][self.numtimesteps-1]
        self.current_time_index = 0
        self.current_time = self.grd[0][self.time_var][self.current_time_index]

        self._np_var = np.zeros((self.numvars, self.numtimesteps, self.numgeom))
        print(f'numtimesteps: {self.numtimesteps} and Start date: {self.str_start}')

    def run_weights(self):

        for index, tvar in enumerate(self.var):
            grid = self.grd[index]

            timestep = self.current_time_index
            print(f'Processing timestep: {timestep}', flush=True)

            val_interp = np.zeros(self.numgeom)
            val_flat_interp = grid[tvar].values[timestep, :, :].flatten(order='K')

            for i in np.arange(len(self.geomindex)):
                try:
                    weight_id_rows = self.unique_geom_ids.get_group(self.geomindex[i])
                    tw = weight_id_rows.w.values
                    tgid = weight_id_rows.grid_ids.values
                    tmp = getaverage(val_flat_interp[tgid], tw)
                    if np.isnan(tmp):
                        val_interp[i] = np_get_wval(val_flat_interp[tgid], tw, self.geomindex[i])
                    else:
                        val_interp[i] = tmp
                except KeyError:
                    val_interp[i] = netCDF4.default_fillvals['f8']

                if i % 10000 == 0:
                    print(f'    Processing {tvar} for hru {i}', flush=True)

            self._np_var[index, timestep, :] = val_interp[:]

        self.current_time_index += 1
        # self.current_time = self.grd[0][self.time_var][self.current_time_index]

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

    @property
    def current_date(self):
        return self.current_date

    @property
    def num_timesteps(self):
        return self.numtimesteps

    @property
    def current_mapped_data(self):
        return self._np_var[:, self.current_time_index, :]

    def write_file(self, elev_file, punits=0, datetag=None, filename=None, append=False):
        if datetag is None:
            datetag = str(datetime.datetime.now().strftime('%Y_%m_%d'))

        if not append:
            ncfile = netCDF4.Dataset(
                self.opath / (self.fileprefix + 'climate_' + datetag.strftime("%Y_%m_%d")
                              + '.nc'),
                mode='w', format='NETCDF4_CLASSIC'
            )

            def getxy(pt):
                return pt.x, pt.y

            centroidseries = self.gdf1.geometry.centroid.to_crs(epsg=4327)
            tlon, tlat = [list(t) for t in zip(*map(getxy, centroidseries))]

            # Global Attributes
            ncfile.Conventions = 'CF-1.8'
            ncfile.featureType = 'timeSeries'
            ncfile.history = ''

            sp_dim = len(self.gdf1.index)
            # Create dimensions

            ncfile.createDimension('geomid', size=sp_dim)  # hru_id
            ncfile.createDimension('time', size=None)  # unlimited axis (can be appended to).

            # Create Variables
            time = ncfile.createVariable('time', 'f4', ('time',))
            time.long_name = 'time'
            time.standard_name = 'time'
            time.units = 'days since ' + self.str_start
            time.calendar = 'standard'
            time[:] = np.arange(0, self.current_time_index, dtype=np.float)

            hru = ncfile.createVariable('geomid', 'i', ('geomid',))
            hru.cf_role = 'timeseries_id'
            hru.long_name = 'local model hru id'
            hru[:] = np.asarray(self.gdf1.index)

            lat = ncfile.createVariable('hru_lat', np.dtype(np.float32).char, ('geomid',))
            lat.long_name = 'Latitude of HRU centroid'
            lat.units = 'degrees_north'
            lat.standard_name = 'hru_latitude'
            lat[:] = tlat

            lon = ncfile.createVariable('hru_lon', np.dtype(np.float32).char, ('geomid',))
            lon.long_name = 'Longitude of HRU centroid'
            lon.units = 'degrees_east'
            lon.standard_name = 'hru_longitude'
            lon[:] = tlon

            # crs = ncfile.createVariable('crs', np.dtype(np.int))
            # crs.GeoTransform = self.grd[0].crs.GeoTransform
            # # crs.NAME = self.grd[0].crs.NAME
            # crs.grid_mapping_name = self.grd[0].crs.grid_mapping_name
            # crs.inverse_flattening = self.grd[0].crs.inverse_flattening
            # crs.long_name = self.grd[0].crs.long_name
            # crs.longitude_of_prime_meridian = self.grd[0].crs.longitude_of_prime_meridian
            # crs.semi_major_axis = self.grd[0].crs.semi_major_axis
            # crs.spatial_ref = self.grd[0].crs.spatial_ref

        else:
            ncfile = netCDF4.Dataset(
                self.opath / (self.fileprefix + 'climate_' + str(datetime.datetime.now().strftime('%Y%m%d'))
                              + '.nc'),
                mode='a', format='NETCDF_CLASSIC'
            )

        for index, tvar in enumerate(self.var_output):
            vartype = self.grd[index][self.var[index]].dtype
            ncvar = ncfile.createVariable(tvar, vartype, ('time', 'geomid'))
            ncvar.fill_value = netCDF4.default_fillvals['f8']
            ncvar.long_name = self.grd[index][self.var[index]].long_name
            ncvar.standard_name = self.grd[index][self.var[index]].standard_name
            ncvar.description = self.grd[index][self.var[index]].description
            # ncvar.grid_mapping = 'crs'
            ncvar.units = self.grd[index][self.var[index]].units
            if tvar in ['tmax', 'tmin']:
                if punits == 1:
                    conv = units.degC
                    ncvar[:, :] = units.Quantity(self._np_var[index, 0:self.current_time_index, :], ncvar.units)\
                        .to(conv).magnitude
                    ncvar.units = conv.format_babel()
                else:
                    conv = units.degF
                    # ncvar[:,:] = ((self._np_var[index, 0:self.current_time_index, :]-273.15)*1.8)+32.0
                    ncvar[:, :] = units.Quantity(self._np_var[index, 0:self.current_time_index, :], ncvar.units)\
                        .to(conv).magnitude
                    ncvar.units = conv.format_babel()
            elif tvar == 'prcp':
                if punits == 1:
                    conv = units('mm')
                    ncvar[:, :] = units.Quantity(self._np_var[index, 0:self.current_time_index, :], ncvar.units)\
                        .to(conv).magnitude
                    ncvar.units = conv.units.format_babel()
                else:
                    conv = units('inch')
                    ncvar[:, :] = units.Quantity(self._np_var[index, 0:self.current_time_index, :], ncvar.units)\
                        .to(conv).magnitude
                    ncvar.units = conv.units.format_babel()
                # else units are already  in mm
                # ncvar[:,:] = np.multiply(self._np_var[index, 0:self.current_time_index, :], conv.magnitude)
            else:
                ncvar[:, :] = self._np_var[index, 0:self.current_time_index, :]
                ncvar.units = self.grd[index][self.var[index]].units

        elevf = gpd.read_file(elev_file, layer='hru_elev')
        elev = elevf['hru_elev'].values

        if all(x in self.var_output for x in ['tmax', 'tmin', 'shum']):
            tmax_ind = self.var_output.index('tmax')
            tmin_ind = self.var_output.index('tmin')
            shum_ind = self.var_output.index('shum')

        print(f'tmaxind: {tmax_ind}, tminind: {tmin_ind}, shumind: {shum_ind}')

        rel_h = np.zeros((self.current_time_index, self.numgeom))
        for j in np.arange(np.int(self.numgeom)):
            pr = mpcalc.height_to_pressure_std(units.Quantity(elev[j], "m"))
            for i in np.arange(np.int(self.current_time_index)):
                tmax = units.Quantity(self._np_var[tmax_ind, i, j], units.kelvin)
                tmin = units.Quantity(self._np_var[tmin_ind, i, j], units.kelvin)
                spch = units.Quantity(self._np_var[shum_ind, i, j], "kg/kg")
                rhmax = mpcalc.relative_humidity_from_specific_humidity(pr, tmax, spch)
                rhmin = mpcalc.relative_humidity_from_specific_humidity(pr, tmin, spch)
                rel_h[i, j] = (rhmin.magnitude + rhmax.magnitude)/2.0

        ncvar = ncfile.createVariable('humidity', rel_h.dtype, ('time', 'geomid'))
        ncvar.units = "1"
        ncvar.fill_value = netCDF4.default_fillvals['f8']
        ncvar[:, :] = rel_h[0:self.current_time_index, :]

        ncfile.close()
