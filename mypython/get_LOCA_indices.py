import numpy as np
from netCDF4 import Dataset
from netcdftime import utime
import json

def get_lls():
    f_name = 'LOCA_lls.nc'
    ds = Dataset(f_name, 'r')
    lats = ds.variables['lat'][:]
    lons = ds.variables['lon'][:]
    return lats, lons

def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx

def is_float_convertible(val):
    try:
        np.isnan(val)
        return False
    except:
        pass

    try:
        float(val)
        return True
    except:
        return False

def write_netcdf(INDICES, var_name, model, rcp, year, lons, lats, doys, out_dir):
    print 'Writing netcdf file'
    out_file = out_dir + var_name + '_' + model + '_' + rcp + '_5th_Indices_WUSA_' + str(year) + '.nc'
    DS = Dataset(out_file, 'w', format='NETCDF3_64BIT')
    DS.description = '''
           At each livneh gridpoint in the Western United States
           (BBOX: : -125, 31, -102, 49.1) for each DOY in Winter season(DJF),
           the number degrees below the 5th percentile at that point are recorded
           '''

    #Define the dimensions
    nlons = lons.shape[0]
    nlats = lats.shape[0]
    ndays = doys.shape[0]
    DS.createDimension('latitude', nlats)
    DS.createDimension('longitude', nlons)
    DS.createDimension('day_in_season', ndays)

    #Define the variables
    lat = DS.createVariable('lat', 'f4', ('latitude',), fill_value=1e20)
    lon = DS.createVariable('lon', 'f4', ('longitude',), fill_value=1e20)
    # FIX ME: save doys as ints
    # OverflowError: Python int too large to convert to C long
    doy = DS.createVariable('doy', 'f4', ('day_in_season',), fill_value=1e20)
    ind = DS.createVariable('index', 'i4', ('day_in_season', 'latitude', 'longitude'), fill_value=1e20)
    ind.units = 'DegC below 5th precentile'

    #Populate variable
    lat[:] = lats
    lon[:] = lons
    doy[:] = doys
    ind[:,:,:] = INDICES
    DS.close()

def compute_indices(year_data, perc_data):
    y_data = year_data
    np.divide(y_data[y_data == 1e20]= np.nan, 1000.0)
    diff = np.absolute(np.rint(np.subtract(perc_data,year_data)))
    ind =  np.where(np.less_equal(y_data, perc_data), diff, 0)
    return ind


def get_indices(years, var_name, perc_file, data_dir, out_dir):
    '''
    For each year in years and day in season,  we compute the cold index
    at each lon, lat: the number of degrees below the 5th percentile
    :param years: list of years (ints)
    :param var_name: variable name
    :param perc_file: file containing the percentiles for each lon, lat and day of season
    :param out_dir: output directory
    :return: num years files containing the index data at each lon, lat and day of season
             out_dir + '5th_Indices_WUSA_' + str(year)+ '.nc',
    '''
    if var_name == 'tmin':
        loca_var_name = 'tasmin'
    if var_name == 'tmax':
        loca_var_name = 'tasmax'

    # Read percentile file
    PD = Dataset(perc_file, 'r')
    lons = PD.variables['lon'][:]
    lats = PD.variables['lat'][:]
    lat_min = lats[np.argmin(lats)]
    lat_max = lats[np.argmax(lats)]
    latbounds = [lat_min, lat_max]
    lon_min = lons[np.argmin(lons)]
    lon_max = lons[np.argmax(lons)]
    lonbounds = [lon_min, lon_max]
    # FIX ME:in  get_percentiles needed--- does not save the doys as ints!!!
    doys = np.array([int(d) for d in list(PD.variables['doy'][:])])
    num_days = doys.shape[0]
    start_doy = doys[0]
    perc_data = PD.variables['percentile'][:,:,:]
    PD.close()

    year_change = False
    if start_doy + num_days > 365:
        end_doy = num_days - (365 - start_doy)
        year_change = True
    else:
        end_doy = start_doy + num_days

    # Read data file to find all lats, lons and set the bounds
    # according to lats, lons in percentile file
    all_lats, all_lons = get_lls()
    # latitude lower and upper index
    latli = np.argmin(np.abs(all_lats - latbounds[0]))
    latui = np.argmin(np.abs(all_lats - latbounds[1]))
    # longitude lower and upper index
    lonli = np.argmin(np.abs(all_lons - lonbounds[0]))
    lonui = np.argmin(np.abs(all_lons - lonbounds[1]))

    # print latbounds, lonbounds
    # print all_lats[latli], all_lats[latui], all_lons[lonli], all_lons[lonui]
    del all_lats, all_lons

    for year_idx, year in enumerate(years):
        print('PROCESSING YEAR ' + str(year))
        p_year = year - 1
        c_year = year
        # Get this year's data
        f_name = data_dir + str(c_year) + '.h5'
        if year_change:
            end_doy_temp = 365
        else:
            end_doy_temp = end_doy
        try:
            this_year_data = np.array(h5py.File(f_name, 'r')[loca_var_name])
            this_year_data = this_year_data[0:end_doy, latli:latui, lonli:lonui]
        except:
            # Last year reached
            break
        # Get last year's data if year change
        if year_change:
            # get December data from previous year
            f_name = data_dir + str(p_year) + '.h5'
            try:
                # FIX ME: How to deal with fill values
                last_year_data = np.array(h5py.File(f_name, 'r')[loca_var_name])
                last_year_data = last_year_data[start_doy:365, latli:latui, lonli:lonui]
            except:
                # On to next year
                continue
            year_data = np.concatenate((last_year_data, this_year_data), axis=0)
        else:
            year_data = this_year_data
        # del this_year_data, last_year_data
        print('GETTING INDICES')
        INDICES = np.empty([num_days, lats.shape[0], lons.shape[0]], dtype=int)
        for doy_idx in range(num_days):
            INDICES[doy_idx] =  compute_indices(year_data[doy_idx], perc_data[doy_idx])

        write_netcdf(INDICES, var_name, model, rcp, year, lons, lats, doys, out_dir)
########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951, 2012)
    # years = range(1951, 1953)
    data_dir = '/media/DataSets/loca/'
    in_file_name = 'percentiles.nc'
    out_dir = 'RESULTS/loca/'
    perc_file = out_dir + in_file_name
    for model in LOCA_CMIP5_MODELS.keys():
        for rcp in rcps:
            data_dir = loca_dir + model + '/' + rcp + '/'
            for var_name in ['tmin', 'tmax']:
                perc_file = out_dir + var_name + '_' +  model + '_' + rcp + '_percentiles.nc'
                out_file = var_name + '_' + model + '_' + rcp + '_percentiles.nc'
                get_indices(years, var_name, perc_file, data_dir, out_dir)
