import numpy as np
from netCDF4 import Dataset
from netcdftime import utime
import json

def get_lls(var_name, year):
    f_name = 'DATA/' + var_name + '.' + str(year) + '.nc'
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

def write_netcdf(INDICES, year, lons, lats, doys, out_dir):
    print 'Writing netcdf file'
    DS = Dataset(out_dir + '5th_Indices_WUSA_' + str(year) + '.nc', 'w', format='NETCDF3_64BIT')
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
    ind = DS.createVariable('index', 'i4', ('day_in_season', 'latitude', 'longitude'), fill_value=-9999)
    ind.units = 'DegC below 5th precentile'

    #Populate variable
    lat[:] = lats
    lon[:] = lons
    doy[:] = doys
    ind[:,:,:] = INDICES
    DS.close()

def compute_indices(year_data, perc_data):
    y_data = year_data
    y_data[y_data == 1e20] = np.nan
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
    all_lats, all_lons = get_lls('tmin', years[0])
    # latitude lower and upper index
    latli = np.argmin(np.abs(all_lats - latbounds[0]))
    latui = np.argmin(np.abs(all_lats - latbounds[1]))
    # longitude lower and upper index
    lonli = np.argmin(np.abs(all_lons - lonbounds[0]))
    lonui = np.argmin(np.abs(all_lons - lonbounds[1]))

    # print latbounds, lonbounds
    # print all_lats[latli], all_lats[latui], all_lons[lonli], all_lons[lonui]
    del all_lats, all_lons

    # PCA_data: [[doy1_year1_indices for eachll], [doy1_year1_indices for eachll], ....[doy90_year_last_indices for each ll]]
    # PCA_data = [[] for i in range(len(years) * num_days)]
    PCA_data = []
    for year_idx, year in enumerate(years):
        print('PROCESSING YEAR ' + str(year))
        p_year = year - 1
        c_year = year
        # Get this year's data
        f_name = data_dir + var_name + '.' + str(c_year) + '.nc'
        if year_change:
            end_doy_temp = 365
        else:
            end_doy_temp = end_doy
        try:
            this_year_data = Dataset(f_name, 'r').variables[var_name][start_doy:end_doy_temp, latli:latui+1, lonli:lonui+1]
        except:
            # Last year reached
            break
        # Get last year's data if year change
        if year_change:
            # get December data from previous year
            f_name = data_dir + var_name + '.' + str(p_year) + '.nc'
            try:
                # FIX ME: How to deal with fill values
                last_year_data = Dataset(f_name, 'r').variables[var_name][0:end_doy, latli:latui+1, lonli:lonui+1]
                #last_year_data = list(Dataset(f_name, 'r').variables[var_name][333:365,lat_idx,lon_idx])
                #last_year_data = [float(v) for v in last_year_data]
                # get the valid lats
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

        write_netcdf(INDICES, year, lons, lats, doys, out_dir)
########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951, 2012)
    # years = range(1951, 1953)
    var_name = 'tmin'
    data_dir = 'DATA/'
    in_file_name = 'percentiles.nc'
    out_dir = 'RESULTS/LIVNEH/'
    perc_file = out_dir + in_file_name
    get_indices(years, var_name, perc_file, data_dir, out_dir)
