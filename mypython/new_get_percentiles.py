#!/usr/local/pythonenv/csc/bin/python
import numpy as np
from netCDF4 import Dataset
import json

def read_vars(var_name, year):
    f_name = 'DATA/' + var_name + '.' + str(year) + '.nc'
    ds = Dataset(f_name, 'r')
    print ds.variables


def get_lls(var_name, year):
    f_name = 'DATA/' + var_name + '.' + str(year) + '.nc'
    ds = Dataset(f_name, 'r')
    lats = ds.variables['lat'][:]
    lons = ds.variables['lon'][:]
    return lats, lons

def write_netcdf(PCTLS, lons, lats, doys, out_file):
    print 'Writing netcdf file'
    DS =  Dataset(out_file, 'w', format='NETCDF3_64BIT')
    DS.description = '''
        At each livneh gridpoint in the Western United States
        (BBOX: : -125, 31, -102, 49.1) the 5th percentile of tmin
        is computed. The base period used was 1951 - 2005
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
    pctl = DS.createVariable('percentile', 'i4', ('day_in_season', 'latitude', 'longitude'), fill_value=-9999)
    pctl.units = 'Deg Celsius'

    #Populate variable
    lat[:] = lats
    lon[:] = lons
    doy[:] = doys
    pctl[:,:,:] = PCTLS
    DS.close()

def compute_percentile(a):
    return int(round(np.percentile(a, 5)))

def get_percentiles(num_days, start_doy, latbounds, lonbounds, var_name, years, data_dir, out_file):
    year_change = False
    if start_doy + num_days > 365:
        end_doy = num_days - (365 - start_doy)
        doys = np.concatenate((np.arange(start_doy,365), np.arange(0,end_doy)))
        year_change = True
    else:
        end_doy = start_doy + num_days
        doys = np.arange(start_doy, end_doy)

    all_lats, all_lons = get_lls('tmin', years[0])
    # latitude lower and upper index
    latli = np.argmin( np.abs( all_lats - latbounds[0] ) )
    latui = np.argmin( np.abs( all_lats - latbounds[1] ) )
    # longitude lower and upper index
    lonli = np.argmin( np.abs( all_lons - lonbounds[0] ) )
    lonui = np.argmin( np.abs( all_lons - lonbounds[1] ) )
    fh = open(out_file, 'w+')
    fh.close()
    lats = np.array([])
    lons = np.array([])
    DOY_DATA = np.array([])
    for year_idx, year in enumerate(years):
        print('PROCESSING YEAR ' + str(year))
        p_year = year - 1
        c_year = year
        # Get this years_data
        f_name = 'DATA/' + var_name + '.' + str(c_year) + '.nc'
        if year_change:
            end_doy_temp = 365
        else:
            end_doy_temp = end_doy
        try:
            this_year_data = Dataset(f_name, 'r').variables[var_name][start_doy:end_doy_temp, latli:latui, lonli:lonui]
            if lats.size == 0:
                lats = Dataset(f_name, 'r').variables['lat'][latli:latui]
            if lons.size == 0:
                lons = Dataset(f_name, 'r').variables['lon'][lonli:lonui]
        except:
            # Last year reached
            break
        # Initialize DOY_DATA
        if lats.size != 0  and lons.size != 0 and DOY_DATA.size == 0:
            DOY_DATA = np.empty([len(years), num_days, lats.shape[0], lons.shape[0]])

        if year_change:
            # get data from previous year
            f_name = data_dir + var_name + '.' + str(p_year) + '.nc'
            try:
                # FIX ME: How to deal with fill values
                last_year_data = Dataset(f_name, 'r').variables[var_name][0:end_doy, latli:latui, lonli:lonui]
                #last_year_data = list(Dataset(f_name, 'r').variables[var_name][333:365,lat_idx,lon_idx])
                #last_year_data = [float(v) for v in last_year_data]
                # get the valid lats
            except:
                # On to next year
                continue
            year_data = np.concatenate((last_year_data, this_year_data), axis=0)
        else:
            year_data = this_year_data

        del this_year_data, last_year_data
        DOY_DATA[year_idx] = year_data
        del year_data

    print('GETTING PERCENTILES')
    PCTLS = np.empty([num_days, lats.shape[0], lons.shape[0]])
    for doy_idx in range(num_days):
        PCTLS[doy_idx] = np.apply_along_axis(compute_percentile, 0, DOY_DATA[:,doy_idx,:,:])
    del DOY_DATA
    print PCTLS
    write_netcdf(PCTLS, lons, lats, doys, out_file)

########
#M A I N
########
if __name__ == '__main__' :
    # read_vars('tmin', 2011)
    out_file = 'TEST_percentiles.nc'
    data_dir = 'DATA/'
    out_dir = 'RESULTS/'
    years = range(1951, 1953)
    num_days = 90
    start_doy = 334
    var_name = 'tmin'
    latbounds = [39 , 39.2 ]
    lonbounds = [ -119 + 360 , -118.8 + 360] # degrees east ?
    get_percentiles(num_days, start_doy, latbounds, lonbounds, var_name, years, data_dir, out_dir + out_file)
    # check_percentiles(out_file)

