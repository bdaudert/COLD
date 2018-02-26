from netCDF4 import Dataset
import os, glob
import numpy as np

def check_percentiles(netfile):
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    perc = ds.variables['percentile'][:,:,:]
    print ds.variables['percentile']
    print('FILE: ' + str(netfile))
    print('Num Lats ' + str(len(lats)))
    print('Num Lons ' + str(len(lons)))
    cnt = 0
    cnt_miss = 0
    cnt_cold = 0
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            if any(perc[:,lat_idx, lon_idx]) != -9999:
                # print indices[lat_idx][lon_idx]
                cnt += 1
            if any(perc[:,lat_idx, lon_idx]) == -9999:
                cnt_miss += 1
            if any(perc[:,lat_idx, lon_idx]) < -25:
                cnt_cold +=1
    print perc[0][20]
    print perc.shape
    print('Number of non-zero index arrays: ' + str(cnt))
    print('Number of missing index arrays: ' + str(cnt_miss))
    print('Number of brrrrrr....: ' + str(cnt_cold))
    # print indices[10][11]

def check_indices(netfile):
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    indices = ds.variables['index']
    print('FILE: ' + str(netfile))
    print('Num Lats ' + str(len(lats)))
    print('Num Lons ' + str(len(lons)))
    print('FILL VALUE INDICES: ')
    print np.where(indices == -9999)
    print('NAN : ')
    print np.where(np.isnan(indices))
    print("GREATER 0")
    print indices[np.where(indices <= 0)]
    #print indices[:,0:10,0:10]
    '''
    cnt = 0
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            if any(indices[lat_idx][lon_idx]) != 0:
                # print indices[lat_idx][lon_idx]
                cnt+=1
    print('Number of non-zero index arrays: ' + str(cnt))
    #print indices[10][11]
    '''
def find_bbox(netfile):
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    lat_min = 9999; lon_min = 9999
    lat_max = -9999; lon_max = -9999
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            if lat < lat_min:
                lat_min = lat
            if lat > lat_max:
                lat_max = lat
            if lon < lon_min:
                lon_min = lon
            if lon > lon_max:
                lon_max = lon
    return lat_min, lat_max, lon_min, lon_max

def count_lls(netfile):
    ds = Dataset(netfile, 'r')
    return len(ds.variables['lon'][:]), len(ds.variables['lat'][:])

def write_ll_file(netfile, outfile):
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]

    print 'Writing netcdf file'
    DS = Dataset(outfile, 'w', format='NETCDF3_64BIT')
    DS.description = '''LOCA latitudes and longitudes'''

    #Define the dimensions
    nlons = lons.shape[0]
    nlats = lats.shape[0]
    DS.createDimension('latitude', nlats)
    DS.createDimension('longitude', nlons)

    #Define the variables
    lat = DS.createVariable('lat', 'f4', ('latitude',), fill_value=1e20)
    lon = DS.createVariable('lon', 'f4', ('longitude',), fill_value=1e20)
    #Populate variable
    lat[:] = lats
    lon[:] = lons
    DS.close()


########
#M A I N
########
if __name__ == '__main__' :
    netfiles = filter(os.path.isfile, glob.glob('RESULTS/LIVNEH/' + '5th_Indices_WUSA_*.nc'))
    for netfile in netfiles[0:1]:
        check_indices(netfile)
    #netfile = 'RESULTS/LIVNEH/percentiles.nc'
    #check_percentiles(netfile)

    '''
    netfile = 'tasmin_day_CMCC-CM_historical_r1i1p1_19500101-19501231.LOCA_2016-04-02.16th.nc'
    ds = Dataset(netfile, 'r')
    print ds.variables['lon']
    # write_ll_file(netfile, 'LOCA_lls.nc')
    '''

    '''
    netfile = 'LOCA_lls.nc'
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    print np.min(lons)
    '''
