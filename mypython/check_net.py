from netCDF4 import Dataset
import os, glob

def check_percentiles(netfile):
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    perc = ds.variables['percentile'][:,:,:]
    print('FILE: ' + str(netfile))
    print('Num Lats ' + str(len(lats)))
    print('Num Lons ' + str(len(lons)))
    cnt = 0
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            if any(perc[lat_idx][lon_idx]) != 1e20:
                # print indices[lat_idx][lon_idx]
                cnt += 1
    print perc
    print perc.shape
    print('Number of non-zero index arrays: ' + str(cnt))
    # print indices[10][11]

def check_indices(netfile):
    ds = Dataset(netfile, 'r')
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    indices = ds.variables['index'][:,:,:]
    print('FILE: ' + str(netfile))
    print('Num Lats ' + str(len(lats)))
    print('Num Lons ' + str(len(lons)))
    cnt = 0
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            if any(indices[lat_idx][lon_idx]) > 0:
                # print indices[lat_idx][lon_idx]
                cnt+=1
    print('Number of non-zero index arrays: ' + str(cnt))
    #print indices[10][11]

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
########
#M A I N
########
if __name__ == '__main__' :
    '''
    netfiles = filter(os.path.isfile, glob.glob('DATA/INDICES_CA/' + '5th_Indices_WUSA_*.nc'))
    for netfile in netfiles[0:1]:
        # check_results(netfile)
        lat_min, lat_max, lon_min, lon_max = find_bbox(netfile)
        print(lat_min, lat_max, lon_min, lon_max)
    '''
    path = 'RESULTS/'
    netfile = '5th_Indices_WUSA_1951.nc'
    check_indices(path + netfile)
