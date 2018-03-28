import sys
import numpy as np
# GEOSPATIAL STUFF
# from osgeo import gdal, osr, ogr
from netCDF4 import Dataset
import json

def get_plotting_data(var_name, years, data_dir):
    doy_data_all = []
    doy_data_nz = []
    hist_data = np.array([])
    for yr_idx, year in enumerate(years):
        net_file = data_dir + var_name + '_5th_Indices_WUSA_' + str(year) + '.nc'
        print net_file
        try:
            ds = Dataset(net_file, 'r')
        except:
            continue
        doys = ds.variables['doy'][:]
        lats = ds.variables['lat'][:]
        lons = ds.variables['lon'][:]
        num_lats = lats.shape[0]
        num_lons = lons.shape[0]
        num_doys = doys.shape[0]
        if lats.size != 0  and lons.size != 0 and hist_data.size == 0:
            hist_data = np.empty([len(years), num_doys, lats.shape[0], lons.shape[0]])
        hist_data[yr_idx] = ds.variables['index'][:,:,:]

    # Ave over years
    hist_data = np.mean(hist_data, axis=0)
    for doy_idx in range(hist_data.shape[0]):
        dat = np.reshape(hist_data[doy_idx], num_lats * num_lons)
        dat_nz = dat[np.where(dat > 0)]
        mn = round(np.mean(dat), 4)
        mn_nz = round(np.mean(dat_nz), 4)
        doy_data_all.append(mn)
        doy_data_nz.append(mn_nz)

    # Ave over doys
    hist_data = np.mean(hist_data, axis = 0)
    hist_data_all = np.reshape(hist_data, num_lats * num_lons)
    hist_data_nz = hist_data[np.where(hist_data > 0)]
    hist_data_all = [round(v,4) for v in list(hist_data_all)]
    hist_data_nz = [round(v,4) for v in list(hist_data_nz)]
    return doy_data_all, doy_data_nz, hist_data_all, hist_data_nz
########
#M A I N
########
if __name__ == '__main__' :
    data_dir = 'RESULTS/livneh/'
    years = range(1951,2012)
    for var_name in ['tmin', 'tmax']:
        doy_data_all, doy_data_nz, hist_data_all, hist_data_nz = get_plotting_data(var_name, years, data_dir)
        print doy_data_all
        print doy_data_nz
        print len(hist_data_all)
        print len(hist_data_nz)
        with open(data_dir + var_name + '_ind_aves_doys.json', 'w') as outfile:
            json.dump(doy_data_all, outfile)
        with open(data_dir + var_name + '_ind_aves_doys_nz.json', 'w') as outfile:
            json.dump(doy_data_nz, outfile)
        with open(data_dir + var_name + '_hist_aves_all_doys_all_locs.json', 'w') as outfile:
            json.dump(hist_data_all, outfile)
        with open(data_dir + var_name + '_hist_aves_all_doys_all_locs_nz.json', 'w') as outfile:
            json.dump(hist_data_nz, outfile)
