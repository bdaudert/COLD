import sys
import numpy as np
# GEOSPATIAL STUFF
# from osgeo import gdal, osr, ogr
from netCDF4 import Dataset
import json

def get_lls(year, data_dir):
    f_name = data_dir + '5th_Indices_WUSA_' + str(year) +'.nc'
    ds = Dataset(f_name, 'r')
    lats = ds.variables['lat']
    lons = ds.variables['lon']
    return lats, lons

def get_plotting_data(years, data_dir):
    doy_data_all = []
    doy_data_nz = []
    hist_data_all = []
    hist_data_nz = []
    year_data_all = []
    year_data_nz = []
    for yr_idx, year in enumerate(years):
        net_file = data_dir + '5th_Indices_WUSA_' + str(year) + '.nc'
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
        data_all = ds.variables['index'][:,:,:]
        for doy_idx in range(num_doys):
            doy_vals_all = np.reshape(data_all[doy_idx,:,:],(num_lats*num_lons,))
            doy_vals_nz = doy_vals_all[np.where(doy_vals_all > 0)]
            year_data_all.append(doy_vals_all)
            year_data_nz.append(doy_vals_nz)
            hist_data_all.append(round(np.mean(doy_vals_all), 4))
            if doy_vals_nz.shape[0] > 0:
                hist_data_nz.append(round(np.mean(doy_vals_nz), 4))
            else:
                hist_data_nz.append(-9999)
    for doy_idx in range(num_doys):
        mn = round(np.mean(year_data_all[doy_idx]), 4)
        doy_data_all.append(mn)
        mn = round(np.mean(year_data_nz[doy_idx]), 4)
        doy_data_nz.append(mn)
    return doy_data_all, doy_data_nz, hist_data_all, hist_data_nz
########
#M A I N
########
if __name__ == '__main__' :
    data_dir = 'RESULTS/LIVNEH/'
    years = range(1951,2012)
    doy_data_all, doy_data_nz, hist_data_all, hist_data_nz = get_plotting_data(years, data_dir)

    with open(data_dir + 'ind_aves_doys.json', 'w') as outfile:
        json.dump(doy_data_all, outfile)
    with open(data_dir + 'ind_aves_doys_non_zero.json', 'w') as outfile:
        json.dump(doy_data_nz, outfile)
    with open(data_dir + 'hist_aves_all_doys_all_locs.json', 'w') as outfile:
        json.dump(hist_data_all, outfile)
    with open(data_dir + 'hist_aves_all_doys_all_locs_non_zero.json', 'w') as outfile:
        json.dump(hist_data_nz, outfile)
