import sys
import numpy as np
import json
#GEOSPATIAL STUFF
# from osgeo import gdal, osr, ogr
from netCDF4 import Dataset

def get_plotting_data(years, data_dir):
    year_sums_all = []
    year_sums_nz = []
    hist_data_all = []
    hist_data_nz = []
    # fort box plots
    ll_data = {}
    count = 0
    for year_idx, year in enumerate(years):
        net_file = data_dir + '5th_Indices_WUSA_' + str(year) + '.nc'
        print net_file
        try:
            ds = Dataset(net_file, 'r')
        except:
            year_sums.append(None)
            continue
        num_lons = ds.variables['lon'][:].shape[0]
        num_lats = ds.variables['lat'][:].shape[0]
        indices = ds.variables['index'][:]
        ll_sums = 0
        ll_sums_non_zero = 0
        indices_all = indices[np.where(indices != -9999)]
        indices_nz = indices[np.where(indices >= 1)]
        sum_all = np.sum(indices_all)
        sum_nz = np.sum(indices_nz)
        hist_data_all.append(sum_all)
        hist_data_nz.append(sum_nz)
        year_sums_all.append(sum_all/(num_lats*num_lons))
        year_sums_nz.append(sum_nz/(num_lats*num_lons))

        '''
        for lat_idx, lat in enumerate(lats):
            for lon_idx, lon in enumerate(lons):
                all_indices = [i for i in indices[lat_idx][lon_idx] if i != -9999]
                non_zero_indices = [i for i in indices[lat_idx][lon_idx] if i >= 1]
                s = sum(all_indices)
                s_non_zero = sum(non_zero_indices)
                if non_zero_indices:
                    hist_data_nz.append(s_non_zero)
                    count+=1
                ll_str = str(lon) + ',' + str(lat)
                if count <= 5 and len(ll_data.keys()) < 5 and year_idx == 0 and non_zero_indices:
                    ll_str = str(lon) + ',' + str(lat)
                    ll_data[ll_str] = [s_non_zero]
                if ll_str in ll_data.keys() and year_idx != 0:
                    ll_data[ll_str].append(s_non_zero)
                ll_sums += s
                ll_sums_non_zero += s_non_zero
                hist_data_all.append(s)
        year_sums_all.append(ll_sums / (len(lats) * len(lons)))
        year_sums_nz.append(ll_sums_non_zero / count);
        '''
    # return year_sums_all, year_sums_nz, hist_data_all, hist_data_nz, ll_data
    return year_sums_all, year_sums_nz, hist_data_all, hist_data_nz
########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951,2012)
    data_dir = 'RESULTS/'
    # year_sums, year_sums_non_zero, hist_data, hist_data_non_zero, box_plots_5locs = get_plotting_data(years, data_dir)
    year_sums, year_sums_non_zero, hist_data, hist_data_non_zero = get_plotting_data(years, data_dir)
    '''
    with open(data_dir + 'box_plots_sums_5locs.json', 'w') as outfile:
        json.dump(box_plots_5locs, outfile)
    '''
    with open(data_dir + 'ind_sums_years.json', 'w') as outfile:
        json.dump(year_sums, outfile)
    with open(data_dir + 'ind_sums_years_non_zero.json', 'w') as outfile:
        json.dump(year_sums_non_zero, outfile)
    with open(data_dir + 'hist_sum_all_years_and_locs.json', 'w') as outfile:
        json.dump(hist_data, outfile)
    with open(data_dir + 'hist_sum_all_years_and_locs_non_zero.json', 'w') as outfile:
        json.dump(hist_data_non_zero, outfile)
    print year_sums
    print year_sums_non_zero

