import sys
import numpy as np
import json
#GEOSPATIAL STUFF
from netCDF4 import Dataset

def get_plotting_data(years, data_dir):
    box_data = {}
    # fort box plots
    for year_idx, year in enumerate(years):
        net_file = data_dir + '5th_Indices_WUSA_' + str(year) + '.nc'
        print net_file
        try:
            ds = Dataset(net_file, 'r')
        except:
            continue
        num_lons = ds.variables['lon'][:].shape[0]
        num_lats = ds.variables['lat'][:].shape[0]
        num_doys = ds.variables['doy'][:].shape[0]
        indices = ds.variables['index'][:,:,:]
        ll_sums = 0
        b_data = []
        for doy_idx in range(num_doys):
            doy_vals_all = np.reshape(indices[doy_idx,:,:],(num_lats*num_lons,))
            doy_vals_nz = doy_vals_all[np.where(doy_vals_all > 0)]
            b_data.append(doy_vals_nz.shape[0])
        '''
        for lat_idx, lat in enumerate(lats):
            for lon_idx, lon in enumerate(lons):
                num_days = len([i for i in indices[lat_idx][lon_idx] if i > 0])
                b_data.append(num_days)
        '''
        box_data[year] = b_data
    return box_data
########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951,1979)
    data_dir = 'RESULTS/'
    box_data = get_plotting_data(years, data_dir)
    with open(data_dir + 'box_plots_num_days_all_years_and_locs.json', 'w') as outfile:
        json.dump(box_data, outfile)

