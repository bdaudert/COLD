import sys
import numpy as np
#GEOSPATIAL STUFF
from osgeo import gdal, osr, ogr
from netCDF4 import Dataset
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import json

def get_lls(year):
    f_name = 'DATA/RESULTS/5th_Indices_WUSA_' + str(year) + '.nc'
    ds = Dataset(f_name, 'r')
    lats = ds.variables['lat']
    lons = ds.variables['lon']
    return lats, lons

def array_to_raster(lats, lons, data_array, output_path, nodata=-9999):
    '''
    Args:
        lats: list of latitudes
        lons: list of longitues
        data_array: array of values, one per lat, lon combination
        output_path-path: path to output file
        nodata: nodata value
    '''
    ## set spatial res and size of lat, lon grids
    yres = abs(lats[1] - lats[0])
    xres = abs(lons[1] - lons[0])
    output_cols = len(lons)
    output_rows = len(lats)
    ulx = min(lons) - (xres / 2.)
    uly = max(lats) + (yres / 2.)

    ## Build the output raster file
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(
        output_path, output_cols, output_rows, 1, gdal.GDT_Float32)

    ## Convert array to float32 and set nodata value
    data_array = data_array.reshape((output_rows,output_cols)).astype(np.float32)
    data_array[data_array == nodata] = np.nan

    ## This assumes the projection is Geographic lat/lon WGS 84
    output_osr = osr.SpatialReference()
    output_osr.ImportFromEPSG(4326)
    output_ds.SetProjection(output_osr.ExportToWkt())
    output_ds.SetGeoTransform([ulx, xres, 0, uly, 0, -yres])
    output_band = output_ds.GetRasterBand(1)
    output_band.WriteArray(data_array.reshape((output_rows,output_cols)))
    output_band.SetNoDataValue(float(np.finfo(np.float32).min))
    output_ds = None

def make_tiff(in_file, out_file):
    ll_strs = []
    lons = []
    lats = []
    indices = []
    index_dict = {}
    with open(in_file, 'r') as f:
        for line in f.readlines():
            ll_vals = line.split(': ')
            ll_str = ll_vals[0]
            index = ll_vals[1].strip('\n\r').split(',')[-1]

            ll_strs.append(ll_str)
            # lon = round(float(ll_str.split(',')[0]), 5)
            # lat = round(float(ll_str.split(',')[1]), 5)
            lon = ll_str.split(',')[0]
            lat = ll_str.split(',')[1]
            lons.append(lon)
            lats.append(lat)
            #indices.append(round(float(index), 5))
            indices.append(index)
            index_dict[ll_str] =  round(float(index), 5)
    '''
    print len(lons), len(lats), len(indices)
    array_to_raster(np.array(lats), np.array(lons), np.array(indices), out_file,nodata=1e20)
    '''
    # Make data array
    data_array = []
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            ll_str = lon + ',' + lat
            if ll_str in index_dict.keys():
                data_array.append(index_dict[ll_str])
            else:
                data_array.append(-9999)
    # Convert lon,lats to floats
    lons = [round(float(l), 5) for l in lons]
    lats = [round(float(l), 5) for l in lats]
    # Rasterize
    array_to_raster(np.array(lats), np.array(lons), np.array(data_array), out_file,nodata=-9999)


def visualize_raster(tiff_file):
    tiff_file = 'DATA//RESULTS/test_indices.tiff'
    try:
        tif = gdal.Open(tiff_file)
        tifArray = tif.ReadAsArray()
    except:
        print 'Cannot open file ' + data_dir + file_name
        sys.exit(0)

    band = tif.GetRasterBand(1)
    bandArray = band.ReadAsArray()
    print type(bandArray)
    print bandArray.shape
    print bandArray.size
    imgplot = plt.imshow(bandArray)
    plt.show()

def get_plotting_data(years):
    doy_data = []
    lats, lons = get_lls(1999)
    l_data = [[] for doy_idx in range(90)]
    # for box plots at 5 locations
    ll_data = {}
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(lons):
            if lat_idx in [2,4,6,8,10] and lon_idx == lat_idx + 1:
                ll_str = str(lon) + ',' + str(lat)
                if ll_str not in ll_data.keys():
                    ll_data[ll_str] = []
                else:
                    pass
            else:
                ll_str = None
            if not ll_str: continue
            print('Processing Lon, Lat ' + ll_str)
            year_data = [[] for doy_idx in range(90)]
            for year in years:
                net_file = 'DATA/RESULTS/5th_Indices_WUSA_' + str(year) + '.nc'
                try:
                    ds = Dataset(net_file, 'r')
                except:
                    continue
                for doy_idx in range(90):
                    index = ds.variables['index'][lat_idx][lon_idx][doy_idx]
                    if index != -9999:
                        year_data[doy_idx].append(index)
            #For each doy, average over all the years
            for doy_idx in range(90):
                mn = np.mean(year_data[doy_idx])
                l_data[doy_idx].append(round(mn,4))
                ll_data[ll_str].append(round(mn,4))
    del year_data
    #Average over all lon, lats
    for doy_idx in range(90):
        doy_data.append(np.mean(l_data[doy_idx]))
    return doy_data, ll_data
########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951,2012)
    ts_data, ll_data = get_plotting_data(years)
    doys = range(1,91)
    x = np.array(doys)
    y = np.array(ts_data)
    fig = plt.figure()
    plt.plot(x,y)
    fig.savefig('DATA/RESULTS/index_over_season.png')
