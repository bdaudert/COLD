#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import json
from osgeo import gdal, osr, ogr

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
    output_ds.SetGeoTransform([ulx, xres, 0, uly, 0, yres])
    output_band = output_ds.GetRasterBand(1)
    output_band.WriteArray(data_array)
    output_band.SetNoDataValue(float(np.finfo(np.float32).min))
    output_ds = None

def get_array_data(var_name, years, data_dir):
    all_percs = np.array([])
    DS = Dataset(data_dir + var_name +  '_percentiles.nc', 'r')
    lons = np.array([l -360 for l in DS.variables['lon'][:]])
    lats = DS.variables['lat'][:]
    percs = DS.variables['percentile'][:,:]
    return percs, lons, lats

if __name__ == '__main__':
    years = range(1951, 2012)
    data_dir = 'RESULTS/livneh/'
    for var_name in ['tmin', 'tmax']:
        percs, lons, lats = get_array_data(var_name, years, data_dir)
        # out_file = data_dir + var_name + '_percentiles.tif'
        out_file = data_dir + var_name + '_percentiles.tif'
        array_to_raster(lats, lons, percs, out_file, nodata=-9999)
