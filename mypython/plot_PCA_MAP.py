import sys
import numpy as np
#GEOSPATIAL STUFF
from osgeo import gdal, osr, ogr
from netCDF4 import Dataset
import matplotlib
import json

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

########
#M A I N
########
if __name__ == '__main__' :
    years = range(1951,2012)
    years = range(1951, 2012)
    data_dir = 'RESULTS/'
    in_file =  data_dir + 'PCA_raw_data.json'
    PCA_data = json.loads(in_file)
    lons = PCA_data['lons']
    lats = PCA_data['lats']
    for c_idx in range(PCA_data['comp_matrix_r'].shape[0]):
        out_file = data_dir + 'PCA_comp_' +  str(c_idx) + '.tiff'
        comp = PCA_data['comp_matrix_r'][c_idx]
        array_to_raster(lats, lons, comp, out_file,nodata=-9999)
