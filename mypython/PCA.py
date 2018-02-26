#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib.dates as mdates
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.decomposition import PCA
import json
import datetime as dt
from osgeo import gdal, osr, ogr
# import pandas as pd

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

def ortho_rotation(lam, method='varimax',gamma=None,
                   eps=1e-6, itermax=100):
    """
    Return orthogal rotation matrix
    TODO: - other types beyond
    """
    if gamma == None:
        if method == 'varimax':
            gamma = 1.0
        if method == 'quartimax':
            gamma = 0.0

    nrow, ncol = lam.shape
    R = np.eye(ncol)
    var = 0

    for i in range(itermax):
        lam_rot = np.dot(lam, R)
        tmp = np.diag(np.sum(lam_rot ** 2, axis=0)) / nrow * gamma
        u, s, v = np.linalg.svd(np.dot(lam.T, lam_rot ** 3 - np.dot(lam_rot, tmp)))
        R = np.dot(u, v)
        var_new = np.sum(s)
        if var_new < var * (1 + eps):
            break
        var = var_new

    return R

def varimax(Phi, gamma=1.0, q=20, tol=1e-6):
    from numpy import eye, asarray, dot, sum, diag
    from numpy.linalg import svd
    p, k = Phi.shape
    R = eye(k)
    d = 0
    for i in xrange(q):
        d_old = d
        Lambda = dot(Phi, R)
        u, s, vh = svd(dot(Phi.T, asarray(Lambda) ** 3 - (gamma / p) *
                       dot(Lambda, diag(diag(dot(Lambda.T, Lambda))))))
        R = dot(u, vh)
        d = sum(s)
        if d / d_old < tol:
            break
    return dot(Phi, R)

def datetime_to_date(dtime, seperator):
    '''
    yyyy-mm-dd
    yyyy/mm/dd
    yyyy:mm:dd
    yyyymmdd
    '''
    if type(dtime) != dt.datetime:
        return '0000' + str(seperator) + '00' + str(seperator) + '00'
    try:y = str(dtime.year)
    except:y = '0000'

    try:m =str(dtime.month)
    except:m = '00'
    if len(m) == 1:m = '0' + m

    try:d =str(dtime.day)
    except:d = '00'
    if len(d) == 1:d = '0' + d
    return y + str(seperator) + m + str(seperator) + d

def advance_date(date_dt, days, back_or_forward):
    if back_or_forward == 'forward':
        d_dt_new = date_dt + dt.timedelta(days=int(days))
    if back_or_forward == 'back':
        d_dt_new = date_dt - dt.timedelta(days=int(days))
    return d_dt_new

def get_PCA_data(years, data_dir):
    '''
    Row are time steps, columns are locations.
    For our domain
    PCA_data.shape = (5490, 103968)
    where 5490 = 61 years * 90 days
    and 103968 is number of locations in domain
    '''
    lons = np.array([])
    lats = np.array([])
    PCA_data = []
    for year in years:
        DS = Dataset(data_dir + '5th_Indices_WUSA_' + str(year) + '.nc', 'r')
        if lons.shape[0] == 0:
            lons = np.array([l -360 for l in DS.variables['lon'][:]])
            lats = DS.variables['lat'][:]
            num_lons = lons.shape[0]
            num_lats = lats.shape[0]
        indices = DS.variables['index'][:,:,:]
        doys = DS.variables['doy'][:]
        for doy_idx in range(doys.shape[0]):
            d_list = np.reshape(indices[doy_idx,:,:], num_lons * num_lats)
            d_list[d_list == -9999] = 0
            PCA_data.append(d_list)
    return np.array(PCA_data), lons, lats

def scale_linear_bycolumn(rawpoints, high=1.0, low=0.0):
    mins = np.min(rawpoints, axis=0)
    maxs = np.max(rawpoints, axis=0)
    rng = maxs - mins
    print('RNG ' + str(rng))
    print(maxs - rawpoints)
    return high - (((high - low) * (maxs - rawpoints)) / rng)

########
# M A I N
########


if __name__ == '__main__':
    start_time = dt.datetime.now()
    print('Start Time: ' + str(start_time))

    years = range(1951, 2012)
    data_dir = 'RESULTS/LIVNEH/'

    print('Extracting PCA data from index files')

    #PCA_data, lons, lats = get_PCA_data(years, data_dir)
    PCA_data, lons, lats = get_PCA_data(years, data_dir)
    print('DATA MATRIX ' + str(PCA_data.shape))

    # Scale PCA data
    # Gives zero/NaN/infinity error
    #PCA_data = scale_linear_bycolumn(PCA_data, high=1.0, low=0.0)

    #print('Minutes elapsed: ' + str((dt.datetime.now() - start_time).total_seconds() / 60.0))
    print('Computing component matrix')

    pca = PCA(n_components=3)
    X_pca = pca.fit_transform(PCA_data)
    #X_pca = pca.fit(PCA_data)
    print('X_pca ' + str(X_pca.shape))
    components = pca.components_.transpose()

    print('COMPONENTS ' + str(components.shape))
    #print('Minutes elapsed: ' + str((dt.datetime.now() - start_time).total_seconds() / 60.0))
    print('Computing rotated components')

    rotated_components = varimax(components).transpose()
    #print('Minutes elapsed: ' + str((dt.datetime.now() - start_time).total_seconds() / 60.0))
    print('ROTATED COMPONENTS ' + str(rotated_components.shape))

    dates_dt = []
    dates_ts = []
    for year in years:
        dates_dt.append(dt.datetime(year,12,1))
        dates_ts.append(datetime_to_date(dates_dt[-1], '-'))
        for doy_idx in range(1,90):
            dates_dt.append(advance_date(dates_dt[-1],1, 'forward'))
            dates_ts.append(datetime_to_date(dates_dt[-1], '-'))

    rotated_components_list = []
    for c_idx in range(3):
        print('Working component ' + str(c_idx + 1))
        #comp = [round(c, 4) for c in rotated_components[c_idx]]
        #rotated_components_list.append(comp)
        #Only get time series for first threee components
        print('Finding time series data')
        ts = rotated_components[c_idx].dot(PCA_data.T)
        hc_ts_data = []
        for i in range(ts.shape[0]):
            hc_ts_data.append([dates_ts[i], round(ts[i], 4)])
        #print('Minutes elapsed: ' + str((dt.datetime.now() - start_time).total_seconds() / 60.0))
        print('Saving time series data')
        # Save time series data
        with open(data_dir + 'PCA_COMPONENT_' + str(c_idx+1) + '_TS.json', 'w') as outfile:
            json.dump(hc_ts_data, outfile)
        print('Saving tiff')
        #Generate tiff for component
        out_file = data_dir + 'PCA_COMPONENT_' + str(c_idx + 1) + '_MAP.tiff'
        array_to_raster(lats, lons, rotated_components[c_idx], out_file, nodata=-9999)
    '''
    # Save PCA data for plotting OR write tiff here
    pca_json = {
        'fracs':[round(f, 4) for f in fracs],
        'rotated_components':rotated_components_list,
        'lons': [round(l, 4) for l in list(lons)],
        'lats': [round(l, 4) for l in list(lats)]
    }
    with open(data_dir + 'PCA_1_3_raw_data.json', 'w') as outfile:
        json.dump(pca_json, outfile)
    print('Minutes elapsed: ' + str((dt.datetime.now() - start_time).total_seconds() / 60.0))
    '''

    '''
    dates = mdates.date2num(dates_dt)
    plt.plot_date(dates, ts)
    plt.show()
    plt.savefig('test_PCA.png')
    '''