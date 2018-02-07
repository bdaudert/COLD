#!/usr/bin/python
#from matplotlib.mlab import PCA
# import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import json


def compute_PCA(Phi):
    # Phi =  2D np.array, rows doy index vals, columns locs
    # numpy array: [[row1], [row2]]
    # construct your numpy array of data
    results = PCA(Phi)

    # this will return an array of variance percentages for each component
    fracs = results.fracs

    # this will return a 2d array of the data projected into PCA space
    component_matrix = results.Y
    return fracs, component_matrix

# Rotated PCA??
# https://stackoverflow.com/questions/17628589/perform-varimax-rotation-in-python-using-numpy


def varimax(Phi, gamma=1, q=20, tol=1e-6):
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

def get_PCA_data(years, data_dir):
    PCA_data = np.array([])
    for year in years:
        DS = Dataset(data_dir + '5th_Indices_WUSA_' + str(year) + '.nc', 'r')
        num_lons = DS.variables['lon'][:].shape[0]
        num_lats = DS.variables['lat'][:].shape[0]
        indices = DS.variables['index'][:,:,:]
        doys = DS.variables['doy'][:]
        for doy_idx in range(doys.shape[0]):
            np.append(PCA_data, np.reshape(indices[doy_idx,:,:], num_lons * num_lats))
    return PCA_data
########
# M A I N
########


if __name__ == '__main__':
    '''
    in_file = 'DATA/TEST_DATA/PCA_index_data.json'
    with open(in_file, 'r') as f:
        Phi = np.array(json.load(f))
    '''
    years = range(1951, 2012)
    data_dir = 'RESULTS/'
    Phi = get_PCA_data(years, data_dir)

    fracs, comp_matrix = compute_PCA(Phi)
    '''
    Phi_r = varimax(Phi)
    fracs_r, comp_matrix_r = compute_PCA(Phi_r)
    '''
    comp_matrix_r = varimax(comp_matrix)

    print comp_matrix_r
    print fracs

    '''
    comp_1 = comp_matrix_r[0]
    ts = Phi.dot(comp_1)
    print ts
    print ts.shape
    dates_dt = []
    years = range(1950,2011)
    for year in years:
        dates_dt.append(dt.datetime(year,12,1))
        for doy_idx in range(1,90):
            dates_dt.append(advance_date(dates_dt[-1],1, 'forward'))
    dates = mdates.date2num(dates_dt)
    plt.plot_date(dates, ts)
    plt.show()
    '''
