#!/usr/bin/python
import matplotlib.dates as mdates
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import json
import datetime as dt
# import pandas as pd


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


def advance_date(date_dt, days, back_or_forward):
    if back_or_forward == 'forward':
        d_dt_new = date_dt + dt.timedelta(days=int(days))
    if back_or_forward == 'back':
        d_dt_new = date_dt - dt.timedelta(days=int(days))
    return d_dt_new
########
# M A I N
########


if __name__ == '__main__':
    in_file = '/Users/bdaudert/WEBPAGES/COLD/media/json/TEST/PCA_index_data.json'
    with open(in_file, 'r') as f:
        Phi = np.array(json.load(f))
    fracs, comp_matrix = compute_PCA(Phi)
    '''
    Phi_r = varimax(Phi)
    fracs_r, comp_matrix_r = compute_PCA(Phi_r)
    '''
    comp_matrix_r = varimax(comp_matrix)

    '''
    comp_1 = np.array([item[0] for item in comp_matrix_r])
    print comp_1.shape
    ts = comp_1.dot(Phi)
    print len(ts)
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
    x = []
    y = []
    z = []
    for item in comp_matrix_r:
        x.append(item[0])
        y.append(item[1])
        z.append(item[2])

    plt.close('all') # close all latent plotting windows
    fig1 = plt.figure() # Make a plotting figure
    ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
    pltData = [x,y,z]
    # make a scatter plot of blue dots from the data
    ax.scatter(pltData[0], pltData[1], pltData[2], 'bo')

    # make simple, bare axis lines through space:
    # 2 points make the x-axis line at the data extrema along x-axis
    xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0))
    # make a red line for the x-axis.
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
    # 2 points make the y-axis line at the data extrema along y-axis
    yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0))
    # make a red line for the y-axis.
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
    # 2 points make the z-axis line at the data extrema along z-axis
    zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2])))
    # make a red line for the z-axis.
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')

    # label the axes
    ax.set_xlabel("pc 1")
    ax.set_ylabel("pc 2")
    ax.set_zlabel("pc 3")
    ax.set_title("Rotated PCA first three components")
    plt.show() # show the plot
    '''
