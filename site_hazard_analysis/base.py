import numpy as np
import pandas as pd
import geopandas as gpd
import utm
import matplotlib.pyplot as plt
import seaborn as sns
import json
import matplotlib as mpl
import os
import subprocess
import h5py

from scipy import stats
from scipy.optimize import minimize
from functools import partial
from shapely.geometry import Point, LineString, shape
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from matplotlib.colors import to_hex


# set all single line variables to be displayed, not just the last line
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

from IPython.display import display


def set_plot_formatting():
    # set up plot formatting
    SMALL_SIZE = 15
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 25

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def grayscale_version(color):
    conversion = np.array([0.299, 0.587, 0.114])

    grayscale = np.repeat(np.matmul(conversion, color), 3)

    return grayscale


def check_colors(colors, linewidth):
    fig, ax = plt.subplots(1, 1, figsize=(3, 2))
    for i in range(len(colors)):
        _ = ax.plot([0, .5], [i, i], color=colors[i], linewidth=linewidth)
        _ = ax.plot([.5, 1], [i, i], color=grayscale_version(colors[i][:-1]), linewidth=linewidth)
        _ = ax.axis('off')


def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,0,zback]])


def utm_conversion(lat, long):
    """
    Converts latitude and longitude coordinates to meters

    Notes:
           if the locations are not all in the same UTM zone,
                  this code will throw an error

    Created by Anne Hulsey, 2/20/2019

    Parameters
    ----------
    lat, long: numpy arrays
        paired arrays of latitude and longitude, in WGS84 coordinate system

    Returns
    -------
    x,y: numpy arrays
        corresponding arrays in meters, with the origin based on the UTM zone
    """

    # check whether the lat and long are a list or single site
    if len(lat.shape) == 0:
        n_sites = 1
        lat = [lat]
        long = [long]
    else:
        n_sites = len(lat)

    # convert WGS84 coordinates to meters
    x = np.zeros(n_sites)
    y = np.zeros(n_sites)
    zone = np.zeros(n_sites)
    for lat, long, i in zip(lat, long, range(len(x))):
        [x[i], y[i], zone[i], _] = utm.from_latlon(lat, long)

    # throw an error if the UTM zones are not the same
    if any(i != zone[0] for i in zone):
        raise ValueError('locations are not in the same UTM zone')

    return x, y
