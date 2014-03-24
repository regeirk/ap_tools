# -*- coding: utf-8 -*-
#
# Description: Utilities to work with datasets, mostly
#              for local, personal use.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com
#
# Obs:

__all__ = ['topo_subset']

import numpy as np
from netCDF4 import Dataset
from oceans.datasets import get_indices

def topo_subset(llcrnrlon=-42, urcrnrlon=-35, llcrnrlat=-23,
                urcrnrlat=-14, tfile='/home/andre/lado/general_data/smith_sandwell/topo30sec.grd'):
    """
    Get a subset from an etopo1, etopo2 or Smith and Sandwell topography file.

    OBS: Modified from oceans.datasets.etopo_subset() FUNCTION BY Filipe Fernandes.
    """
    topo = Dataset(tfile, 'r')

    if 'smith_sandwell' in tfile.lower() or 'etopo1' in tfile.lower():
        lons = topo.variables["lon"][:]
        lats = topo.variables["lat"][:]
    elif 'etopo2' in tfile.lower():
        lons = topo.variables["x"][:]
        lats = topo.variables["y"][:]
    else:
        np.disp('Unknown topography file.')
        return

    res = get_indices(llcrnrlat, urcrnrlat, llcrnrlon, urcrnrlon, lons, lats)

    lon, lat = np.meshgrid(lons[res[0]:res[1]], lats[res[2]:res[3]])

    bathy = topo.variables["z"][int(res[2]):int(res[3]),
                                int(res[0]):int(res[1])]

    return lon, lat, bathy