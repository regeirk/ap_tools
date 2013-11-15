# -*- coding: utf-8 -*-
#
# Description: Utilities to read, manipulate and plot
#              satellite data.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com
#
# Obs: get_aviso(), sav_aviso() and plt_aviso()
#      were modified from original code written by
#      Filipe Fernandes (http://ocefpaf.tiddlyspot.com/).

__all__ = ['get_aviso',
           'sav_aviso',
           'plt_aviso']

import h5py
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date


# obs: Check with product is best for you at:
# http://opendap.aviso.oceanobs.com/thredds/catalog.html
def get_aviso(lon_min, lon_max, lat_min, lat_max, start, end, url):
    r"""Enter lon_min, lon_max, lat_min, lat_max, start, end, url."""

    nc = Dataset(url)

    lon, lat = nc.variables['NbLongitudes'][:], nc.variables['NbLatitudes'][:]

    time = nc.variables['time']
    time = num2date(time[:], time.units, calendar='standard')

    mask_lon = np.logical_and(lon >= lon_min, lon <= lon_max)
    mask_lat = np.logical_and(lat >= lat_min, lat <= lat_max)
    mask_time = np.logical_and(time >= start, time <= end)
    lon, lat, time = lon[mask_lon], lat[mask_lat], time[mask_time]

    data = nc.variables['Grid_0001'][mask_time, mask_lon, mask_lat]
    time_string = [d.strftime("%Y/%m/%d") for d in time]

    return dict(time=time_string, lon=lon, lat=lat, data=data,
                fill_value=data.fill_value)


def sav_aviso(fname, aviso):
    r"""Enter an AVISO dict with lon, lat data, and a file name."""
    with h5py.File(fname, 'w') as f:
        f['time'] = aviso['time']
        f['lon'] = aviso['lon']
        f['lat'] = aviso['lat']
        f['data'] = aviso['data']
        f['fill_value'] = aviso['fill_value']
    return None


def plt_aviso(fname, time_idx=0):
    r"""Enter am AVISO hdf5 file for simple Quick-n-Dirty plotting."""
    fig, ax = plt.subplots()
    with h5py.File(fname, 'r') as f:
        time = f['time'].value
        lon = f['lon'].value
        lat = f['lat'].value
        data = ma.masked_equal(f['data'].value, f['fill_value'].value)
    ax.pcolor(lon - 360, lat, data[time_idx, :, :].T)
    ax.set_ylabel("Latitude [degrees]")
    ax.set_xlabel("Longitude [degrees]")
    ax.set_title("%s" % time[time_idx])
    return fig, ax











if __name__ == '__main__':
    plt.close('all')
    from datetime import datetime

    lon_min, lon_max = -50., -40.
    lat_min, lat_max = -32., -22.
    # start, end = datetime(2010, 9, 1), datetime(2010, 10, 1)
    start, end = datetime(2011,1,3),datetime(2011,1,3)

    # Convert from -180--180 to 0-360.
    lon_min %= 360.
    lon_max %= 360.

    user, password = "aviso-users", "grid2010"
    uri = "opendap.aviso.oceanobs.com/thredds/dodsC/"

    # Change here for your product: "Delayed Time Data, Global, Absolute
    # Dynamic Topography, Updated Data - Daily, Absolute Dynamic Topography".
    path = "dataset-duacs-dt-upd-global-merged-msla-h-daily"

    # Final URL.
    url = "http://%s:%s@%s%s" % (user, password, uri, path)

    # Get.
    aviso = get_aviso(lon_min, lon_max, lat_min, lat_max, start, end, url)
    # Save.
    sav_aviso('test.h5', aviso)
    # Plot.
    fig, ax = plt_aviso('test.h5', time_idx=0)

    plt.show()