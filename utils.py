# -*- coding: utf-8 -*-
#
# Description: Loose scripts and functions for personal use.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

__all__ = ['rot_vec',
           'lon180to360',
           'lon360to180',
           'mnear',
		   'denan',
		   'wind2stress',
		   'maxwindow',
		   'gen_dates',
		   'fmt_isobath',
		   'extract_npz',
		   'mat2npz',
		   'bb_map',
		   'dots_dualcolor']

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
from times import doy2datetime,datetime2doy
from dateutil import rrule,parser
from scipy.io import loadmat,savemat
from glob import glob
from mlabwrap import MatlabPipe
from netCDF4 import Dataset, num2date
from pandas import Panel
from gsw import distance

def rot_vec(u, v, angle=-45, degrees=True):
	"""
	USAGE
	-----
	u_rot,v_rot = rot_vec(u,v,angle=-45.,degrees=True)

	Returns the rotated vector components (`u_rot`,`v_rot`)
	from the zonal-meridional input vector components (`u`,`v`).
	The rotation is done using the angle `angle` positive counterclockwise
	(trigonometric convention). If `degrees` is set to `True``(default),
	then `angle` is converted to radians.
	is

	Example
	-------
	>>> from matplotlib.pyplot import quiver
	>>> from ap_tools.utils import rot_vec
	>>> u = -1.
	>>> v = -1.
	>>> u2,v2 = rot_vec(u,v, angle=-30.)
	"""
	u,v = map(np.asanyarray, (u,v))
	if degrees:
		angle = angle*np.pi/180. # Degrees to radians.

	u_rot = +u*np.cos(angle) + v*np.sin(angle) # Usually the across-shore component.
	v_rot = -u*np.sin(angle) + v*np.cos(angle) # Usually the along-shore component.
	
	return u_rot,v_rot

def lon180to360(lon):
	"""
	Converts longitude values in the range [-180,+180]
	to longitude values in the range [0,360].
	"""
	lon = np.asanyarray(lon)
	return (lon + 360.0) % 360.0

def lon360to180(lon):
	"""
	Converts longitude values in the range [0,360]
	to longitude values in the range [-180,+180].
	"""
	lon = np.asanyarray(lon)
	return ((lon + 180.) % 360.) - 180.

def mnear(x, y, x0, y0):
	"""
	USAGE
	-----
	xmin,ymin = mnear(x, y, x0, y0)

	Finds the the point in a (lons,lats) line
	that is closest to a given (lon0,lat0) point.
	"""
	x,y,x0,y0 = map(np.asanyarray, (x,y,x0,y0))
	point = (x0,y0)

	d = np.array([])
	for n in xrange(x.size):
		xn,yn = x[n],y[n]
		dn = distance((xn,x0),(yn,y0)) # Calculate distance point-wise.
		d = np.append(d,dn)

	idx = d.argmin()

	return x[idx],y[idx]

# def subset(x, y, arr, xmin, xmax, ymin, ymax):
# 	"""
# 	USAGE
# 	-----
# 	x,y,z = subset(X,Y,Z,xmin,xmax,ymin,ymax)

# 	Returns a subset of an array `Z(X,Y)`
# 	and its subsetted axes `x` and `y`.
# 	"""
# 	if np.logical_and(x.ndim==1,y.ndim==1):
# 		x,y = np.meshgrid(x,y)

# 	fx = np.logical_and(x >= xmin, x <= xmax)
# 	fy = np.logical_and(y >= ymin, y <= ymax)
# 	arr = arr[fx,fy]
# 	x = x[fx,fy]
# 	y = y[fx,fy]
# 	return x,y,arr

def denan(arr):
	"""
	USAGE
	-----
	denaned_arr = denan(arr)

	Remove the NaNs from an array.
	"""
	f = np.isnan(arr)
	return arr[~f]

def maxwindow():
	"""
	Maximizes current window.
	"""
	backend = matplotlib.backends.backend
	mng = plt.get_current_fig_manager()

	if backend=='TkAgg':                  # Tk window manager.
		mng.resize(*mng.window.maxsize())
	elif backend=='WXAgg':                # Wx window manager.
		mng.frame.Maximize(True)
	elif backend=='QT4Agg':               # Qt window manager.
		mng.window.showMaximized()

def wind2stress(u, v, formula='large_pond1981-modified'):
	"""
	USAGE
	-----
	taux,tauy = wind2stress(u, v, formula='mellor2004')

	Converts u,v wind vector components to taux,tauy
	wind stress vector components.
	"""
	rho_air = 1.226            # kg/m3
	mag = np.sqrt(u**2+v**2)   # m/s
	Cd = np.zeros( mag.shape ) # Drag coefficient.

	if formula=='large_pond1981-modified':
		# Large and Pond (1981) formula
		# modified for light winds, as
		# in Trenberth et al. (1990).
		f=mag<=1.
		Cd[f] = 2.18e-3
		f=np.logical_and(mag>1.,mag<3.)
		Cd[f] = (0.62+1.56/mag[f])*1e-3
		f=np.logical_and(mag>=3.,mag<10.)
		Cd[f] = 1.14e-3
		f=mag>=10.
		Cd[f] = (0.49 + 0.065*mag[f])*1e-3
	elif formula=='mellor2004':
		Cd = 7.5e-4 + 6.7e-5*mag
	else:
		np.disp('Unknown formula for Cd.')
		pass

	# Computing wind stress [N/m2]
	taux = rho_air*Cd*mag*u
	tauy = rho_air*Cd*mag*v

	return taux,tauy

def gen_dates(start, end, dt='day', input_datetime=False):
	"""
	Returns a list of datetimes within the date range
	from `start` to `end`, at a `dt` time interval.

	`dt` can be 'second', 'minute', 'hour', 'day', 'week',
	'month' or 'year'.

	If `input_datetime` is False (default), `start` and `end`
	must be a date in string form. If `input_datetime` is True,
	`start` and `end` must be datetime objects.

	Note
	----
	Modified from original function
	by Filipe Fernandes (ocefpaf@gmail.com).

	Example
	-------
	>>> from ap_tools.utils import gen_dates
	>>> from datetime import datetime
	>>> start = '1989-08-19'
	>>> end = datetime.utcnow().strftime("%Y-%m-%d")
	>>> gen_dates(start, end, dt='day')
	"""
	DT = dict(second=rrule.SECONDLY,
		      minute=rrule.MINUTELY,
		      hour=rrule.HOURLY,
		      day=rrule.DAILY,
		      week=rrule.WEEKLY,
		      month=rrule.MONTHLY,
		      year=rrule.YEARLY)

	dt = DT[dt]

	if input_datetime: # Input are datetime objects. No parsing needed.
		dates = rrule.rrule(dt, dtstart=start, until=end)
	else:              # Input in string form, parse into datetime objects.
		dates = rrule.rrule(dt, dtstart=parser.parse(start), until=parser.parse(end))
	return list(dates)

def fmt_isobath(cs, fontsize=8, fmt='%g', inline=True, inline_spacing=7, manual=True):
	"""
	Formats the labels of isobath contours. `manual` is set to `True` by default,
	but can be `False`, or a tuple/list of tuples with the coordinates of the labels.
	All options are passed to plt.clabel().
	"""
	isobstrH = plt.clabel(cs, fontsize=fontsize, fmt=fmt, inline=inline, inline_spacing=inline_spacing, manual=manual)
	for ih in range(0, len(isobstrH)): # Appends 'm' for meters at the end of the label.
		isobstrh = isobstrH[ih]
		isobstr = isobstrh.get_text()
		isobstr = isobstr.replace('-','') + ' m'
		isobstrh.set_text(isobstr)

def extract_npz(fname):
	"""
	USAGE
	-----
	varlist = extract_npz(fname)

	Extract variables stored in a .npz file,
	and returns them in a list.
	"""
	d = np.load(fname)
	arrs = []
	for arr in d.iterkeys():
		exec("arrs.append(d['%s'])"%arr)
	return arrs

def mat2npz(matname):
	"""
	USAGE
	-----
	mat2npz(matname)

	Extract variables stored in a .mat file,
	and saves them in a .npz file.
	"""
	d = loadmat(matname)
	_ = d.pop('__header__')
	_ = d.pop('__globals__')
	_ = d.pop('__version__')
	npzname = matname[:-4] + '.npz'
	np.savez(npzname,**d)
	return None

def bb_map(lons, lats, projection='merc', resolution='i'):
	"""
	USAGE
	-----
	m = bb_map(lons, lats, **kwargs)

	Returns a Basemap instance with lon,lat bounding limits
	inferred from the input arrays `lons`,`lats`.
	Coastlines, countries, states, parallels and meridians
	are drawn, and continents are filled.
	"""
	lons,lats = map(np.asanyarray, (lons,lats))
	lonmin,lonmax = lons.min(),lons.max()
	latmin,latmax = lats.min(),lats.max()

	m = Basemap(llcrnrlon=lonmin,
				urcrnrlon=lonmax,
				llcrnrlat=latmin,
				urcrnrlat=latmax,
				projection=projection,
				resolution=resolution)

	plt.ioff() # Avoid showing the figure.
	m.fillcontinents(color='0.9', zorder=9)
	m.drawcoastlines(zorder=10)
	m.drawstates(zorder=10)
	m.drawcountries(linewidth=2.0, zorder=10)
	m.drawmapboundary(zorder=9999)
	m.drawmeridians(np.arange(np.floor(lonmin), np.ceil(lonmax), 1), linewidth=0.15, labels=[1, 0, 1, 0], zorder=12)
	m.drawparallels(np.arange(np.floor(latmin), np.ceil(latmax), 1), linewidth=0.15, labels=[1, 0, 0, 0], zorder=12)
	plt.ion()
	return m

def dots_dualcolor(x, y, z, thresh=20., color_low='b', color_high='r', marker='o', markersize=5):
	"""
	USAGE
	-----
    dots_dualcolor(x, y, z, thresh=20., color_low='b', color_high='r')

	Plots dots colored with a dual-color criterion,
	separated by a threshold value.
	"""
	ax = plt.gca()
	# Below-threshold dots.
	f=z<=thresh
	ax.plot(x[f], y[f], lw=0, marker=marker, ms=markersize, mfc=color_low, mec=color_low)
	# Above-threshold dots.
	f=z>thresh
	ax.plot(x[f], y[f], lw=0, marker=marker, ms=markersize, mfc=color_high, mec=color_high)






def helloworld():
  np.disp('Hello there, Mr. World.')

if __name__=='__main__':
  helloworld()
  import doctest
  doctest.testmod()