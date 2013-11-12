# -*- coding: utf-8 -*-
#
# Description: Loose scripts and functions for personal use.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

__all__ = ['rot_vec',
           'lon180to360'
           'lon360to180'
           'mnear',
		   'denan',
		   'wind2stress',
		   'maxwindow',
		   'gen_dates',
		   'fmt_isobath',
		   'extract_npz',
		   'mat2npz',
		   'bb_map',
		   'wind_subset',
		   'dl_goes']

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

def rot_vec(u, v, angle=-45):
	"""
	USAGE
	-----
	u_rot,v_rot = rot_vec(u,v,angle=-45.)

	Returns the rotated vector components (`u_rot`,`v_rot`)
	from the zonal-meridional input vector components (`u`,`v`).
	The rotation is done using the angle `ang` in degrees,
	positive counterclockwise (trigonometric convention).

	Example
	-------
	>>> from matplotlib.pyplot import quiver
	>>> from ap_tools.utils import rot_vec
	>>> u = -1.
	>>> v = -1.
	>>> u2,v2 = rot_vec(u,v, angle=-30.)
	"""
	u,v = map(np.asanyarray, (u,v))
	ang = angle*np.pi/180. # Degrees to radians.

	u_rot = +u*np.cos(ang) + v*np.sin(ang) # Usually the across-shore component.
	v_rot = -u*np.sin(ang) + v*np.cos(ang) # Usually the along-shore component.
	
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
	rho_air = 1.226               # kg/m3
	mag = np.sqrt(u**2+v**2)      # m/s
	Cd = np.zeros( mag.size ) # Drag coefficient.

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

def fmt_isobath(cs):
	"""Formats the labels of isobath contours."""
	isobstrH = plt.clabel(cs, fontsize=8, fmt='%g', inline_spacing=7, manual=True)
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
	m.drawmeridians(np.arange(np.floor(lonmin), np.ceil(lonmax), 1), linewidth=0.15, labels=[1, 0, 1, 0],zorder=12)
	m.drawparallels(np.arange(np.floor(latmin), np.ceil(latmax), 1), linewidth=0.15, labels=[1, 0, 0, 0],zorder=12)
	plt.ion()
	return m

def wind_subset(times=(datetime(2009,12,28), datetime(2010,1,10)),
	dt='daily', llcrnrlon=-45, urcrnrlon=-35, llcrnrlat=-25, urcrnrlat=-18):
	"""
	USAGE
	-----
	u,v = wind_subset(...)

	Gets wind vectors from the NCDC/NOAA Blended Seawinds L4 product,
	given a given lon, lat, time bounding box.

	Input
	-----
	* times: A tuple with 2 datetime objects marking the
	  start and end of the desired subset.
	  If dt=='clim', this should be an integer indicating the
	  desired month, from 1 (January) through 12 (December).
	  This climatology is built from the 1995-2005 data.

	* dt: Time resolution wanted. Choose `six-hourly`,
	  `daily` (default) or `monthly`.

	* llcrnrlon, urcrnrlon: Longitude band wanted.
	* llcrnrlat, urcrnrlat: Latitude band wanted.

	Returns
	-------
	u,v: pandas Panels containing the data subsets.

	Example
	-------
	TODO.
	"""
	head = 'http://www.ncdc.noaa.gov/thredds/dodsC/oceanwinds'
	tail = {'six-hourly':'6hr','daily':'dly','monthly':'mon','clim':'clm-agg'}

	url = head + tail[dt]
	nc = Dataset(url)

	lon = nc.variables.pop('lon')[:]
	lon = lon360to180(lon) # Longitude is in the [0,360] range, BB is in the [-180,+180] range.
	lat = nc.variables.pop('lat')[:]
	time = nc.variables.pop('time')
	time = num2date(time[:],time.units,calendar='standard')

	# Subsetting the data.
	maskx = np.logical_and(lon >= llcrnrlon, lon <= urcrnrlon)
	masky = np.logical_and(lat >= llcrnrlat, lat <= urcrnrlat)

	# Ignoring more precise time units than
	# the resolution in the desired dataset.
	nt = time.size
	if dt=='six-hourly':
		for i in xrange(nt):
			time[i] = time[i].replace(minute=0,second=0)
	elif dt=='daily':
		for i in xrange(nt):
			time[i] = time[i].replace(hour=0,minute=0,second=0)
	elif dt=='monthly':
		for i in xrange(nt):
			time[i] = time[i].replace(day=1,hour=0,minute=0,second=0)
	elif dt=='clim':
		time[1] = time[1].replace(month=2) # Fixing February.
		aux = []; [aux.append(Time.month) for Time in time]
		time = np.copy(aux); del aux
	else:
		pass

	# If the time wanted is a single month/day/hour,
	# get the nearest time available in the dataset.
	if dt=='clim':
		maskt=time==times
	else:
		if np.size(times)<2:
			idxt = np.abs(time-times).argmin()
			maskt=time==time[idxt]
		else:
			maskt = np.logical_and(time >= times[0], time <= times[1])

	np.disp('Downloading subset...')
	lon, lat, time = lon[maskx], lat[masky], time[maskt]
	u = nc.variables['u'][maskt,0,masky,maskx]
	v = nc.variables['v'][maskt,0,masky,maskx]
	np.disp('Downloaded.')

	lat,lon,time,u,v = map(np.atleast_1d, [lat,lon,time,u,v])

	# Sets the panels.
	U = Panel(data=u, items=time, major_axis=lat, minor_axis=lon)
	V = Panel(data=v, items=time, major_axis=lat, minor_axis=lon)
	
	# Sorting major axis (longitude) to fix the crossing of the 0º meridian.
	U = U.sort_index(axis=2)
	V = V.sort_index(axis=2)

    # Retrieves the u,v arrays to fix masked values.
	u,v,lon,lat = map(np.atleast_1d, (U.values,V.values,U.minor_axis.values,U.major_axis.values))
	fill_value = -9999.
	u = np.ma.masked_equal(u, fill_value)
	v = np.ma.masked_equal(v, fill_value)

	# Resets the Panels with the sorted and masked u,v arrays.
	U = Panel(data=u, items=time, major_axis=lat, minor_axis=lon)
	V = Panel(data=v, items=time, major_axis=lat, minor_axis=lon)

	return U,V

def dl_goes(time=datetime(2013,9,13), dt=24, dest_dir='./'):
	"""
	USAGE
	-----
	dl_goes(time=datetime(2013,9,13))

	Downloads full GOES SST images from the PODAAC FTP (12 MB each). Uses wget.

	ftp://podaac-ftp.jpl.nasa.gov/allData/goes/L3/goes_6km_nrt/americas/.

	* `time` is a datetime object or a list of datetime objects containing the desired times.

	* `dt` is an integer representing the time resolution desired. Choose from `24` (default),
	`3` or `1`. If `dt` is either 1 or 3, all images for each day in `time` will be downloaded.

	* `dest_dir` is the directory in which the downloaded data will be saved.

	TODO
	----
	Find an openDAP link for this dataset (WITH the bayesian cloud filter).
	"""

	if type(time)!=list:
		time = [time]

	original_dir = os.getcwd()  # Store original working directory.
	if os.path.isdir(dest_dir): # If dest_dir already exists.
		os.chdir(dest_dir)
	else:                       # Create it if it does not exist.
		os.makedirs(dest_dir)
		os.chdir(dest_dir)

	for date in time: # Getting files for each day in the list.
		yyyy = str(date.year)
		dd = int(datetime2doy(date)) # Get the julian day.
		dd = str(dd).zfill(3)
		head = 'ftp://podaac-ftp.jpl.nasa.gov/OceanTemperature/goes/L3/goes_6km_nrt/americas/%s/%s/' %(yyyy,dd)
		filename = 'sst%s?_%s_%s' % (str(dt),yyyy,dd) # dt can be 1, 3 or 24 (hourly, 3-hourly or daily).
		url = head + filename                         # The 'b' character is only for 2008-present data.
		cmd = "wget --tries=inf %s" %url
		os.system(cmd) # Download file.

	os.chdir(original_dir) # Return to the original working directory.
	np.disp("Done downloading all files.")

	return None

# TODO.
# def mur_subset(times=(datetime(2009,12,28), datetime(2010,1,10)),
# 	llcrnrlon=-45, urcrnrlon=-35, llcrnrlat=-25, urcrnrlat=-18, return_panels=True):
# 	"""
# 	USAGE
# 	-----
# 	sst = mur_subset(..., return_panels=True)

# 	*OR*

# 	[lon,lat,time,sst] = mur_subset(..., return_panels=False)

# 	Gets sst fields from the Multi-scale Ultra-high Resolution (MUR)
# 	L4 product, given a given lon, lat, time bounding box.

# 	http://mur.jpl.nasa.gov/

# 	Input
# 	-----
# 	* times: A tuple with 2 datetime objects marking the
# 	  start and end of the desired subset.

# 	* dt: Time resolution wanted. Choose `six-hourly`,
# 	  `daily` (default) or `monthly`.

# 	* llcrnrlon, urcrnrlon: Longitude band wanted.
# 	* llcrnrlat, urcrnrlat: Latitude band wanted.

# 	* return_panels: Whether or not to return data subsets
# 	  as a pandas Panels. Returns lon,lat,time,u,v numpy
# 	  arrays if `False`. Defalts to `True`.

# 	Returns
# 	-------
# 	u,v: pandas Panels containing the data subsets.

# 	*OR*

# 	lon,lat,time,u,v: 1D numpy arrays containing the data subset.

# 	Example
# 	-------
# 	TODO.
# 	"""
# 	head = 'http://www.ncdc.noaa.gov/thredds/dodsC/oceanwinds'

# 	url = head + tail[dt]
# 	nc = Dataset(url)

# 	llcrnrlon, urcrnrlon = map(lon180to360, (llcrnrlon, urcrnrlon))
# 	lon = lon180to360(nc.variables.pop('lon')[:])
# 	lat = nc.variables.pop('lat')[:]
# 	time = nc.variables.pop('time')
# 	time = num2date(time[:], time.units, calendar='standard')

# 	# Subsetting the data.
# 	maskx = np.logical_and(lon >= llcrnrlon, lon <= urcrnrlon)
# 	masky = np.logical_and(lat >= llcrnrlat, lat <= urcrnrlat)

# 	# If the time wanted is a single month/day/hour,
# 	# get the nearest time available in the dataset.
# 	if np.size(times)<2:
# 		idxt = np.abs(time-times[0]).argmin()
# 		maskt = time==time[idxt]
# 	else:
# 		maskt = np.logical_and(time >= times[0], time <= times[1])

# 	np.disp('Downloading wind subset...')
# 	lon, lat, time = lon[maskx], lat[masky], time[maskt]
# 	sst = nc.variables['analyzed_sst'][maskt,0,masky,maskx]
# 	np.disp('Done.')

# 	lon = lon360to180(lon)
# 	lat,lon,time,u,v = map(np.atleast_1d, [lat,lon,time,sst])

# 	# Returns data as pandas Panels or as numpy arrays.
# 	if return_panels:
# 		sst = Panel(data=sst, items=time, major_axis=lat, minor_axis=lon)
# 		return sst
# 	else:
# 		return lon,lat,time,sst

def helloworld():
  np.disp('Hello there, Mr. World.')

if __name__=='__main__':
  helloworld()
  import doctest
  doctest.testmod()