# -*- coding: utf-8 -*-
#
# Description: Loose scripts and functions for personal use.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
from glob import glob
from mlabwrap import MatlabPipe
from oceans.ff_tools import wrap_lon180, wrap_lon360
from netCDF4 import Dataset, num2date

def blendedseawinds_subset(times=[datetime(2010,1,1),datetime(2010,1,2),datetime(2010,1,3)],
	dt='daily', llcrnrlon=-45, urcrnrlon=-38, llcrnrlat=-25, urcrnrlat=-18):
	"""
	Get wind vectors from the NCDC-NOAA Blended Seawinds L4 product,
	at a given lon, lat, time bounding box.

	USAGE
	-----
	[lon,lat,time,u,v] = blendedseawinds_subset(...)

	Input
	-----
	times: Datetime object or list of datetime objects.

	dt: Time resolution wanted. Choose `six-hourly`,
	`daily` (default) or `monthly`.

	llcrnrlon, urcrnrlon: Longitude band wanted.
	llcrnrlat, urcrnrlat: Latitude band wanted.

	Returns
	-------
	lon,
	lat,
	time,
	u,
	v

	Example
	-------
	TODO.
	"""

	head = 'http://www.ncdc.noaa.gov/thredds/dodsC/oceanwinds'
	tails = {'six-hourly':'6hr','daily':'dly','monthly':'mon'}

	url = head + tails[dt]
	nc = Dataset(url)

	llcrnrlon, urcrnrlon = map(wrap_lon360, (llcrnrlon, urcrnrlon))
	lon = wrap_lon360(nc.variables.pop('lon')[:])
	lat = nc.variables.pop('lat')[:]
	time = nc.variables.pop('time')
	time = num2date(time[:], time.units, calendar='standard')

	# if dt=='daily':
	# 	tcorr = times
	# 	times += tcorr
	# elif dt=='monthly':


	# Subsettind the data.
	maskx = np.logical_and(lon >= llcrnrlon, lon <= urcrnrlon)
	masky = np.logical_and(lat >= llcrnrlat, lat <= urcrnrlat)

	times = np.asanyarray(times) # Wanted times.

	# If the time wanted is an interval of months/days/hours.
	if times.ndim>0:
		maskt = np.logical_and(time >= times[0], time <= times[-1])
	# If the time wanted is a single month/day/hour,
	# get the nearest time available.
	else:
		idxt = np.abs(time-times).argmin()
		maskt = time[idxt]

	lon, lat, time = lon[maskx], lat[masky], time[maskt]
	u = nc.variables['u'][maskt,0,masky,maskx]
	v = nc.variables['v'][maskt,0,masky,maskx]

	# Returns pandas Panel .
	if return_panel:
		data = Panel4D(subset, major_axis=lat, minor_axis=lon,
		               labels=np.atleast_1d(times),
		               items=np.atleast_1d(depth))
	# Returns numpy arrays.
	else:
		return lon,lat,time,u,v

# ABORTED.
# def goes_subset(time=datetime(2010,1,2), data_dir='/home/andre/.goes_data/', dt=24, cloud_thresh=2.0, matlab_process_path='/usr/local/matlabr2008b/bin/matlab', matlab_version='2008b'):
# 	"""Downloads a subset of GOES SST data from the PODAAC FTP.
# 	ftp://podaac-ftp.jpl.nasa.gov/allData/goes/L3/goes_6km_nrt/americas/."""

# 	if type(time)!=list:
# 		time = [time]

#     # Getting files for each day.
# 	for date in time:
# 		yyyy = str(date.year)
# 		dd = str(date.day).zfill(3)
# 		head = 'ftp://podaac-ftp.jpl.nasa.gov/OceanTemperature/goes/L3/goes_6km_nrt/americas/%s/%s/' %(yyyy,dd)
# 		filename = 'sst%sb_%s_%s' % (str(dt),yyyy,dd) # Interval can be 1, 3 or 24 (hourly, 3-hourly or daily).
# 		url = head + filename                         # The 'b' character is only for 2008-present data.
# 		cmd = "wget -r --tries=inf %s" %url
# 		original_dir = os.getcwd()
# 		if os.path.isdir(data_dir):
# 			os.chdir(data_dir)
# 		else:
# 			os.makedirs(data_dir)
# 			os.chdir(data_dir)

# 		os.system(cmd) # Download file.
# 		filelist = glob(filename); filelist.sort()
# 		lz = len(filelist)*len(time)
# 		# Run read_goes.m and subset SST array.
# 		for fname in filelist:
# 			matlab = MatlabPipe(matlab_process_path=matlab_process_path, matlab_version=matlab_version)
# 			matlab.open(print_matlab_welcome=False)
# 			o=matlab.eval("pwd;ls")
# 			cmd = "[sst,x,y] = read_goes('%s', %s);" %(fname,str(cloud_thresh))
# 			print cmd
# 			o=matlab.eval(cmd)
# 			o=matlab.eval("sst = sst';") # !
# 			o=matlab.eval("[xg,yg] = meshgrid(x,y);")
# 			o=matlab.eval("x_aux = (x>=-44.8) & (x<=-37.8); y_aux = (y>=-25.2) & (y<=-19.8);")
# 			o=matlab.eval("x_aux2 = (x>=-42) & (x<=-41); y_aux2 = (y>=-22.8) & (y<=-22);")
# 			o=matlab.eval("shp = [sum(y_aux) sum(x_aux)];")
# 			o=matlab.eval("shp2 = [sum(y_aux2) sum(x_aux2)];")
# 			o=matlab.eval("fcst = (xg>=-44.8) & (xg<=-37.8) & (yg>=-25.2) & (yg<=-19.8);")
# 			o=matlab.eval("fcst2 = (xg>=-42) & (xg<=-41) & (yg>=-22.8) & (yg<=-22);")
# 			o=matlab.eval("sst1 = reshape(sst(fcst),shp);")
# 			o=matlab.eval("sst2 = reshape(sst(fcst2),shp2);")
# 			o=matlab.eval("sst = sst1; clear sst1")
# 			o=matlab.eval("x = reshape(xg(fcst),shp);")
# 			o=matlab.eval("y = reshape(yg(fcst),shp);")
# 			o=matlab.eval("x2 = reshape(xg(fcst2),shp2);")
# 			o=matlab.eval("y2 = reshape(yg(fcst2),shp2);")
# 			# o=matlab.eval("save temp.mat x y x2 y2 sst sst2")
# 			# o=matlab.eval("save temp.mat x y sst")
# 			matlab.close()
# 	os.chdir(original_dir)

# 	return None

def helloworld():
  np.disp('Hello there, Mr. World.')

if __name__=='__main__':
  helloworld()