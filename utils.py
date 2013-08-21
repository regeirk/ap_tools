# -*- coding: utf-8 -*-
#
# Description: Loose scripts and functions of personal use.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta

def goes_subset(time=datetime(2010,1,2), data_dir='/home/andre/.goes_data/', dt=24, cloud_thresh=2.0):
	"""Downloads a subset of GOES SST data from the PODAAC FTP.
	ftp://podaac-ftp.jpl.nasa.gov/allData/goes/L3/goes_6km_nrt/americas/."""

	if type(time)!=list:
		time = [time]

    # Getting files for each day.
	for date in time:
		yyyy = str(date.year)
		dd = str(date.day).zfill(3)
		head = 'ftp://podaac-ftp.jpl.nasa.gov/OceanTemperature/goes/L3/goes_6km_nrt/americas/%s/%s/' %(yyyy,dd)
		filename = 'sst%s?_%s_%s' % (str(dt),yyyy,dd) # Interval can be 1, 3 or 24 (hourly, 3-hourly or daily).
		url = head + filename
		cmd = "wget --tries=inf %s" %url
		original_dir = os.getcwd()
		if os.path.isdir(data_dir):
			os.chdir(data_dir)
		else:
			os.makedirs(data_dir)
			os.chdir(data_dir)
		# Downloading files.
		os.system(cmd)
		# read_goes.py()
		# cmd = "mv %s" %filename
	os.chdir(original_dir)

	return None

class GOES(object):
	"""Container class for GOES SST data."""
	def __init__(self):
		return None		

def helloworld():
  np.disp('Hello there, Mr. World.')

if __name__=='__main__':
  helloworld()