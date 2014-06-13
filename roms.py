# -*- coding: utf-8 -*-
#
# Description: Loose scripts and functions to manipulate ROMS fileds.
#              works with the 'romslab' (https://github.com/rsoutelino/romslab)
#              and 'pyroms' (https://github.com/kshedstrom/pyroms) modules.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

__all__ = ['maxvel_ke',
           'make_flat_ini']

import numpy as np
from scipy.interpolate import interp1d
from romslab import RomsGrid, RunSetup
from netCDF4 import Dataset
from pyroms.vgrid import s_coordinate, s_coordinate_2, s_coordinate_4

def vel_ke(avgfile, verbose=False):
	"""
	USAGE
	-----
	t, avgvel, maxvel, KEavg = vel_ke(avgfile, verbose=False)

	Calculates domain-averaged kinetic energy, mean and maximum velocity
	for each time record of a ROMS *.avg or *.his file.
	"""
	nc = Dataset(avgfile)

	try:
		U = nc.variables['u']
		V = nc.variables['v']
		uvrho = False
	except KeyError:
		U = nc.variables['u_eastward']
		V = nc.variables['v_northward']
		uvrho = True

	nt = U.shape[0]
	t = nc.variables['ocean_time'][:]
	t = t - t[0]

	avgvel = np.array([])
	maxvel = np.array([])
	KEavg = np.array([])
	for ti in xrange(nt):
		tp = ti + 1
		print "Processing time record %s of %s"%(tp,nt)
		uu = U[ti,:]
		vv = V[ti,:]

		if not uvrho:
			# Calculate u and v at PSI-points.
			u = 0.5*(uu[:,1:,:] + uu[:,:-1,:])
			v = 0.5*(vv[:,:,1:] + vv[:,:,:-1])
		else:
			# U and V both at RHO-points, no reshaping necessary.
			u = uu
			v = vv
			pass

		# Maximum velocity and domain-averaged kinetic energy.
		u2 = u.ravel()**2
		v2 = v.ravel()**2
		mag = np.sqrt(u2+v2)
		magavg = mag.mean()
		magmax = mag.max()

		ke = 0.5*(u2 + v2)
		ke = ke.mean()

		if verbose:
			print "meanvel, maxvel, avgKE = %.2f m/s, %.2f m/s, %f m2/s2"%(magavg,magmax,ke)

		KEavg = np.append(KEavg,ke)
		avgvel = np.append(avgvel,magavg)
		maxvel = np.append(maxvel,magmax)

	return t, avgvel, maxvel, KEavg

def make_flat_ini(grdfile, setup_file, profile_z, profile_T, profile_S, return_all=False):
	"""
	USAGE
	-----
	z_rho, temp, salt = make_flat_ini(grdfile, setup_file, profile_z, profile_T, profile_S, return_all=False)

	or

	z_rho, temp, salt, u, v, ubar, vbar, zeta = make_flat_ini(grdfile, setup_file, profile_z, profile_T, profile_S, return_all=True)

	-----
	Creates a ROMS initial conditions file with flat stratification. Requires a ROMS grid file, a ROMS
	setup file (read by the romslab.RunSetup function) and a single T,S profile that is repeated across the domain.

	INPUTS
	------
	grdfile:    A string with the path to the ROMS grid file.
	setup_file: A string with the path to the ROMS setup file.
	profile_z:  The depth coordinates of the T,S profiles (positive and increasing downward).
	profile_T:  The source temperature profile used to build the flat stratification field.
	profile_S:  The source salinity profile used to build the flat stratification field.

	OUTPUTS
	-------
	z_rho:                              Array with The depths of the RHO-points of the ROMS grid
										(negative and decreasing downward).
	temp, salt, u, v, ubar, vbar, zeta: The ROMS initial fields. temp and salt are have flat isolines,
	                                    and u, v, ubar, vbar and zeta are zero everywhere.
	"""
	profile_z, profile_T, profile_S = map(np.asanyarray, (profile_z, profile_T, profile_S))

	# Getting grid parameters from *.setup file.
	Run = RunSetup(setup_file)
	N = np.int(Run.klevels)
	Run.vtransform = np.int(Run.vtransform)
	Run.vstretching = np.int(Run.vstretching)

	# Read grid file.
	grd = RomsGrid(grdfile)

	hraw = grd.ncfile.variables['hraw'][:]
	h = grd.ncfile.variables['h'][:]

	# get the vertical ROMS grid.
	if Run.vstretching == 1:
	    scoord = s_coordinate(h, Run.theta_b, Run.theta_s, Run.tcline, N, hraw=hraw, zeta=None)
	elif Run.vstretching == 2:
	    scoord = s_coordinate_2(h, Run.theta_b, Run.theta_s, Run.tcline, N, hraw=hraw, zeta=None)
	elif Run.vstretching == 4:
	    scoord = s_coordinate_4(h, Run.theta_b, Run.theta_s, Run.tcline, N, hraw=hraw, zeta=None)
	else:
	    raise Warning, 'Unknow vertical stretching: Vtrans = %s'%Run.vstretching

	# Array of depths of RHO-points.
	z_rho = scoord.z_r[:]
	zroms = -np.flipud(z_rho)

	# Initializing temp and salt arrays.
	etamax, ximax = h.shape
	rhodim3 = (N, etamax, ximax)
	temp = np.empty(rhodim3)
	salt = np.empty(rhodim3)

	for j in xrange(etamax):
		for i in xrange(ximax):

			roms_z = zroms[:,j,i]

			# If at some point in the ROMS grid the deepest RHO-point lies
			# deeper than the deepest available value in the source T,S profiles,
			# extrapolate that deeper value to the deepest RHO-point depth before
			# interpolating.
			if roms_z[-1] > profile_z[-1]:
				profile_z = np.append(profile_z, roms_z[-1])
				profile_T = np.append(profile_T, profile_T[-1])
				profile_S = np.append(profile_S, profile_S[-1])

			# Vertically interpolate source temperature profile.
			I = interp1d(profile_z, profile_T, kind='linear', bounds_error=True)
			roms_T = I(roms_z)
			temp[:,j,i] = np.flipud(roms_T) # ATTENTION: Flipping array according to the
											# orientation of the sigma axis [-0.95 ... -0.05]

			# Vertically interpolate source salinity profile.
			I = interp1d(profile_z, profile_S, kind='linear', bounds_error=True)
			roms_S = I(roms_z)
			salt[:,j,i] = np.flipud(roms_S) # ATTENTION: Flipping array according to the
											# orientation of the sigma axis [-0.95 ... -0.05]

	# Creating free-surface elevation and velocity fields at rest.
	rhodim2 = (etamax, ximax)
	udim2 = (etamax, ximax-1)
	vdim2 = (etamax-1, ximax)
	udim3 = (N, etamax, ximax-1)
	vdim3 = (N, etamax-1, ximax)

	zeta = np.zeros(rhodim2)
	ubar = np.zeros(udim2)
	vbar = np.zeros(vdim2)
	u = np.zeros(udim3)
	v = np.zeros(vdim3)

	# Adding time axis to all arrays.
	zeta,ubar,vbar,u,v,temp,salt = map(np.expand_dims, (zeta,ubar,vbar,u,v,temp,salt), [0]*7)

	if return_all:
		return z_rho, temp, salt, u, v, ubar, vbar, zeta
	else:
		return z_rho, temp, salt









if __name__=='__main__':
  import doctest
  doctest.testmod()