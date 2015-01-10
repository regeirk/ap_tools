# -*- coding: utf-8 -*-
#
# Description: Loose scripts and functions to manipulate ROMS fileds.
#              works with the 'romslab' (https://github.com/rsoutelino/romslab)
#              and 'pyroms' (https://github.com/kshedstrom/pyroms) modules.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

__all__ = ['vel_ke',
           'ape',
           'make_flat_ini']

import numpy as np
from scipy.interpolate import interp1d
try:
	from seawater import pres, pden, grav
except:
	try:
		from seawater.csiro import pres, pden, grav
	except:
		pass		
# from gsw import SA_from_SP, CT_from_pt, rho, grav

try:
	from romslab import RomsGrid, RunSetup
except:
	try:
		from romslab.romslab import RomsGrid, RunSetup
	except:
		pass

from netCDF4 import Dataset
import pyroms
from pyroms.vgrid import s_coordinate, s_coordinate_2, s_coordinate_4

def vel_ke(avgfile, verbose=False):
	"""
	USAGE
	-----
	t, avgvel2, maxvel2, KEavg2, avgvel3, maxvel3, KEavg3 = vel_ke(avgfile, verbose=False)

	Calculates barotropic and baroclinic domain-averaged kinetic energies
	and mean/maximum velocities for each time record of a ROMS *.avg or *.his file.
	"""
	nc = Dataset(avgfile)

	try:
		Ubar = nc.variables['ubar']
		Vbar = nc.variables['vbar']
		uvrho2 = False
	except KeyError:
		Ubar = nc.variables['ubar_eastward']
		Vbar = nc.variables['vbar_northward']
		uvrho2 = True

	try:
		U = nc.variables['u']
		V = nc.variables['v']
		uvrho3 = False
	except KeyError:
		U = nc.variables['u_eastward']
		V = nc.variables['v_northward']
		uvrho3 = True

	t = nc.variables['ocean_time'][:]
	t = t - t[0]
	nt = t.size

	avgvel2 = np.array([])
	maxvel2 = np.array([])
	avgvel3 = np.array([])
	maxvel3 = np.array([])
	KEavg2 = np.array([])
	KEavg3 = np.array([])
	for ti in xrange(nt):
		tp = ti + 1
		print "Processing time record %s of %s"%(tp,nt)
		uubar = Ubar[ti,:]
		vvbar = Vbar[ti,:]
		uu = U[ti,:]
		vv = V[ti,:]

		if not uvrho2:
			# Calculate ubar and vbar at PSI-points.
			ubar = 0.5*(uubar[1:,:] + uubar[:-1,:])
			vbar = 0.5*(vvbar[:,1:] + vvbar[:,:-1])
		else:
			# Ubar and Vbar both at RHO-points, no reshaping necessary.
			ubar = uubar
			vbar = vvbar

		if not uvrho3:
			# Calculate u and v at PSI-points.
			u = 0.5*(uu[:,1:,:] + uu[:,:-1,:])
			v = 0.5*(vv[:,:,1:] + vv[:,:,:-1])
		else:
			# U and V both at RHO-points, no reshaping necessary.
			u = uu
			v = vv

		# Mean and maximum barotropic/baroclinic velocities and domain-averaged barotropic/baroclinic kinetic energy.
		ubar2 = ubar.ravel()**2
		vbar2 = vbar.ravel()**2
		mag2 = np.sqrt(ubar2+vbar2)		
		u2 = u.ravel()**2
		v2 = v.ravel()**2
		mag3 = np.sqrt(u2+v2)
		magavg2 = mag2.mean()
		magmax2 = mag2.max()
		magavg3 = mag3.mean()
		magmax3 = mag3.max()

		ke2 = 0.5*(ubar2 + vbar2)
		ke2 = ke2.mean()
		ke3 = 0.5*(u2 + v2)
		ke3 = ke3.mean()

		if verbose:
			print "meanvel2D, maxvel2D, avgKE2D = %.2f m/s, %.2f m/s, %f m2/s2"%(magavg2,magmax2,ke2)
			print "meanvel3D, maxvel3D, avgKE3D = %.2f m/s, %.2f m/s, %f m2/s2"%(magavg3,magmax3,ke3)

		avgvel2 = np.append(avgvel2,magavg2)
		maxvel2 = np.append(maxvel2,magmax2)
		KEavg2 = np.append(KEavg2,ke2)
		avgvel3 = np.append(avgvel3,magavg3)
		maxvel3 = np.append(maxvel3,magmax3)
		KEavg3 = np.append(KEavg3,ke3)

	return t, avgvel2, maxvel2, KEavg2, avgvel3, maxvel3, KEavg3

def pe(avgfile, grdfile, gridid=None, maskfile='/media/Armadillo/bkp/lado/MSc/work/ROMS/plot_outputs3/msk_shelf.npy', normalize=False, verbose=True):
	"""
	USAGE
	-----
	t, pe = pe(avgfile, grdfile, gridid=None, maskfile='/media/Armadillo/bkp/lado/MSc/work/ROMS/plot_outputs3/msk_shelf.npy', normalize=False, verbose=True):

	Calculates Potential Energy (PE) change integrated within a control volume
	for each time record of a ROMS *.avg or *.his file. The change is computed relative to
	the initial conditions, i.e., rhop(x,y,z,t=ti) = rho(x,y,z,t=ti) - rho(x,y,z,t=t0).

                                          [-g*(rhop^2)]
	PE = Integrated in a control volume V [-----------]     # [J]
                                          [ 2*drhodz  ]

	If 'normalize' is set to 'True', then PE/V (mean PE density [J/m3]) is returned instead.

	Reference:
	----------
	Cushman-Roisin (1994): Introduction to Geophysical Fluid Dynamics, page 213,
	Combination of Equations 15-29 and 15-30.
	"""
	print "Loading outputs and grid."

	## Get outputs.
	avg = Dataset(avgfile)

	## Load domain mask.
	if maskfile:
		mask = np.load(maskfile)
		if type(mask[0,0])==np.bool_:
			pass
		else:
			mask=mask==1.

	## Getting mask indices.
	ilon,ilat = np.where(mask)

	## Get grid, time-dependent free-surface and topography.
	zeta = avg.variables['zeta']
	grd = pyroms.grid.get_ROMS_grid(gridid, zeta=zeta, hist_file=avgfile, grid_file=grdfile)

	## Get time.
	t = avg.variables['ocean_time'][:]
	t = t - t[0]
	nt = t.size

	## Get grid coordinates at RHO-points.
	lonr, latr = avg.variables['lon_rho'][:], avg.variables['lat_rho'][:]

	## Get grid spacings at RHO-points.
	## Find cell widths.
	dx = grd.hgrid.dx             # Cell width in the XI-direction.
	dy = grd.hgrid.dy             # Cell width in the ETA-direction.
	if maskfile:
		dA = dx[mask]*dy[mask]
	else:
		dA = dx*dy

	## Get temp, salt.
	temp = avg.variables['temp']
	salt = avg.variables['salt']

	## Find cell heights (at ti=0).
	zw = grd.vgrid.z_w[0,:]             # Cell depths (at ti=0).
	if maskfile:
		dz = zw[1:,ilat,ilon] - zw[:-1,ilat,ilon] # Cell height.
	else:
		dz = zw[1:,:] - zw[:-1,:]
	dz = 0.5*(dz[1:,:] + dz[:-1,:])               # Cell heights at W-points.

	## Get pres, g and pden (at ti=0).
	p0 = -zw # Approximation, for computational efficiency.
	p0 = 0.5*(p0[1:,:]+p0[:-1,:])

	if maskfile:
		rho0 = pden(salt[0,:,ilat,ilon],temp[0,:,ilat,ilon],p0[:,ilat,ilon],pr=0.)
	else:
		rho0 = pden(salt[0,:],temp[0,:],p0,pr=0.)

	if maskfile:
		g = grav(latr[mask])
	else:
		g = grav(latr)

	drho0 = rho0[1:,:] - rho0[:-1,:]
	rho0z = drho0/dz # Background potential density vertical gradient.

	PE = np.array([])
	for ti in xrange(nt):
		tp = ti + 1
		print "Processing time record %s of %s"%(tp,nt)

		if maskfile:
			rhoi = pden(salt[ti,:,ilat,ilon],temp[ti,:,ilat,ilon],p0[:,ilat,ilon],pr=0.)
		else:
			rhoi = pden(salt[ti,:],temp[ti,:],p0,pr=0.)

		rhop = rhoi - rho0                                  # Density anomaly, i.e., rho(x,y,z,t=ti) - rho(x,y,z,t=0)
		rhop = 0.5*(rhop[1:,:] + rhop[:-1,:])

		## Find cell heights.
		zw = grd.vgrid.z_w[ti,:]                      # Cell depths (at ti=0).
		if maskfile:
			dz = zw[1:,ilat,ilon] - zw[:-1,ilat,ilon] # Cell height.
		else:
			dz = zw[1:,:] - zw[:-1,:]

		## Find cell volumes.
		print dx.shape,dy.shape,dz.shape
		dV = dA*dz # [m3]
		dV = 0.5*(dV[1:,:]+dV[:-1,:])

		## Gravitational Available Potential Energy density (energy/volume).
		print g.shape
		print rhop.shape
		print rho0z.shape
		pe = -g*(rhop**2)/(2*rho0z) # [J/m3]

		## Do volume integral to calculate Gravitational Available Potential Energy of the control volume.
		Pe = np.sum(pe*dV) # [J]

		if normalize:
			V = dV.sum()
			Pe = Pe/V
			print ""
			print "Total volume of the control volume is %e m3."%V
			print "Normalizing PE by this volume, i.e., mean PE density [J/m3]."
			print ""

		if verbose:
			if normalize:
				print "PE = %e J/m3"%Pe
			else:
				print "PE = %e J"%Pe

		PE = np.append(PE, Pe)

	return t, PE

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
