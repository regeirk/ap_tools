# -*- coding: utf-8 -*-
#
# Description: Functions to calculate common dynamical quantities.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

from __future__ import division

__all__ = ['pgf']

import numpy as np
from gsw import grav

def pgf(z, y, x, eta, rho, pa=0., rho0=1025., geographic=True):
	"""
	USAGE
	-----
	Py, Px = pgf(z, y, x, eta, rho, pa=0., rho0=1025., geographic=True)

	Calculates total horizontal pressure gradient force per unit mass [m/s2]
	(barotropic + baroclinic + barometric components), i.e.,

	P(z,y,x) = -(1/rho0)*grad_H(pa) -g*grad_H(eta) + (g/rho0)*Integral{grad_H(rho)}.

	'rho0' is the reference potential density used in the Boussinesq Approximation
	(defaults to 1025. kg/m3), 'g' is the gravitational acceleration, 'pa(y,x)' is
	the atmospheric pressure in [N/m2] (defults to 0.), 'eta(y,x)' is the free surface elevation in [m],
	'rho(z,y,x)' is the potential density and 'grad_H' is the horizontal gradient operator.
	The Integral(rho) is calculated from z' = eta(x,y) down through z' = z.

	The coordinate arrays (z,y,x) are distances in the (vertical,meridional,zonal)
	directions. The vertical axis originates at the surface (z = 0), i.e.,
	rho[0,y,x] = rho_surface and
	rho[-1,y,x] = rho_bottom.

	If geographic==True (default), (y,x) are assumed to be
	(latitude,longitude) and are converted to meters before
	computing (dy,dx). If geographic==False, (y,x) are assumed to be in meters.
	"""
	z, y, x, eta, rho = map(np.asanyarray, (z, y, x, eta, rho))

	ny, nx = eta.shape                       # Shape of the (x,y,u,v) arrays.
	if z.ndim==1:
		z = np.expand_dims(z, 1)
		z = np.expand_dims(z, 1)
		z = np.tile(z, (1, ny, nx))

	## Calculating grid spacings.
	if geographic:
		dlat, _ = np.gradient(y)
		_, dlon = np.gradient(x)
		deg2m = 111120.0                     # [m/deg]
		dx = dlon*deg2m*np.cos(y*np.pi/180.) # [m]
		dy = dlat*deg2m                      # [m]
	else:
		dy, _ = np.gradient(y)
		_, dx = np.gradient(x)

	dz, _, _ = np.gradient(z)
	dz = np.abs(dz)

	## Get gravitational acceleration.
	if geographic:
		g = grav(y)
	else:
		g = 9.81

	## pa (x,y) derivatives.
	if pa:
		dpay, dpax = np.gradient(pa)
		dpady = dpay/dy
		dpadx = dpax/dx

	## eta (x,y) derivatives.
	detay, detax = np.gradient(eta)
	detady = detay/dy
	detadx = detax/dx

	## rho (x,y) derivatives.
	_, drhoy, drhox = np.gradient(rho)
	drhody = drhoy/dy
	drhodx = drhox/dx

	## Barometric pressure gradient force per unit mass.
	if pa==0.:
		PGF_bm_y, PGF_bm_x = np.zeros((ny, nx)), np.zeros((ny, nx))
	else:
		PGF_bm_y = -dpady/rho0
		PGF_bm_x = -dpadx/rho0

	## Barotropic pressure gradient force per unit mass.
	PGF_bt_y = -g*detady
	PGF_bt_x = -g*detadx

	## Vertical integration from z' = eta(x,y) through z' = z.
	Iy = np.cumsum(drhody*dz)
	Ix = np.cumsum(drhodx*dz)

	## Baroclinic pressure gradient force per unit mass.
	PGF_bc_y = +g*Iy/rho0
	PGF_bc_x = +g*Ix/rho0

	## Total pressure gradient force per unit mass.
	PGF_y = PGF_bm_y + PGF_bt_y + PGF_bc_y
	PGF_x = PGF_bm_x + PGF_bt_x + PGF_bc_x

	return PGF_y, PGF_x

if __name__=='__main__':
  import doctest
  doctest.testmod()