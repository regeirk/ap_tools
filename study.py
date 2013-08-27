# -*- coding: utf-8 -*-
#
# Description: Coursework-related scripts.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

__all__ = ['dynmodes']

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from romslab import near
from oceans.plotting import rstyle

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

def dynmodes(n=6, lat0=5., plot=False, model='Fratantoni_etal1995'):
	"""
	Computes the discrete eigenmodes (dynamical modes)
	for a quasi-geostrophic ocean with n isopycnal layers.
	Rigid lids are placed at the surface and the bottom.

	Inputs:
	-------
	n:    Number of layers.
	lat0: Reference latitude.
	H:    List of rest depths of each layer.
    S:    List of potential density values for each layer.
	"""
	omega = 2*np.pi/86400.              # [rad s^{-1}]
	f0 = 2*omega*np.sin(lat0*np.pi/180) # [s^{-1}]
	f02 = f0**2                         # [s^{-2}]
	g = 9.81                            # [m s^{-1}]
	rho0 = 27.                          # [kg m^{-3}]

	if model=='Fratantoni_etal1995':
		H = np.array([80.,170.,175.,250.,325.,3000.])
		S = np.array([24.97,26.30,26.83,27.12,27.32,27.77])
		tit = ur'Modelo de seis camadas para a CNB de Fratantoni \textit{et al.} (1995)'
		figname = 'vmodes_fratantoni_etal1995'
	elif model=='Bub_Brown1996':
		H = np.array([150.,440.,240.,445.,225.,2500.])
		S = np.array([24.13,26.97,27.28,27.48,27.74,27.87])
		tit = ur'Modelo de seis camadas para a CNB de Bub e Brown (1996)'
		figname = 'vmodes_bub_brown1996'

	# Normalized density jumps.
	E = (S[1:]-S[:-1])/rho0
	# Rigid lids at the surface and the bottom,
	# meaning infinite density jumps.
	E = np.hstack( (np.inf,E,np.inf) )

	# Building the tridiagonal matrix.
	A = np.zeros( (n,n) )
	for i in xrange(n):
		A[i,i] = -f02/(E[i+1]*g*H[i]) -f02/(E[i]*g*H[i]) # The main diagonal.
		if i>0:
			A[i,i-1] = f02/(E[i]*g*H[i])
		if i<(n-1):
			A[i,i+1] = f02/(E[i+1]*g*H[i])

	# get eigenvalues and convert them
	# to internal deformation radii
	lam,v = eig(A)
	lam = np.abs(lam)

	# Baroclinic def. radii in km:
	uno = np.ones( (lam.size,lam.size) )
	Rd = 1e-3*uno/np.sqrt(lam)
	Rd = np.unique(Rd)

	np.disp("Deformation radii [km]:")
	[np.disp(int(r)) for r in Rd]

	# orthonormalize eigenvectors, i.e.,
	# find the dynamical mode vertical structure functions
	F = np.zeros( (n,n) )

	for i in xrange(n):
		mi = v[:,i] # The vertical structure of the i-th vertical mode.
		fac = np.sqrt(np.sum(H*mi*mi)/np.sum(H))
		F[:,i] = 1/fac*mi

	F=-F
	F[:,0] = np.abs(F[:,1])
	z = np.linspace(0,np.sum(H),1000)

	nz = z.size
	Fi = np.zeros( (nz,n) )
	for i in xrange(n):
		for j in xrange(nz):
			idx = near(H,z[j])[0][0]
			Fi[j,i] = F[idx,i]

	z=-z
	# Plot the vertical modes.
	if plot:
		plt.close('all')
		kw = dict(fontsize=15, fontweight='black')
		fig,ax = plt.subplots()
		ax.hold(True)
		for i in xrange(n):
			l = 'Modo %s'%str(i+1)
			ax.plot(Fi[:,i], z, label=l)
		xl,xr = ax.set_xlim(-5,5)
		H2 = -np.cumsum(H)
		ax.hlines(H2,xl,xr,linestyle='dashed')
		ax.hold(False)
		ax.set_title(tit, **kw)
		ax.set_xlabel(ur'Autofunção [adimensional]', **kw)
		ax.set_ylabel(ur'Profundidade [m]', **kw)
		ax.legend(loc='best')
		rstyle(ax)
		fmt='png'
		fig.savefig('/home/andre/'+figname+'.'+fmt, format=fmt, bbox='tight')

if __name__=='__main__':
	print None