# -*- coding: utf-8 -*-
#
# Description: Coursework-related scripts.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

__all__ = []

import numpy as np
import matplotlib.pyplot as plt

def dynmodes(n=6,lat0=5.,H=[80,170,175,250,325,3000],S=[24.97,26.30,26.83,27.12,27.32,27.77]):
	"""Computes the discrete eigenmodes (dynamical modes)
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
	rho0 = 1027.                        # [kg m^{-3}]

if __name__=='__main__':
	print None