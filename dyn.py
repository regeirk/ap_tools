# -*- coding: utf-8 -*-
#
# Description: Functions to calculate common dynamical quantities.
#
# Author:      André Palóczy Filho
# E-mail:      paloczy@gmail.com

from __future__ import division

__all__ = ['pgrd']

import numpy as np

def pgrd():
	"""
	USAGE
	-----
	P = pgrd(x, y, z, eta, rho)

	"""
	x,y,z,eta,rho = map(np.asanyarray,(x,y,z,eta,rho)

	return 1

if __name__=='__main__':
  import doctest
  doctest.testmod()