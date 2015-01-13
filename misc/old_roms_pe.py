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
