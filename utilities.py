## Utilities for wake model

import numpy as np
import time
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import wakeModels
import copy
import time
import numba as nb
import pandas as pd

#==================================================================FUNCTIONS============================================================================
# 1. visualizeHorizontal - plot the freestream flow field at hub height 
# 2. visualizeLidar - plot the output of the Stuttgart lidar model implemented
# 3. rotatedCoordinates - rotate the turbines to 270 degrees and evaluate the wake model in this way.  Makes coding the wake models simpler
# 4. avgVelocity - compute the average velocity across a rotor
# 5. combineWakes - combine the wakes based linear superposition or sum of squares
# 6. computeAddedTI - compute the turbulence intensity added to the wake based on turbine operation and ambient turbulence intensity
# 7. computeOverlap - compute overlap of an upstream turbine wake (based on a specified threshold)
# 8. computePower - compute the power at each turbine based on the effective velocity determined using the wake model
# 9. initializeFlowField - initialize the flow field used in the 3D model based on shear using the power log law
# 10. determineCpCt - interpolation function based on Cp, Ct, and wind speed
# 11. outputUpts - compute the output velocity at the points specified in inputData
# ======================================================================================================================================================

def visualizeHorizontal(xLen,yLen,zLen,Ufield,inputData):

	# plot the freestream flow field at hub height 

	# future extension (bold lines based on axial induction)
	print('Plotting Flow Field Horizontal...\n')

	Uinf 			= inputData['windSpeed']
	windDirection 	= inputData['windDirection']
	tilt 			= inputData['tilt']
	D 				= inputData['rotorDiameter']
	yaw 			= inputData['yawAngles']
	xTurb 			= inputData['turbineX']
	yTurb 			= inputData['turbineY']
	zTurb 			= inputData['turbineZ']

	# number of turbines
	nTurbs = len(xTurb)

	# rotation angle for the turbines the wind turbines
	RotAng = -np.radians(windDirection -270.)

	# plot horizontal flow field
	plt.figure(figsize=(30,10))
	plt.contourf(xLen,yLen,Ufield[0,:,:],50,cmap='coolwarm',vmin=3.0,vmax=Uinf)
	cb = plt.colorbar()
	cb.ax.tick_params(labelsize=15)

	if inputData['Lidar']:

		yLidar = inputData['yLidar']
		zLidar = inputData['zLidar']
		xLidar = inputData['xLidar']

		numLocs = len(np.unique(xLidar))
		xLocs = np.zeros(numLocs+1)
		ymin = np.zeros(numLocs+1)
		ymax = np.zeros(numLocs+1)
		cols = yLidar.columns
		ymin[0] = yTurb[inputData['turbineLidar']]
		ymax[0] = yTurb[inputData['turbineLidar']]
		xLocs[0] = xTurb[inputData['turbineLidar']]
		for i in range(1,numLocs+1):
			xLocs[i] = xTurb[inputData['turbineLidar']] + np.unique(xLidar)[i-1]
			ymin[i] = yTurb[inputData['turbineLidar']] + np.min(yLidar[cols[i-1]])
			ymax[i] = yTurb[inputData['turbineLidar']] + np.max(yLidar[cols[i-1]])
			plt.plot([xLocs[i],xLocs[i]],[ymin[i],ymax[i]],'k')

		plt.plot(xLocs,ymin,'k',linewidth=3)
		plt.plot(xLocs,ymax,'k',linewidth=3)
		plt.plot([xTurb[inputData['turbineLidar']],xTurb[inputData['turbineLidar']] + 700],[yTurb[inputData['turbineLidar']],yTurb[inputData['turbineLidar']]],'k--',linewidth=3)

	for idx in range(nTurbs):

		if zTurb[idx] == zTurb[0]:

			yawR = np.radians(yaw[idx])

			# x coordinate
			xRot1 = xTurb[idx] - (-(D[idx]/2))*np.sin(RotAng+yawR)
			xRot2 = xTurb[idx] - ((D[idx]/2))*np.sin(RotAng+yawR)

			# z coordinate
			yRot1 = yTurb[idx] + (-(D[idx]/2))*np.cos(RotAng+yawR)
			yRot2 = yTurb[idx] + ((D[idx]/2))*np.cos(RotAng+yawR)

			plt.plot([xRot1,xRot2],[yRot1,yRot2],'k',linewidth=3)

	plt.axis('equal')
	plt.xlabel('x(m)',fontsize=15)
	plt.ylabel('y(m)',fontsize=15)
	plt.tick_params(which='both',labelsize=15)
	strTitle = 'Horizontal'
	plt.title('Horizontal',fontsize=15)

def visualizeLidar(xTurb,yTurb,X,Y,Z,Ufield,inputData):

	# plot the output of the Stuttgart lidar model implemented

	Uinf = inputData['windSpeed']
	
	# parameters from inputData
	idxTurb = inputData['turbineLidar']
	zTurb = inputData['turbineZ'][idxTurb]

	# Lidar parameters
	xLocs = np.unique(inputData['xLidar']) + xTurb
	numLocs = len(xLocs)

	# plot Lidar panels
	plt.figure(figsize=(20,3))
	nSamplesX = X.shape[2]
	nSamplesY = Y.shape[1]
	nSamplesZ = Z.shape[0]

	cols = inputData['yLidar'].columns

	for kk in range(numLocs):
		plt.subplot(1,numLocs,kk+1)
		yPlot = np.unique(inputData['yLidar'][cols[kk]]) + yTurb
		zPlot = np.unique(inputData['zLidar'][cols[kk]]) + zTurb

		Uplot = np.zeros((len(yPlot),len(zPlot)))

		# find the Ufield points that correspond to yPlot and zPlot
		for ii in range(len(zPlot)):
			for jj in range(len(yPlot)):
				if X[ii,jj,kk] == xLocs[kk] and Y[ii,jj,kk] == yPlot[jj] and Z[ii,jj,kk] == zPlot[ii]:
					Uplot[ii,jj] = Ufield[ii,jj,kk]

		plt.contourf(-yPlot,zPlot,Uplot,50,cmap='coolwarm',vmin=2.0,vmax=Uinf)
		plt.xlim([-60-yTurb,60-yTurb])
		plt.ylim([40,160])
		plt.colorbar()

def rotatedCoordinates(inputData,X,Y,Z):

	# this function rotates the coordinates so that the flow field is now in the frame of reference of the 270 degree wind direction.
	# this makes computing wakes and wake overlap much simpler

	windDirection = inputData['windDirection']
	xTurb = inputData['turbineX']
	yTurb = inputData['turbineY']
	nTurbs = len(xTurb)

	RotAng 	= np.radians(windDirection - 270.)          # rotate the wind frame of reference for easier computation
	xCenter = np.mean(xTurb)                       		# find the center of the xy-domain, i.e. the mean
	yCenter = np.mean(yTurb)                       		# find the center of the xy-domain, i.e. the mean

	if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
	    # initialize output vectors
		Xrot 	= np.zeros(X.shape)
		Yrot 	= np.zeros(Y.shape)		

		# number of samples in the x and y coordinates
		nSamplesX = X.shape[2]
		nSamplesY = Y.shape[1]
		nSamplesZ = Z.shape[0]

		# demean the X and Y grids and the x and y turbine coordinates in order to rotate
		X 		= X - xCenter
		Y 		= Y - yCenter
	
		# rotate the grid
		for i in range(nSamplesZ):
			for j in range(nSamplesY):
				for k in range(nSamplesX):
					Xrot[i,j,k] = X[i,j,k]*np.cos(RotAng) - Y[i,j,k]*np.sin(RotAng) 
					Yrot[i,j,k] = X[i,j,k]*np.sin(RotAng) + Y[i,j,k]*np.cos(RotAng)
					#print(Xrot[i,j,k],Yrot[i,j,k],X[i,j,k],Y[i,j,k])

	xTurb 	= xTurb - xCenter
	yTurb 	= yTurb - yCenter
	xTurbRot = np.zeros(nTurbs)
	yTurbRot = np.zeros(nTurbs)

	# rotate the turbine coordinates
	for i in range(nTurbs):
		xTurbRot[i] = xTurb[i]*np.cos(RotAng) - yTurb[i]*np.sin(RotAng)
		yTurbRot[i] = xTurb[i]*np.sin(RotAng) + yTurb[i]*np.cos(RotAng)

	# put the mean back in and return the X, Y domains and the rotated turbine locations
	if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
		return Xrot+xCenter,Yrot+yCenter,xTurbRot+xCenter,yTurbRot+yCenter
	else:
		return xTurbRot+xCenter,yTurbRot+yCenter

def avgVelocity(X,Y,Z,Ufield,xTurb,yTurb,zTurb,D,rotorPts,inputData,turbI):

	# this function averages the velocity across the rotor.  The mean wind speed across the rotor is used.  This may be updated in the future.

	tilt = inputData['tilt'][turbI]
	yaw  = inputData['yawAngles'][turbI]

	xCenter = xTurb
	zR = (D/2.)*np.cos(np.radians(tilt))
	yR = (D/2.)*np.cos(np.radians(yaw))
	xR = (D/2.)*np.sin(np.radians(yaw)) + (D/2.)*np.sin(np.radians(tilt))
	yPts = np.linspace(yTurb-yR,yTurb+yR,rotorPts)
	zPts = np.linspace(zTurb-zR,zTurb+zR,rotorPts)

	# this function averages the velocity across the whole rotor
	# number of points in the X and Y domain
	nSamplesX = X.shape[2]
	nSamplesY = Y.shape[1]
	nSamplesZ = Z.shape[0]

	# define the rotor plane along with the points associated with the rotor
	Xtmp = []
	Ytmp = []
	Ztmp = []
	Utmp = []

	# if X and Y are large, it is essential to resample the domains to only include points that we care about
	# this significantly speeds up the process
	resampStart = time.time()
	for i in range(nSamplesZ):
		for j in range(nSamplesY):
			for k in range(nSamplesX):
				if X[i,j,k] >= (xTurb - D/4.) and X[i,j,k] <= (xTurb + D/4.):
					dist = np.sqrt( (Y[i,j,k]-yTurb)**2 + (Z[i,j,k]-zTurb)**2 )
					if dist <= (D/2.):
						Xtmp.append(X[i,j,k])
						Ytmp.append(Y[i,j,k])
						Ztmp.append(Z[i,j,k])
						Utmp.append(Ufield[i,j,k])
	resampEnd = time.time()

	# interpolate the points to find the average velocity across the rotor 
	interpStart = time.time()
	if len(Xtmp) == 0 or len(Ytmp) == 0 or len(Ztmp) == 0:
		print('Too few points for outputFlowField, please run again with more points')

	utmp = []

	for i in range(rotorPts):
		for j in range(rotorPts):
			dist = np.sqrt( (yPts[i]-yTurb)**2 + (zPts[j]-zTurb)**2 )
			if dist <= (D/2.):
				utmp.append(griddata((Xtmp,Ytmp,Ztmp),Utmp,(xCenter,yPts[i],zPts[j]), method='nearest'))

	interpEnd = time.time()

	return np.mean(utmp)

def combineWakes(inputData,Uinf,Ueff,Ufield,tmp):

	# this function allows for different wake superpositions

	# freestream linear superposition
	if inputData['combineWakes'] == 0:
		vel = Uinf - ((Uinf - Ufield) + tmp)

	# local velocity linear superposition
	elif inputData['combineWakes'] == 1:
		vel = Uinf - ((Ueff - Ufield) + tmp)

	# sum of squares freestream superposition
	elif inputData['combineWakes'] == 2:
		vel = Uinf - np.sqrt((Uinf - Ufield)**2 + tmp**2)

	# sum of squares local velocity superposition
	elif inputData['combineWakes'] == 3:
		vel = Ueff - np.sqrt((Ueff - Ufield)**2 + tmp**2)
	else:
		vel = Ufield

	return vel

def computeAddedTI(X,Y,Z,inputData,xTurb,yTurb,zTurb,turbI,turbineID):
	# this function is necessary for computing the TI that impacts the specified turbine
	# loop through all of the turbines and compute only the velocity at the rotor of interest
	# computes areaoverlap and TI contribution
	# the TI is addition is limited to a specified distance upstream (default = 15*D)

	D 					= inputData['rotorDiameter']
	Uinf 				= inputData['windSpeed']
	kd 					= inputData['wakeDeflection']
	ke 					= inputData['wakeExpansion']
	ad 					= inputData['ad']
	bd 					= inputData['bd']
	B 					= inputData['NumBlades']
	TSR 				= inputData['TSR']

	nTurbs 				= len(inputData['turbineX'])
	rotorPts 			= inputData['rotorPts']
	rotorPts 		= int(np.round(np.sqrt(rotorPts)))
	TI_0			    = inputData['turbulenceIntensity']
	TI_a 	= inputData['TIa']
	TI_b 	= inputData['TIb']
	TI_c 	= inputData['TIc']
	TI_d 	= inputData['TId']

	# turbine operation
	yaw 				= inputData['yawAngles']
	tilt 				= inputData['tilt']
	Ct 					= inputData['Ct']

	# define axial induction in terms of Ct
	a = np.zeros(nTurbs)
	for i in range(nTurbs):
		a[i] = (0.5/np.cos(yaw[i]*np.pi/180.))*(1-np.sqrt(1-Ct[i]*np.cos(yaw[i]*np.pi/180)))

	# generate compressed grid just for the one turbine (speed up things a bunch)
	Xt = [xTurb[turbI]]
	Yt = np.linspace(yTurb[turbI]-(D[turbI]/2), yTurb[turbI]+(D[turbI]/2), rotorPts)
	Zt = np.linspace(zTurb[turbI]-(D[turbI]/2), zTurb[turbI]+(D[turbI]/2), rotorPts)

	X = np.zeros((len(Zt),len(Yt),len(Xt)))
	Y = np.zeros((len(Zt),len(Yt),len(Xt)))
	Z = np.zeros((len(Zt),len(Yt),len(Xt)))

	# generate grid
	count = 0.
	for k in range(len(Xt)):
		Yt = np.linspace(yTurb[turbI]-(D[turbI]/2), yTurb[turbI]+(D[turbI]/2), rotorPts)
		Zt = np.linspace(zTurb[turbI]-(D[turbI]/2), zTurb[turbI]+(D[turbI]/2), rotorPts)
		for j in range(len(Yt)):
			for i in range(len(Zt)):
				X[i,j,k] = Xt[k]
				Y[i,j,k] = Yt[j]
				Z[i,j,k] = Zt[i]

	# determine the amount of the grid to loop through (depends on if outputting the velocity field)
	nSamplesX 			= X.shape[2]
	nSamplesY 			= Y.shape[1]
	nSamplesZ 			= Z.shape[0]

	# initialize flow field
	Ufield 				= initializeFlowField(X,Y,Z,inputData)
	UfieldOrig			= copy.copy(Ufield)

	# determine the effective velocity generated by that specific turbine
	AreaOverlap			= np.zeros(nTurbs)
	TI_calc 			= np.zeros(nTurbs)

	initWake 			= 0
	n1 = time.time()

	for turbIdx in turbineID:

		Ufield = UfieldOrig
		UfieldOrig			= copy.copy(Ufield)

		# compute the x and y offset due to yaw
		yR = (D[turbIdx]/2.)*np.cos(np.radians(yaw[turbIdx]))
		xR = (D[turbIdx]/2.)*np.sin(np.radians(yaw[turbIdx]))

		dist = np.sqrt( (xTurb[turbI]-xTurb[turbIdx])**2 + (yTurb[turbI]-yTurb[turbIdx])**2 )
		if dist > inputData['TIdistance']:
			continue

		# cycle through the previously defined X and Y domains
		for xIdx in range(nSamplesX):
			for yIdx in range(nSamplesY):
				for zIdx in range(nSamplesZ):

					if X[zIdx,yIdx,xIdx] >= 1.05*abs(xTurb[turbIdx]):

						# use GAUSS wake model
						if inputData['WakeModel'] == 2:
							c 		= wakeModels.GAUSS(X[zIdx,yIdx,xIdx],Y[zIdx,yIdx,xIdx],Z[zIdx,yIdx,xIdx],
													   xTurb[turbIdx],yTurb[turbIdx],zTurb[turbIdx],inputData,turbIdx,TI_0,UfieldOrig[zIdx,yIdx,xIdx])
							velDef 	= c

						# use FLORIS or Jensen wake model
						elif inputData['WakeModel'] == 0 or inputData['WakeModel'] == 1:	

							if (X[zIdx,yIdx,xIdx] > (xTurb[turbIdx]-abs(xR)) and 
							    Y[zIdx,yIdx,xIdx] > (yTurb[turbIdx] - 2*D[turbIdx]) and Y[zIdx,yIdx,xIdx] < (yTurb[turbIdx] + 2*D[turbIdx]) and
							    Z[zIdx,yIdx,xIdx] > (zTurb[turbIdx] - 2*D[turbIdx]) and Z[zIdx,yIdx,xIdx] < (zTurb[turbIdx] + 2*D[turbIdx])):
								yDisp = wakeModels.Jimenez(np.radians(yaw[turbIdx]),Ct[turbIdx],kd[turbIdx],X[zIdx,yIdx,xIdx]-xTurb[turbIdx],D[turbIdx],ad,bd)
								zDisp = wakeModels.Jimenez(np.radians(tilt[turbIdx]),Ct[turbIdx],kd[turbIdx],X[zIdx,yIdx,xIdx]-xTurb[turbIdx],D[turbIdx],ad,bd)
							else:
								yDisp = 0.0
								zDisp = 0.0
								continue

							# define the edges of the wake
							xTurbY = xTurb[turbIdx] + ( (Y[zIdx,yIdx,xIdx]-yTurb[turbIdx])*np.tan(-np.radians(yaw[turbIdx])))
							rWake = ke[turbI]*(X[zIdx,yIdx,xIdx]-(xTurb[turbI]-xR))
							rCenterY = xTurb[turbI] + yDisp
							rCenterZ = xTurb[turbI] + zDisp

							# compute the velocity deficit
							# define the edges of the wake
							rWake = ke[turbI]*(X[zIdx,yIdx,xIdx]-(xTurb[turbI]-xR)) + (D[turbI]/2)
							rCenterY = yTurb[turbI] + yDisp
							rCenterZ = zTurb[turbI] + zDisp

							rtmp = np.sqrt( (Y[zIdx,yIdx,xIdx] - rCenterY)**2 + (Z[zIdx,yIdx,xIdx] - rCenterZ)**2 )
							if (X[zIdx,yIdx,xIdx] >= xTurb[turbI] and rtmp <= rWake):

								# FLORIS model
								if inputData['WakeModel']==1:
									c = wakeModels.FLORIS(X[zIdx,yIdx,xIdx],Y[zIdx,yIdx,xIdx],Z[zIdx,yIdx,xIdx],yDisp,zDisp,xTurb[turbIdx],yTurb[turbIdx],zTurb[turbIdx],inputData,turbIdx)
									velDef = Uinf*2.*a[turbIdx]*c
								
								# Jensen model
								elif inputData['WakeModel']==0:
									c = wakeModels.Jensen(inputData,turbIdx,X[zIdx,yIdx,xIdx],xTurb[turbIdx])
									velDef = Uinf*2.*a[turbIdx]*c

							else:
								velDef = 0.0

						else:
							if initWake == 0:
								print('No wake model specified')
								velDef = 0.0
								initWake = 1
								break

						Ufield[zIdx,yIdx,xIdx] = UfieldOrig[zIdx,yIdx,xIdx] - velDef

		# compute percent overlap
		if turbIdx == turbI:
			AreaOverlap[turbIdx] = 0.0
		else:
			AreaOverlap[turbIdx] = computeOverlap(Ufield,UfieldOrig,inputData,Z)

		# length of near wake calucation based on Bastankah and Porte-Agel
		m = 1/np.sqrt(1-Ct[turbIdx])
		r0 = D[turbIdx]*np.sqrt((m+1)/2)
		drdx_a = 2.5*TI_0 + 0.005
		drdx_m = (1-m)*np.sqrt(1.49 + m)/(9.76*(1+m))
		drdx_lambda = 0.012*B[turbIdx]*TSR[turbIdx]
		drdx = np.sqrt(drdx_a**2 + drdx_m**2 + drdx_lambda**2)
		xn = (np.sqrt(0.214 + 0.144*m)*(1 - np.sqrt(0.134 + 0.124*m))*r0) / ((1-np.sqrt(0.214+0.144*m))*np.sqrt(0.134+0.124*m)*drdx)

		# Niayifar and Porte-Agel, "A new analytical model for wind farm power prediction", 2015
		if turbI == turbIdx or xTurb[turbI] < xTurb[turbIdx]:
			TI_calc[turbIdx] = 0.0

		else:
			if (xTurb[turbI]-xTurb[turbIdx]) > 0:
				TI_calc[turbIdx] = TI_a*(a[turbIdx]**TI_b)*(TI_0**TI_c)*(((xTurb[turbI]-xTurb[turbIdx])/D[turbIdx])**(TI_d))
			else:
				TI_calc[turbIdx] = 0.0
	
	n2 = time.time()
	
	# compute components of TI added
	TI_added = []
	TIidx = []
	Overlap = []
	for i in turbineID:
		if AreaOverlap[i]*TI_calc[i] > 0.0:
			TI_added.append(AreaOverlap[i]*TI_calc[i])
			TIidx.append(i) 
			Overlap.append(AreaOverlap[i])

	return TI_added, TIidx, Overlap

def computeOverlap(Ufield,UfieldOrig,inputData,Z):

	# compute wakeOverlap based on the number of points that are not freestream velocity, i.e. affected by a wake

	idx = Ufield.shape
	numPts = idx[0]*idx[1]*idx[2]
	count = 0.
	for i in range(idx[0]):
		for j in range(idx[1]):
			for k in range(idx[2]):
				if Ufield[i,j,k] < 0.99*UfieldOrig[i,j,k]:
					count = count + 1

	perOverlap = count/numPts
	return perOverlap

def computePower(Ueff,inputData):

	# compute the power at each turbine based on the effective velocity determined using the wake model

	nTurbs 		= len(inputData['turbineX'])
	rho 		= inputData['airDensity']
	D 			= inputData['rotorDiameter']
	pP 			= inputData['pP']
	pT 			= inputData['pT']
	gE 			= inputData['generatorEfficiency']
	CpCtTable	= inputData['TurbineInfo']['CpCtWindSpeed']

	# turbine operation
	yaw 		= inputData['yawAngles']
	tilt 		= inputData['tilt']
	Ct 			= inputData['Ct']
	Cp 			= inputData['Cp']

	powerOut 	= np.zeros(nTurbs)
	for i in range(nTurbs):
		A = np.pi*(D[i]/2.)**2
		Cptmp = Cp[i]*(np.cos(yaw[i]*np.pi/180.)**pP[i])*(np.cos((tilt[i])*np.pi/180.)**pT[i])
		powerOut[i] = 0.5*rho*A*Cptmp*gE[i]*Ueff[i]**3

	return powerOut

def initializeFlowField(X,Y,Z,inputData):

	# initialize the flow field used in the 3D model based on shear using the power log law

	# atmospheric conditions
	shear   = inputData['shear']
	HH 		= inputData['hubHeight'][0]
	Uinf 	= inputData['windSpeed']

	Ufield = np.zeros(X.shape)
	for i in range(X.shape[0]):
		for j in range(X.shape[1]):
			for k in range(X.shape[2]):
				Ufield[i,j,k] = Uinf*(Z[i,j,k]/HH)**shear

	return Ufield

def determineCpCt(GCPCT,windSpeed):

	# interpolation function based on Cp, Ct, and wind speed
    
	fCp = interp1d(GCPCT['wind_speed'],GCPCT['CP'])
	fCt = interp1d(GCPCT['wind_speed'],GCPCT['CT'])

	if windSpeed < np.min(GCPCT['wind_speed']):
		#print('here')
		Cp = np.max(GCPCT['CP'])
		Ct = 0.99
	else:    
		Cp = fCp(windSpeed)
		Ct = fCt(windSpeed)
    
	return Cp,Ct

def outputUpts(inputData,X,Y,Z,Ufield):

	# compute the output velocity at the points specified in inputData

	xPts = inputData['xPts']
	yPts = inputData['yPts']
	zPts = inputData['zPts']

	nSamplesX = X.shape[2]
	nSamplesY = Y.shape[1]
	nSamplesZ = Z.shape[0]

	Upts = np.zeros(len(xPts))
	count = 0
	#print(nSamplesX,nSamplesY,nSamplesZ)
	for kk in range(len(xPts)):
		nextPt = 0
	#while count < len(xPts):
		for i in range(nSamplesZ):
			for j in range(nSamplesY):
				for k in range(nSamplesX):
					if X[i,j,k] == xPts[kk] and Y[i,j,k] == yPts[kk] and Z[i,j,k] == zPts[kk]:
						if nextPt == 0:
							Upts[count] = Ufield[i,j,k]
							count = count + 1
						nextPt = 1

	return Upts







