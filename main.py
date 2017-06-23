## New analytical flow model - main file
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import time
from scipy.interpolate import griddata
from scipy import interpolate
from scipy.optimize import minimize
import utilities
import wakeModels
import OptModules
import copy

#==================================================================FUNCTIONS============================================================================
# 1. WakeModel - Takes care of all the high level computations such as generating grids to compute the wake model 
# 2. windPlant - main function to run the code (uses inputData)
# ======================================================================================================================================================

def WakeModel(inputData):

	# this function uses an engineering model to compute the time averaged wake of a wind turbine
	# you can model the wake using Gaussian, FLORIS, or Jensen parameters

	# wake parameters
	Uinf 			= inputData['windSpeed'] 				# wind speed
	windDirection 	= inputData['windDirection']			# wind direction
	TI_0			= inputData['turbulenceIntensity']		# turbulence intensity
	D 				= inputData['rotorDiameter']			# rotor diameter
	ke 				= inputData['wakeExpansion'] 			# wake expansion coefficient
	kd				= inputData['wakeDeflection']			# wake deflection coefficient
	ad 				= inputData['ad']						# lateral wake deflection a + b*X
	bd 				= inputData['bd']						# lateral wake deflection a + b*X
	rotorPts		= inputData['rotorPts']					# number of rotor points evaluated over the rotor
	rotorPts 		= int(np.round(np.sqrt(rotorPts)))		# evaluate the number of points in the y-z direction (equal in each direction)
	idxLidar 		= inputData['turbineLidar'] 			# turbine with the Lidar

	# Generate grid to evaluate the wake model
	nSamplesX 		= inputData['nSamplesX']				# number of X samples, used for visualization
	nSamplesY 		= inputData['nSamplesY']				# number of Y samples
	nSamplesZ 		= inputData['nSamplesZ']				# number of Z samples
	xLen 			= np.linspace(inputData['xLen'][0],inputData['xLen'][1],nSamplesX) 	# x domain
	yLen 			= np.linspace(inputData['yLen'][0],inputData['yLen'][1],nSamplesY) 	# y domain
	zLen 			= np.linspace(inputData['zLen'][0],inputData['zLen'][1],nSamplesZ)	# z domain
	dx 				= xLen[1]-xLen[0]						# grid spacing in the x direction
	dy 				= yLen[1]-yLen[0]						# grid spacing in the y direction
	nTurbs 			= len(inputData['turbineX'])			# number of turbines
	CpCtTable 		= inputData['TurbineInfo']['CpCtWindSpeed']	# Cp/Ct tables used to describe the turbine

	# turbine operation
	yaw 			= inputData['yawAngles']				# yaw angles for each turbine (deg)
	bladePitch 		= inputData['bladePitch']				# blade pitch angles for each turbine (deg)
	tilt 			= inputData['tilt']						# tilt angles for each turbine (deg)

	Ct 			  	= inputData['Ct']						# Cp for each turbine [-]
	Cp 				= inputData['Cp']						# Ct for each turbine [-]


	xTurb = inputData['turbineX']							# x locations of the turbines
	yTurb = inputData['turbineY']							# y locations
	zTurb = inputData['turbineZ'] 							# z lodcations


	# rotating the coordinates - flow always comes form the west
	xTurbOrig 		= xTurb
	yTurbOrig 		= yTurb
	zTurbOrig 		= zTurb

	# different grids generated based on input selection
	# want output of the flow field or to visualize the flow
	if inputData['outputFlowField'] or inputData['visualizeHorizontal']:
		if inputData['visualizeHorizontal']:
			inputData['nSamplesZ'] 	= 1
			zLen 					= inputData['hubHeight'][0]
			
		Z,Y,X 			   	= np.meshgrid(zLen,yLen,xLen,indexing='ij')
		X,Y,xTurb,yTurb 	= utilities.rotatedCoordinates(inputData,X,Y,Z)

	# Stuttgart lidar specifics
	elif inputData['Lidar']:

		print('Lidar operating on Turbine', inputData['turbineLidar'])

		xLidar = inputData['xLidar']
		yLidar = inputData['yLidar']
		zLidar = inputData['zLidar']

		# generate compressed grid
		X = []
		Y = []
		Z = []
		xTurb, yTurb 	= utilities.rotatedCoordinates(inputData,X,Y,Z)

		cols = yLidar.columns

		numPts = len(np.unique(yLidar[cols[0]]))
		numLocs = len(np.unique(xLidar))

		X = np.zeros((numPts,numPts,numLocs+nTurbs))
		Y = np.zeros((numPts,numPts,numLocs+nTurbs))
		Z = np.zeros((numPts,numPts,numLocs+nTurbs))

		for k in range(numLocs+nTurbs):

			if k < numLocs:
				Xl = xTurb[idxLidar] + np.unique(xLidar)[k]
				Yl = yTurb[idxLidar] + np.unique(yLidar[cols[k]])
				Zl = zTurb[idxLidar] + np.unique(zLidar[cols[k]])
			else:
				Xl = xTurb[k-numLocs]
				Yl = np.linspace(yTurb[k-numLocs]-(D[k-numLocs]/2), yTurb[k-numLocs]+(D[k-numLocs]/2), numPts)
				Zl = np.linspace(zTurb[k-numLocs]-(D[k-numLocs]/2), zTurb[k-numLocs]+(D[k-numLocs]/2), numPts)

			for j in range(len(Zl)):
				for i in range(len(Zl)):
					X[i,j,k] = Xl
					Y[i,j,k] = Yl[j]
					Z[i,j,k] = Zl[i]

	# output generic points of the flow field - Ex: used for custom lidar model
	elif inputData['points']:

		# generate compressed grid with individual x,y,z points
		X = []
		Y = []
		Z = []
		xTurb, yTurb 	= utilities.rotatedCoordinates(inputData,X,Y,Z)

		xPts = inputData['xPts']
		yPts = inputData['yPts']
		zPts = inputData['zPts']

		Xt = np.concatenate([xTurb, xPts])
		Yt = np.linspace(yTurb[0]-(D[0]/2), yTurb[0]+(D[0]/2), rotorPts)
		Zt = np.linspace(zTurb[0]-(D[0]/2), zTurb[0]+(D[0]/2), rotorPts)

		X = np.zeros((len(Zt),len(Yt),len(Xt)))
		Y = np.zeros((len(Zt),len(Yt),len(Xt)))
		Z = np.zeros((len(Zt),len(Yt),len(Xt)))

		# generate grid
		count = 0.
		for k in range(len(Xt)):
			if k < nTurbs:
				Yt = np.linspace(yTurb[k]-(D[k]/2), yTurb[k]+(D[k]/2), rotorPts)
				Zt = np.linspace(zTurb[k]-(D[k]/2), zTurb[k]+(D[k]/2), rotorPts)
			else:
				idx = k - nTurbs
				Yt = [yPts[idx]]
				Zt = [zPts[idx]]
			for j in range(len(Yt)):
				for i in range(len(Zt)):
					X[i,j,k] = Xt[k]
					Y[i,j,k] = Yt[j]
					Z[i,j,k] = Zt[i]


	# if you just want power output and effective velocities at each turbine:
	else:

		# generate compressed grid
		X = []
		Y = []
		Z = []
		xTurb, yTurb 	= utilities.rotatedCoordinates(inputData,X,Y,Z)

		Xt = xTurb 
		Yt = np.linspace(yTurb[0]-(D[0]/2), yTurb[0]+(D[0]/2), rotorPts)
		Zt = np.linspace(zTurb[0]-(D[0]/2), zTurb[0]+(D[0]/2), rotorPts)

		X = np.zeros((len(Zt),len(Yt),len(Xt)))
		Y = np.zeros((len(Zt),len(Yt),len(Xt)))
		Z = np.zeros((len(Zt),len(Yt),len(Xt)))

		count = 0.
		for k in range(len(Xt)):
			Yt = np.linspace(yTurb[k]-(D[k]/2), yTurb[k]+(D[k]/2), rotorPts)
			Zt = np.linspace(zTurb[k]-(D[k]/2), zTurb[k]+(D[k]/2), rotorPts)
			for j in range(len(Yt)):
				for i in range(len(Zt)):
					X[i,j,k] = Xt[k]
					Y[i,j,k] = Yt[j]
					Z[i,j,k] = Zt[i]


	# initialize z location
	zTurb 				= np.array(zTurb)

	# determine the amount of the grid to loop through (depends on if outputting the velocity field)
	nSamplesZ 			= Z.shape[0]	
	nSamplesY 			= Y.shape[1]
	nSamplesX 			= X.shape[2]

	# sort turbine coordinates from front to back
	turbineID		= [i[0] for i in sorted(enumerate(xTurb), key=lambda x:x[1])]

	# initialize flow field with a uniform shear layer
	Ufield 			= utilities.initializeFlowField(X,Y,Z,inputData)
	UfieldOrig		= copy.copy(Ufield)
	Ueff	 		= np.zeros(nTurbs)

	TIfield 		= np.zeros((Z.shape[0],Y.shape[1],X.shape[2]))	

	# make sure a wake model is specified
	initWake		= 0

	# turbine operation
	a 				= np.zeros(nTurbs)
	if inputData['TurbineInfo']['PitchCpCt'] == True:
		fCp,fCt,beta 	= inputData['TurbineInfo']['CpCtPitch']


	for turbI in turbineID:
		
		# compute effective wind speed at each turbine
		# take the average across the rotor disk
		Ueff[turbI] = utilities.avgVelocity(X,Y,Z,Ufield,xTurb[turbI],yTurb[turbI],zTurb[turbI],D[turbI],rotorPts,inputData,turbI)

		# adjust Cp/Ct based on local velocity (both tables were generated using FAST, Jonkman et. al. 2005)

		# this table is based on pitch and wind speed
		if inputData['TurbineInfo']['PitchCpCt'] == True:
			Cp[turbI] = fCp(Ueff[turbI],bladePitch[turbI])
			Ct[turbI] = fCt(Ueff[turbI],bladePitch[turbI])

		# this table is based on wind speed only 
		else:
			Cp[turbI],Ct[turbI] = utilities.determineCpCt(CpCtTable,Ueff[turbI])

		if Ct[turbI] >= 1.0:
			Ct[turbI] = 0.99999999
		inputData['Ct'][turbI] = Ct[turbI]
		inputData['Cp'][turbI] = Cp[turbI]
		a[turbI] = (0.5/np.cos(yaw[turbI]*np.pi/180.))*(1-np.sqrt(1-Ct[turbI]*np.cos(yaw[turbI]*np.pi/180)))

		# compute the x, y, z offset due to yaw
		zR = (D[turbI]/2.)*np.cos(np.radians(tilt[turbI]))
		yR = (D[turbI]/2.)*np.cos(np.radians(yaw[turbI]))
		xR = (D[turbI]/2.)*np.sin(np.radians(yaw[turbI])) + (D[turbI]/2.)*np.sin(np.radians(tilt[turbI]))

		# compute the initial added turbulence at the rotor
		TI_added, TIidx, AreaOverlap	= utilities.computeAddedTI(X,Y,Z,inputData,xTurb,yTurb,zTurb,turbI,turbineID)

		# add turbulence via sum of squares
		tmpTI = TI_0**2
		for i in range(len(TI_added)):
			tmpTI = tmpTI + TI_added[i]**2
		TI 	= np.sqrt(tmpTI)				# total TI added to the downstream turbine

		# cycle through the grid generated above
		for xIdx in range(nSamplesX):
			for yIdx in range(nSamplesY):
				for zIdx in range(nSamplesZ):

					if X[zIdx,yIdx,xIdx] == 0.0 and Y[zIdx,yIdx,xIdx] == 0.0 and Z[zIdx,yIdx,xIdx] == 0.0:
						continue

					else:
						# use Gaussian wake model
						if inputData['WakeModel'] == 2: 

							if X[zIdx,yIdx,xIdx] >= xTurb[turbI] - D[turbI]:
								c 		 	= wakeModels.GAUSS(X[zIdx,yIdx,xIdx],Y[zIdx,yIdx,xIdx],Z[zIdx,yIdx,xIdx],
															   xTurb[turbI],yTurb[turbI],zTurb[turbI],
															   inputData,turbI,TI,UfieldOrig[zIdx,yIdx,xIdx])
								velDef 	 	= c
							
							else:
								velDef = 0.0

							# save value for wake merging 
							if Ufield[zIdx,yIdx,xIdx] != UfieldOrig[zIdx,yIdx,xIdx]:
								tmp = UfieldOrig[zIdx,yIdx,xIdx] - Ufield[zIdx,yIdx,xIdx]
							else:
								tmp = 0.0


						# use FLORIS or Jensen wake model
						elif inputData['WakeModel'] == 0 or inputData['WakeModel'] == 1:	
							if (X[zIdx,yIdx,xIdx] > (xTurb[turbI]-abs(xR)) and 
							   Y[zIdx,yIdx,xIdx] > (yTurb[turbI] - 2*D[turbI]) and 
							   Y[zIdx,yIdx,xIdx] < (yTurb[turbI] + 2*D[turbI]) and
							   Z[zIdx,yIdx,xIdx] > (zTurb[turbI] - 2*D[turbI]) and
							   Z[zIdx,yIdx,xIdx] < (zTurb[turbI] + 2*D[turbI])):
								yDisp = wakeModels.Jimenez(np.radians(yaw[turbI]),Ct[turbI],kd[turbI],X[zIdx,yIdx,xIdx]-xTurb[turbI],D[turbI],ad,bd)
								zDisp = wakeModels.Jimenez(np.radians(tilt[turbI]),Ct[turbI],kd[turbI],X[zIdx,yIdx,xIdx]-xTurb[turbI],D[turbI],ad,bd)
							else:
								yDisp = 0.0
								zDisp = 0.0
								continue

							# define the edges of the wake
							rWake = ke[turbI]*(X[zIdx,yIdx,xIdx]-(xTurb[turbI]-xR)) + (D[turbI]/2)
							rCenterY = yTurb[turbI] + yDisp
							rCenterZ = zTurb[turbI] + zDisp

							rtmp = np.sqrt( (Y[zIdx,yIdx,xIdx] - rCenterY)**2 + (Z[zIdx,yIdx,xIdx] - rCenterZ)**2 )
							if (X[zIdx,yIdx,xIdx] >= xTurb[turbI] and rtmp <= rWake):
								
								# FLORIS model
								if inputData['WakeModel']==1:
									c = wakeModels.FLORIS(X[zIdx,yIdx,xIdx],Y[zIdx,yIdx,xIdx],Z[zIdx,yIdx,xIdx],yDisp,zDisp,xTurb[turbI],yTurb[turbI],zTurb[turbI],inputData,turbI)
									velDef = Uinf*2.*a[turbI]*c
								
								# Jensen model
								elif inputData['WakeModel']==0:
									c = wakeModels.Jensen(inputData,turbI,X[zIdx,yIdx,xIdx],xTurb[turbI])
									velDef = Uinf*2.*a[turbI]*c

							else:
								velDef = 0.0

							# save value for wake merging 
							if Ufield[zIdx,yIdx,xIdx] != UfieldOrig[zIdx,yIdx,xIdx]:
								tmp = UfieldOrig[zIdx,yIdx,xIdx] - Ufield[zIdx,yIdx,xIdx]
							else:
								tmp = 0.0

						else:
							if initWake == 0:
								print('No wake model specified')
								velDef = 0.0
								tmp = 0.0
								initWake = 1
								break

						# update the velocity field
						Ufield[zIdx,yIdx,xIdx] = UfieldOrig[zIdx,yIdx,xIdx] - velDef

						if tmp != 0.0:
							Ufield[zIdx,yIdx,xIdx] = utilities.combineWakes(inputData,UfieldOrig[zIdx,yIdx,xIdx],Ueff[turbI],Ufield[zIdx,yIdx,xIdx],tmp)

	# plot the flow field if specified
	if inputData['visualizeHorizontal']:
		utilities.visualizeHorizontal(xLen,yLen,zLen,Ufield,inputData)

	elif inputData['Lidar']:
		utilities.visualizeLidar(xTurb[idxLidar],yTurb[idxLidar],X,Y,Z,Ufield,inputData)
		
		
	if inputData['outputFlowField']:
		return Ueff,Ufield,xLen,yLen
	elif inputData['Lidar']:
		return Ueff,Ufield,X,Y,Z
	elif inputData['points']:
		Upts = utilities.outputUpts(inputData,X,Y,Z,Ufield)
		return Ueff, Upts
	else:
		return Ueff


def windPlant(inputData):

	# time the execution of this model
	startTime = time.time()
	
	# initialize parameters------------------------------
	windSpeed 		= inputData['windSpeed']
	windDirection 	= inputData['windDirection']
	nTurbs 			= len(inputData['turbineX'])
	yaw 			= inputData['yawAngles']
	tilt 			= inputData['tilt']
	Cp 				= inputData['Cp']
	Ct 				= inputData['Ct']

	powerOut 		= np.zeros(nTurbs)
	outputData 		= dict()
	# ---------------------------------------------------

	# reversed sign convention for Gaussian wake model
	if inputData['WakeModel'] == 2:
		inputData['yawAngles'] 	= np.array(yaw)
		inputData['tilt'] 		= np.array(tilt)

	# ============================================================
	# Outputs based on input selection
	# ============================================================

	# output flow field
	if inputData['outputFlowField']:
		Ueff,Ufield,xLen,yLen 	= WakeModel(inputData)
		outputData['Ufield'] 	= Ufield
		outputData['xLen']	= xLen
		outputData['yLen'] 	= yLen
	# output lidar parameters
	elif inputData['Lidar']:
		Ueff,Ufield,X,Y,Z 		= WakeModel(inputData)
		outputData['Ufield'] 	= Ufield
		outputData['X']			= X
		outputData['Y']			= Y
		outputData['Z']			= Z
	# output velocities at specified points
	elif inputData['points']:
		Ueff,Upts 				= WakeModel(inputData)
		outputData['Upts']		= Upts
	# output effective velocities at each turbine
	else:
		Ueff 					= WakeModel(inputData)		

	powerOut 				= utilities.computePower(Ueff,inputData)
	outputData['powerOut'] 	= powerOut
	outputData['Ueff']      = Ueff

	# initial power output for optimization comparison
	power0 = np.sum(outputData['powerOut'])

	# =============================================================
	# perform optimization(s)
	# =============================================================

	if inputData['axial_opt']:
		CpOpt,CtOpt,bladePitchOpt = OptModules.axialOpt(inputData)
		inputData['Cp'] 		= CpOpt
		inputData['Ct'] 		= CtOpt		
		outputData['Cp_opt'] 	= CpOpt
		outputData['Ct_opt'] 	= CtOpt
		inputData['bladePitch'] = bladePitchOpt
	elif inputData['yaw_opt']:
		fileloc 				= []
		yawOpt 					= OptModules.yawOpt(inputData)
		inputData['yawAngles'] 	= yawOpt
		outputData['yaw_opt'] 	= yawOpt

	# ==============================================================
	# rerun the wake model with the optimized parameters
	# ==============================================================

	if inputData['axial_opt'] or inputData['yaw_opt']:
		Ueff 					= WakeModel(inputData)			

		powerOut 				= utilities.computePower(Ueff,inputData)
		outputData['powerOut'] 	= powerOut
		outputData['Ueff']		= Ueff
		outputData['yawAngles'] = inputData['yawAngles']
		outputData['bladePitch'] = inputData['bladePitch']

		powerOpt = np.sum(outputData['powerOut'])
		powerGain = 100*(powerOpt - power0)/power0
		print('Power gain = ', powerGain)
		if powerGain < 0.0:
			outputData['bladePitch'] = 1.9*np.ones(len(inputData['turbineX']))
			outputData['yawAngles'] = np.zeros(len(inputData['turbineX']))

	return outputData 

	

