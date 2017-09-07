## Utilities for wake model

import numpy as np
import time
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import scipy.io
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import wakeModels
import copy
import time
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

def visualizeCut(xLen,yLen,zLen,Ufield,inputData):


	print('Plotting Flow Field Cut Through Slice at ', inputData['downLocs'], '...\n')

	Uinf 			= inputData['windSpeed'][0]
	windDirection 	= inputData['windDirection'][0]
	tilt 			= -1*inputData['tilt']
	yaw 			= inputData['yawAngles']
	D 				= inputData['rotorDiameter']
	HH 				= inputData['hubHeight']

	xTurb 			= inputData['turbineX']
	yTurb 			= inputData['turbineY']
	zTurb 			= inputData['turbineZ']

	# number of turbines
	nTurbs = len(xTurb)

	# rotation angle for the turbines the wind turbines
	RotAng = -(windDirection -270.)*np.pi/180.

	# plot horizontal flow field
	plt.figure(figsize=(10,7))
	plt.contourf(-yLen,zLen,Ufield[:,:,0],50,cmap='coolwarm',vmin=4.0,vmax=Uinf)
	plt.colorbar()

	# plot rotor 
	yCirc = np.linspace(-D[0]/2,D[0]/2,100)
	zCirc1 = np.sqrt( (D[0]/2)**2 - (yCirc-yTurb[0])**2 ) + zTurb[0]
	zCirc2 = -np.sqrt( (D[0]/2)**2 - (yCirc-yTurb[0])**2 ) + zTurb[0]
	plt.plot(yCirc,zCirc1,'k',linewidth=3)
	plt.plot(yCirc,zCirc2,'k',linewidth=3)

	#for i in range(len(xLen)):
	#	for j in range(len(yLen)):
	#		plt.plot(xLen[i],yLen[j],'k.',linewidth=3)

	#plt.axis('equal')
	plt.xlabel('y(m)',fontsize=15)
	plt.ylabel('z(m)',fontsize=15)
	strTitle = 'Cut Through Slice at ', inputData['downLocs'][0]/D[0],'D'
	plt.title(strTitle,fontsize=15)
	plt.xlim([-2*D[0],2*D[0]])
	plt.ylim([0.0,2*HH[0]])
	plt.axis('equal')
	

def visualizeLidar(xTurb,yTurb,X,Y,Z,Ufield,inputData,vlos):

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
				#print(xLocs[kk],yPlot[jj],zPlot[ii])
				if X[ii,jj,kk] == xLocs[kk] and Y[ii,jj,kk] == yPlot[jj] and Z[ii,jj,kk] == zPlot[ii]:
					Uplot[ii,jj] = Ufield[ii,jj,kk]

		plt.contourf(-yPlot,zPlot,Uplot,50,cmap='coolwarm',vmin=2.0,vmax=Uinf)
		plt.xlim([-60-yTurb,60-yTurb])
		plt.ylim([40,160])
		plt.colorbar()


def getWindCoords(inputData):

	Parameter = scipy.io.loadmat(inputData['LidarParams'])
	Trajectory = scipy.io.loadmat(inputData['LidarTrajectory'])	

	a = Parameter['Parameter'][0][0][0]['a'][0][0][0]
	Yaw = Parameter['Parameter'][0][0][0]['YawAngle'][0][0][0]
	Pitch = Parameter['Parameter'][0][0][0]['PitchAngle'][0][0][0]
	Roll = Parameter['Parameter'][0][0][0]['RollAngle'][0][0][0]
	PositionLinW = Parameter['Parameter'][0][0][0]['PositionLinI'][0][0][0]
	
	x_L = np.concatenate(Trajectory['Trajectory'][0][0]['x_L_AllDistances'])
	y_L = np.concatenate(Trajectory['Trajectory'][0][0]['y_L_AllDistances'])
	z_L = np.concatenate(Trajectory['Trajectory'][0][0]['z_L_AllDistances'])

	# yaw is a rotation around the z-axis
	T_Yaw = [[np.cos(Yaw),-np.sin(Yaw), 0],
	         [np.sin(Yaw), np.cos(Yaw), 0],
	         [0,0,1]]

	T_Pitch = [[np.cos(Yaw), 0, -np.sin(Yaw)],
	           [0,1,0],
	           [np.sin(Yaw), 0, np.cos(Yaw)]]

	T_Roll = [[1,0,0],
	          [0,np.cos(Roll),-np.sin(Roll)],
	          [0,np.sin(Roll), np.cos(Roll)]]

	T = np.dot(np.dot(T_Yaw,T_Pitch),T_Roll)

	x_R = T[0,0]*x_L + T[0,1]*y_L + T[0,2]*z_L
	y_R = T[1,0]*x_L + T[1,1]*y_L + T[1,2]*z_L
	z_R = T[2,0]*x_L + T[2,1]*y_L + T[2,2]*z_L

	#print(PositionLinW)

	x_W = x_R + PositionLinW[0] + inputData['turbineX'][inputData['turbineLidar']]
	y_W = y_R + PositionLinW[1] + inputData['turbineY'][inputData['turbineLidar']]
	z_W = z_R + PositionLinW[2] + inputData['hubHeight'][inputData['turbineLidar']]

	for i in range(len(x_W)):
		if z_W[i] < 0:
			z_W[i] = 0.01

	#print(-y_W,y_W)

	return x_W, y_W, z_W

def getPointsVLOS(x_W,y_W,z_W,inputData):

	# retrieves the line of sight wind speed (vlos) from a wind field

	# input: parameter, wind field
	# output: v_los

	#print(len(x_W))

	Parameter = scipy.io.loadmat(inputData['LidarParams'])
	Trajectory = scipy.io.loadmat(inputData['LidarTrajectory'])

	a 			= Parameter['Parameter'][0][0][0]['a'][0][0][0]
	nWeights 	= len(a)	
	nDataPoints = len(x_W)

	#print(x_W)

	#x_LW = np.zeros(nDataPoints) 
	#y_LW = np.zeros(nDataPoints)
	#z_LW = np.zeros(nDataPoints)
	x_LW = np.ones(nDataPoints)*inputData['turbineX'][inputData['turbineLidar']]
	y_LW = np.ones(nDataPoints)*inputData['turbineY'][inputData['turbineLidar']]
	z_LW = np.ones(nDataPoints)*inputData['hubHeight'][inputData['turbineLidar']]

	# calculation of the normalized laser vector
	LaserVector_Wx = [np.transpose(x_W)-np.transpose(x_LW)][0]
	LaserVector_Wy = [np.transpose(y_W)-np.transpose(y_LW)][0]
	LaserVector_Wz = [np.transpose(z_W)-np.transpose(z_LW)][0]

	#print(LaserVector_Wx)

	#print(np.transpose(x_W),np.transpose(x_LW))

	#print(LaserVector_Wx)

	# np.sqrt(x**2 + y**2 + z**2) -> xhat = x/Norm
	NormedLaserVector_Wx = np.zeros(len(LaserVector_Wx))
	NormedLaserVector_Wy = np.zeros(len(LaserVector_Wy))
	NormedLaserVector_Wz = np.zeros(len(LaserVector_Wz))

	#print(len(x_LW))
	for i in range(nDataPoints):
		NormLaserVector_W = np.sqrt( LaserVector_Wx[i]**2 + LaserVector_Wy[i]**2 + LaserVector_Wz[i]**2 )
		#print(np.sqrt( LaserVector_Wx[i]**2 + LaserVector_Wy[i]**2 + LaserVector_Wz[i]**2 ))
		#print(LaserVector_Wx[i],LaserVector_Wy[i],LaserVector_Wz[i],NormLaserVector_W)
		if NormLaserVector_W == 0:
			NormLaserVector_W = 1
		NormedLaserVector_Wx[i] = LaserVector_Wx[i]/NormLaserVector_W
		NormedLaserVector_Wy[i] = LaserVector_Wy[i]/NormLaserVector_W
		NormedLaserVector_Wz[i] = LaserVector_Wz[i]/NormLaserVector_W

	BackscatterNormedLaserVector_Wx = -NormedLaserVector_Wx
	BackscatterNormedLaserVector_Wy = -NormedLaserVector_Wy
	BackscatterNormedLaserVector_Wz = -NormedLaserVector_Wz

	# Calculation of considered points in the laser beam
	Points_WFx = np.zeros(nDataPoints*nWeights)
	Points_WFy = np.zeros(nDataPoints*nWeights)
	Points_WFz = np.zeros(nDataPoints*nWeights)
	for i in range(nDataPoints):
		#print(x_W[i]*np.ones(nWeights))
		Points_WFx[i*nWeights:((i*nWeights)+nWeights)] = x_W[i]*np.ones(nWeights) + (a*BackscatterNormedLaserVector_Wx[i]*np.ones(nWeights))
		Points_WFy[i*nWeights:((i*nWeights)+nWeights)] = y_W[i]*np.ones(nWeights) + (a*BackscatterNormedLaserVector_Wy[i]*np.ones(nWeights))
		Points_WFz[i*nWeights:((i*nWeights)+nWeights)] = z_W[i]*np.ones(nWeights) + (a*BackscatterNormedLaserVector_Wz[i]*np.ones(nWeights))

	#print(Points_WFx)
	#print(inputData['turbineX'][inputData['turbineLidar']])

	Points_WFy = Points_WFy - inputData['turbineY'][inputData['turbineLidar']]
	Points_WFy = inputData['turbineY'][inputData['turbineLidar']] - Points_WFy

	return Points_WFx,Points_WFy,Points_WFz

def VLOS(x_W,y_W,z_W,inputData,Upts):

	#print(x_W)

	Parameter = scipy.io.loadmat(inputData['LidarParams'])
	Trajectory = scipy.io.loadmat(inputData['LidarTrajectory'])

	a 			= Parameter['Parameter'][0][0][0]['a'][0][0][0]
	nWeights 	= len(a)	
	nDataPoints = len(x_W)
	f_L_d 		= Parameter['Parameter'][0][0][0]['f_L_d'][0][0][0]

	#nWeights = len(Parameter.Lidar.a)
	x_L = np.concatenate(Trajectory['Trajectory'][0][0]['x_L_AllDistances'])
	y_L = np.concatenate(Trajectory['Trajectory'][0][0]['y_L_AllDistances'])
	z_L = np.concatenate(Trajectory['Trajectory'][0][0]['z_L_AllDistances']) 
	#t    = Parameter.t

	# origin of the lidar and the origin of the wind coordinate system 
	#x_LW = np.zeros(nDataPoints)
	#y_LW = np.zeros(nDataPoints)
	#z_LW = np.zeros(nDataPoints)
	x_LW = np.ones(nDataPoints)*inputData['turbineX'][inputData['turbineLidar']]
	y_LW = np.ones(nDataPoints)*inputData['turbineY'][inputData['turbineLidar']]
	z_LW = np.ones(nDataPoints)*inputData['hubHeight'][inputData['turbineLidar']]

	# calculation of the normalized laser vector (vector from the zero in the wind to the trajectory point in the wind)
	# lidar in the wind does not need to be at zero
	LaserVector_Wx = [np.transpose(x_W)-np.transpose(x_LW)][0]
	LaserVector_Wy = [np.transpose(y_W)-np.transpose(y_LW)][0]
	LaserVector_Wz = [np.transpose(z_W)-np.transpose(z_LW)][0]

	NormedLaserVector_Wx = np.zeros(len(LaserVector_Wx))
	NormedLaserVector_Wy = np.zeros(len(LaserVector_Wy))
	NormedLaserVector_Wz = np.zeros(len(LaserVector_Wz))

	for i in range(nDataPoints):
		NormLaserVector_W = np.sqrt( LaserVector_Wx[i]**2 + LaserVector_Wy[i]**2 + LaserVector_Wz[i]**2 )
		#print(LaserVector_Wx[i],LaserVector_Wy[i],LaserVector_Wz[i],NormLaserVector_W)
		NormedLaserVector_Wx[i] = LaserVector_Wx[i]/NormLaserVector_W
		NormedLaserVector_Wy[i] = LaserVector_Wy[i]/NormLaserVector_W
		NormedLaserVector_Wz[i] = LaserVector_Wz[i]/NormLaserVector_W

	BackscatterNormedLaserVector_Wx = NormedLaserVector_Wx
	BackscatterNormedLaserVector_Wy = NormedLaserVector_Wy
	BackscatterNormedLaserVector_Wz = NormedLaserVector_Wz

	#inputData = FLORISInput(Points_WFx,Points_WFy,Points_WFz)
	u_W = Upts
	v_W = np.zeros(nDataPoints*nWeights)
	w_W = np.zeros(nDataPoints*nWeights) 

	RelativeWindVector_W = [u_W,v_W,w_W]
	BackScatterNormed_W  = [BackscatterNormedLaserVector_Wx,BackscatterNormedLaserVector_Wy,BackscatterNormedLaserVector_Wz]

	for i in range(3):
		tmp = BackScatterNormed_W[i]
		for j in range(nDataPoints):
			if j == 0:
				tmp1 = tmp[j]*np.ones(nWeights)
			else:
				tmp1 = np.concatenate((tmp1,tmp[j]*np.ones(nWeights)))	
		if i == 0:
			BackscatterNormedLaserVector_Wx = tmp1	
		elif i == 1:
			BackscatterNormedLaserVector_Wy = tmp1
		else:
			BackscatterNormedLaserVector_Wz = tmp1

	#print(BackscatterNormedLaserVector_Wy)

	BackScatterNormed_W1 = [BackscatterNormedLaserVector_Wx,BackscatterNormedLaserVector_Wy,BackscatterNormedLaserVector_Wz]

	vlos_d = np.multiply(RelativeWindVector_W,BackScatterNormed_W1)
	v_los = np.zeros(nDataPoints)

	for i in range(nDataPoints):
		v_los[i] = np.dot(vlos_d[0,i*nWeights:((i*nWeights)+nWeights)],f_L_d)

	return v_los

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

	if inputData['avgCube']:
		return np.mean(utmp**3)**(1./3.)
	else:
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
				#print(Z[i,j,k],Ufield[i,j,k])

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

	if inputData['points']:
		xPts = inputData['xPts']
		yPts = inputData['yPts']
		zPts = inputData['zPts']
	elif inputData['Lidar']:
		x_W, y_W, z_W = getWindCoords(inputData)
		xTraj,yTraj,zTraj = getPointsVLOS(x_W,y_W,z_W,inputData)
		xPts = xTraj
		yPts = yTraj
		zPts = zTraj

	nSamplesX = X.shape[2]
	nSamplesY = Y.shape[1]
	nSamplesZ = Z.shape[0]

	#print(np.shape(xPts))
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

	#print(np.shape(Upts))

	return Upts







