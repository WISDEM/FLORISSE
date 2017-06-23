## Different Wake models

import numpy as np
import utilities

#==================================================================FUNCTIONS============================================================================
# 1. Jimenez - empirical calculation used to determine the wake deflection due to yaw 
# 2. Jensen - original linear wind plant model that uses a constant velocity in the wake. 
# 3. FLORIS - linear wake model that uses three wake zones (near, far, and mixing) to define the wake.
# 4. GAUSS - Gaussian wake based on self-similarity theory 
# ======================================================================================================================================================

# =================================================================WAKE MODELS=======================================================

def Jimenez(yaw,Ct,kd,x,D,ad,bd):

	# this function defines the angle at which the wake deflects in relation to the yaw of the turbine
	# this is coded as defined in the Jimenez et. al. paper

	# angle of deflection
	xi_init = 1./2.*np.cos(yaw)*np.sin(yaw)*Ct
	xi = (xi_init)/(( 1 + 2*kd*(x/D) )**2)

	# yaw displacement
	yYaw_init = ( xi_init * ( 15*((2*kd*x/D) + 1)**4. + xi_init**2. ) / ( (30*kd/D)*(2*kd*x/D + 1)**5. ) ) - ( xi_init*D*(15 + xi_init**2.) / (30*kd) )

	# corrected yaw displacement with lateral offset
	yYaw = yYaw_init + ( ad + bd*x )

	return yYaw

def Jensen(inputData,turbI,X,xTurb):

	# compute the velocity deficit based on the classic Jensen/Park model, see Jensen 1983
	D 	= inputData['rotorDiameter'][turbI]
	ke 	= inputData['wakeExpansion'][turbI]
	c 	= (D/(2*ke*(X-xTurb) + D))**2

	return c

def FLORIS(X,Y,Z,yDisp,zDisp,xTurb,yTurb,zTurb,inputData,turbI):

	# compute the velocity deficit based on wake zones, see Gebraad et. al. 2016

	# wake parameters
	ke 		= inputData['wakeExpansion'][turbI]
	MU 		= inputData['MU']
	me 		= inputData['me']
	aU 		= inputData['aU']
	bU 		= inputData['bU']
	D  		= inputData['rotorDiameter'][turbI]
	yaw 	= inputData['yawAngles'][turbI]

	# distance from wake centerline
	rY 		= abs( Y-(yTurb+yDisp) )
	rZ		= abs( Z-(zTurb+zDisp) ) 		
	dx 		= X-xTurb

	# wake zone diameters
	d1 		= ( D + 2*ke*me[0]*dx )/2.
	d2 		= ( D + 2*ke*me[1]*dx )/2.
	d3 		= ( D + 2*ke*me[2]*dx )/2.

	# defining wake parameters
	mU_0 	= MU[0] / (np.cos( np.radians(aU + bU*yaw) ))
	mU_1 	= MU[1] / (np.cos( np.radians(aU + bU*yaw) ))
	mU_2	= MU[2] / (np.cos( np.radians(aU + bU*yaw) ))

	# near wake zone
	if rY <= d1:
		c 	= ( D / (D + (2*ke*mU_0*dx)) )**2
	# far wake zone
	elif rY <= d2:
		c 	= ( D / (D + (2*ke*mU_1*dx)) )**2
	# mixing zone
	else:
		c 	= ( D / (D + (2*ke*mU_2*dx)) )**2

	return c

def GAUSS(X,Y,Z,xTurb,yTurb,zTurb,inputData,turbI,TI,Uloc):

	# new analytical wake model based on self-similarity and Gaussian wake model
	# Gaussian Analytical U* Self-Similar Wake model

	# input data
	Uinf 		= inputData['windSpeed']
	D 			= inputData['rotorDiameter'][turbI]
	TI_0		= inputData['turbulenceIntensity']
	ka			= inputData['ka']
	kb 			= inputData['kb']
	alpha		= inputData['alpha']
	beta 		= inputData['beta']	
	HH 			= inputData['hubHeight'][turbI]
	yaw 	    = inputData['yawAngles'][turbI] 						# sign-convention issue
	tilt		= inputData['tilt'][turbI]
	Ct 			= inputData['Ct'][turbI]
	veer 		= inputData['veer']

	ad 			= inputData['ad']
	bd 			= inputData['bd']
	aT 			= inputData['aT']
	bT 			= inputData['bT']

	# sign convention reversed in literature
	yaw 		= -1*yaw
	
	# initial velocity deficits
	uR 			= Uloc*Ct*np.cos(yaw*np.pi/180.)/(2.*(1-np.sqrt(1-(Ct*np.cos(yaw*np.pi/180.)))))
	u0 			= Uloc*np.sqrt(1-Ct)

	# initial Gaussian wake expansion
	sigma_z0 	= D*0.5*np.sqrt( uR/(Uloc + u0) )
	sigma_y0 	= sigma_z0*(np.cos((yaw)*np.pi/180.))*(np.cos(veer*np.pi/180.))

    # quantity that determines when the far wake starts
	x0 			= D*(np.cos(yaw*np.pi/180.)*(1+np.sqrt(1-Ct*np.cos(yaw*np.pi/180.)))) / (np.sqrt(2)*(4*alpha*TI + 2*beta*(1-np.sqrt(1-Ct)))) + xTurb

    # wake expansion parameters
	ky 			= ka*TI + kb 
	kz			= ka*TI + kb

	C0 			= 1 - u0/Uinf
	M0 			= C0*(2-C0)	
	E0 			= C0**2 - 3*np.exp(1./12.)*C0 + 3*np.exp(1./3.)

	# yaw parameters (skew angle and distance from centerline)
	theta_c0 	= 2*((0.3*yaw*np.pi/180.)/np.cos(yaw*np.pi/180.))*(1-np.sqrt(1-Ct*np.cos(yaw*np.pi/180.)))	# skew angle  
	theta_z0 	= 2*((0.3*tilt*np.pi/180.)/np.cos(tilt*np.pi/180.))*(1-np.sqrt(1-Ct*np.cos(tilt*np.pi/180.)))	# skew angle  
	delta0 		= np.tan(theta_c0)*(x0-xTurb)
	delta_z0 	= np.tan(theta_z0)*(x0-xTurb)

	## COMPUTE VELOCITY DEFICIT
	yR = Y - yTurb
	xR = yR*np.tan(yaw*np.pi/180.) + xTurb

	if X > xR:

		if X >= (x0):
			sigma_y = ky*( X - x0 ) + sigma_y0
			sigma_z = kz*( X - x0 ) + sigma_z0
			ln_deltaNum = (1.6+np.sqrt(M0))*(1.6*np.sqrt(sigma_y*sigma_z/(sigma_y0*sigma_z0)) - np.sqrt(M0))
			ln_deltaDen = (1.6-np.sqrt(M0))*(1.6*np.sqrt(sigma_y*sigma_z/(sigma_y0*sigma_z0)) + np.sqrt(M0))
			delta = delta0 + (theta_c0*E0/5.2)*np.sqrt(sigma_y0*sigma_z0/(ky*kz*M0))*np.log(ln_deltaNum/ln_deltaDen) + ( ad + bd*(X-xTurb) )  
			deltaZ = delta_z0 + (theta_z0*E0/5.2)*np.sqrt(sigma_y0*sigma_z0/(ky*kz*M0))*np.log(ln_deltaNum/ln_deltaDen) + ( aT + bT*(X-xTurb) )
		else: 
			sigma_y = (((x0-xR)-(X-xR))/(x0-xR))*0.501*D*np.sqrt(Ct/2.) + ((X-xR)/(x0-xR))*sigma_y0
			sigma_z = (((x0-xR)-(X-xR))/(x0-xR))*0.501*D*np.sqrt(Ct/2.) + ((X-xR)/(x0-xR))*sigma_z0
			delta = ((X-xR)/(x0-xR))*delta0 + ( ad + bd*(X-xTurb) )	
			deltaZ = ((X-xR)/(x0-xR))*delta_z0 + ( aT + bT*(X-xTurb) )

		a = (np.cos(veer*np.pi/180.)**2)/(2*sigma_y**2) + (np.sin(veer*np.pi/180.)**2)/(2*sigma_z**2)
		b = -(np.sin(2*veer*np.pi/180))/(4*sigma_y**2) + (np.sin(2*veer*np.pi/180.))/(4*sigma_z**2)
		c = (np.sin(veer*np.pi/180.)**2)/(2*sigma_y**2) + (np.cos(veer*np.pi/180.)**2)/(2*sigma_z**2)
		totGauss = np.exp( -( a*((Y-yTurb)-delta)**2 - 2*b*((Y-yTurb)-delta)*((Z-zTurb)-deltaZ) + c*((Z-zTurb)-deltaZ)**2 ) ) 

		velDef = Uloc*(1-np.sqrt(1-((Ct*np.cos(yaw*np.pi/180.))/(8.0*sigma_y*sigma_z/D**2)) ) )*totGauss

	else:
		velDef = 0.0

	c = velDef

	return c





