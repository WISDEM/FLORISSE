// Pieter Gebraad, 2014
// FLOw Redirection and Induction Steady-State model
// given effective wind direction, predicts the powers of turbines in a wind farm as a function of yaw and Cp under non-yawed conditions

// Modified by S. Andrew Ning
// April 22, 2014
// changed float -> double
// changed input structure and removed structs so more easily callable from Python
// use M_PI already defined in math library
// switched pointer/malloc/free constructs to arrays with fixed size

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Because some compilers are missing these functions, possibly define here
#ifndef fmin
	double fmin(double x, double y) { return(x < y ? x : y); }
	double fmax(double x, double y) { return(x > y ? x : y); }
#endif

// void FLORISmodel(float * power, int numberOfTurbines, turbine * TURB, inflow INF, modelParams mP, int flagUseMeasPower)
void FLORISmodel(double power[], double effU[], double effU_In[], double X[], double Y[], double axialInd[],
	double yaw[], double rotorDiameter[], double Cp[], double measPower[], int nTurb,
	double effUdXY, double rho,
	double kd, double ad, double bd, double ke, double me1, double me2, double me3,
	double pP, double MU1, double MU2, double MU3, double aU, double bU,
	int flagUseMeasPower){

	#ifndef M_PI
		const double M_PI = 3.1415926535897932385;
	#endif

	// other variables
	double * rotorA  =  malloc(nTurb * sizeof(double));
	double * turbXd  =  malloc(nTurb * sizeof(double));			// downwind-crosswind coordinates turbine
	double * turbYd  =  malloc(nTurb * sizeof(double));
	double * yawD    = malloc(nTurb * sizeof(double));			// yaw angles corrected with wind direction
	int * front  =  malloc(nTurb * sizeof(int));	 // set of front turbines, not influenced by other turbines

   // now declared dynamically
	//double wakeDiameters[nTurb][nTurb][3]; // stores diameters of wake zone q of turbine i at downwind turbine j ([i][j][q])
	//double wakeCenters[nTurb][nTurb];    // stores wake center of turbine i at x-location of turbine j ([i],[j])
	//double wakeOverlap[nTurb][nTurb][3];	// stores area of wake zone q of turbine i that is overlapping with downwind turbine j ([i][j][q])
    double ***wakeDiameters; // stores diameters of wake zone q of turbine i at downwind turbine j ([i][j][q])
	double **wakeCenters;    // stores wake center of turbine i at x-location of turbine j ([i],[j])
	double ***wakeOverlap;	// stores area of wake zone q of turbine i that is overlapping with downwind turbine j ([i]



	////some variables temporarily used in calculation of wake properties
	double CTtilde;													// {\tilde C}_T, used to calculate wake deflection
	double OVdYd, OVr, OVR, OVz, OVL;								// distances, used in calculation of wake overlapping area
	double maxA, wakeEffCoeff,wakeEffCoeffPerZone,coeffScalingYaw;	// used in calculation of effective wind speed at dowstream turbines

	///Declare indices for for loops
	int turb, turbI, zone;

	//Initialize effU to effU_In
	for(turb = 0; turb<nTurb; turb++)
	{
		effU[turb] = effU_In[turb];
	}

	// find rotor areas
	for(turb = 0; turb<nTurb; turb++)
	{
		rotorA[turb] = 0.25*M_PI*rotorDiameter[turb]*rotorDiameter[turb];
	}

	// find down-wind,cross-wind coordinates
	for(turb = 0; turb<nTurb; turb++)
	{
		turbXd[turb] = X[turb]*cos(-effUdXY)-Y[turb]*sin(-effUdXY);
		turbYd[turb] = X[turb]*sin(-effUdXY)+Y[turb]*cos(-effUdXY);
		yawD[turb] = yaw[turb]; //  - effUdXY; //yaw is now defined relative to wind direction anyway
	};

	// allocate arrays to write wake geometry to
	// calculate the wake centers
	// calculate the wake diameters
    wakeDiameters =  malloc(nTurb * sizeof(double **));
	wakeCenters =    malloc(nTurb * sizeof(double *));
	for (turb = 0; turb<nTurb; turb++)
	{
      wakeDiameters[turb] =  malloc(nTurb * sizeof(double *));
		wakeCenters[turb] =  malloc(nTurb * sizeof(double));
		for (turbI = 0;turbI<nTurb; turbI++)
		{
         wakeDiameters[turb][turbI] =  malloc(3 * sizeof(double));
			CTtilde = pow(cos(yawD[turb]),2.0)*sin(yawD[turb])*2*axialInd[turb]*(1-axialInd[turb]);
			if(turbXd[turbI]>turbXd[turb])
			{

				//calculate the wake centers
				wakeCenters[turb][turbI] = turbYd[turb];
				wakeCenters[turb][turbI] = wakeCenters[turb][turbI]+CTtilde*(pow(CTtilde,2.0) + 15*pow((1 + 2*kd*(turbXd[turbI]-turbXd[turb])/rotorDiameter[turb]),4.0))/((30*kd/rotorDiameter[turb])*pow(1 + 2*kd*((turbXd[turbI]-turbXd[turb])/rotorDiameter[turb]),5.0));
				wakeCenters[turb][turbI] = wakeCenters[turb][turbI]-CTtilde*(pow(CTtilde,2.0) + 15)/((30*kd/rotorDiameter[turb]));
				wakeCenters[turb][turbI] = wakeCenters[turb][turbI]+ad+bd*(turbXd[turbI]-turbXd[turb]);
				//calculate the wake diameters
				wakeDiameters[turb][turbI][0] = fmax(rotorDiameter[turb]+2*ke*me1*(turbXd[turbI]-turbXd[turb]),0);
				wakeDiameters[turb][turbI][1] = fmax(rotorDiameter[turb]+2*ke*me2*(turbXd[turbI]-turbXd[turb]),0);
				wakeDiameters[turb][turbI][2] = fmax(rotorDiameter[turb]+2*ke*me3*(turbXd[turbI]-turbXd[turb]),0);
			};
		};
	};

	//calculate the wake overlapping areas
   wakeOverlap =  malloc(nTurb * sizeof(double **));
	for (turb = 0;turb<nTurb; turb++)
	{
      wakeOverlap[turb] =  malloc(nTurb * sizeof(double *));
		for (turbI = 0;turbI<nTurb; turbI++)
		{
         wakeOverlap[turb][turbI] =  malloc(3 * sizeof(double));
			OVdYd = wakeCenters[turb][turbI]-turbYd[turbI];
			OVr = rotorDiameter[turbI]/2;
			for (zone = 0;zone<3; zone++)
			{
				OVR = wakeDiameters[turb][turbI][zone]/2;
				OVdYd = fabs(OVdYd);
				if(OVdYd!=0)
				{
					OVL = (-pow(OVr,2.0)+pow(OVR,2.0)+pow(OVdYd,2.0))/(2.0*OVdYd);
				}
				else
				{
					OVL = 0;
				};

				OVz = pow(OVR,2.0)-pow(OVL,2.0);
				if(OVz>0)
				{
					OVz = sqrt(OVz);
				}
				else
				{
					OVz = 0;
				};

				if(OVdYd<(OVr+OVR))
				{
					if(OVL<OVR && (OVdYd-OVL)<OVr)
					{
						wakeOverlap[turb][turbI][zone] = pow(OVR,2.0)*acos(OVL/OVR) + pow(OVr,2.0)*acos((OVdYd-OVL)/OVr) - OVdYd*OVz;
					}
					else if (OVR>OVr)
					{
						wakeOverlap[turb][turbI][zone] = M_PI*pow(OVr,2.0);
					}
					else
					{
						wakeOverlap[turb][turbI][zone] = M_PI*pow(OVR,2.0);
					};
				}
				else
				{
					wakeOverlap[turb][turbI][zone] = 0;
				};
			};
		};
	};

	for (turb = 0;turb<nTurb; turb++)
	{
		for (turbI = 0;turbI<nTurb; turbI++)
		{
			wakeOverlap[turb][turbI][2] = wakeOverlap[turb][turbI][2]-wakeOverlap[turb][turbI][1];
			wakeOverlap[turb][turbI][1] = wakeOverlap[turb][turbI][1]-wakeOverlap[turb][turbI][0];
		};
	};

	// find the front turbines, not influenced by other turbines
	for (turbI = 0;turbI<nTurb; turbI++)
	{
		power[turbI] = 0.0;
		front[turbI] = 1;
		for (turb = 0;turb<nTurb; turb++)
		{
			if(wakeOverlap[turb][turbI][0]>0 || wakeOverlap[turb][turbI][1]>0 || wakeOverlap[turb][turbI][2]>0)
			{
				front[turbI] = 0;
			};
		};
		// find the effective wind speeds for the front turbines by inverting the Cp-power relation
		if(front[turbI]==1)
		{
			if(flagUseMeasPower==1)
			{
				effU[turbI] = pow(measPower[turbI]/(0.5*rho*rotorA[turbI]*Cp[turbI]*pow(cos(yawD[turbI]),pP)),1/3.0);
				power[turbI] = measPower[turbI];
			}
			else
			{
				power[turbI] = pow(effU[turbI],3.0) * (0.5*rho*rotorA[turbI]*Cp[turbI]*pow(cos(yawD[turbI]),pP));
			};

		};
	};

	// find effective wind speeds at downstream turbines, then predict power downstream turbine
	for (turbI = 0;turbI<nTurb; turbI++)
	{
		if(front[turbI] == 0)
		{
			wakeEffCoeff = 0;
			// 1) find the front turbine that a turbine has the largest overlap with, and use this as inflow speed
			maxA = 0;
			for (turb = 0;turb<nTurb; turb++)
			{
				if(front[turb] == 1)
				{
					if ((wakeOverlap[turb][turbI][0] + wakeOverlap[turb][turbI][1] + wakeOverlap[turb][turbI][2]) > maxA)
					{
						maxA = wakeOverlap[turb][turbI][0] + wakeOverlap[turb][turbI][1] + wakeOverlap[turb][turbI][2];
						effU[turbI] = effU[turb];
					};
				};
			};
			// 2) multiply the inflow speed with the wake coefficients to find effective wind speed at turbine;
			for (turb = 0;turb<nTurb; turb++)
			{
				wakeEffCoeffPerZone = 0;
				coeffScalingYaw = cos(aU*M_PI/180+bU*yawD[turb]);
				for (zone = 0;zone<3; zone++)
				{
					wakeOverlap[turb][turbI][zone] = fmin(wakeOverlap[turb][turbI][zone]/rotorA[turbI],1);
				};
				if (turbXd[turbI]>turbXd[turb])
				{
					wakeEffCoeffPerZone = wakeEffCoeffPerZone + pow((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke*(MU1/coeffScalingYaw)*(turbXd[turbI]-turbXd[turb])),2.0)*wakeOverlap[turb][turbI][0];
					wakeEffCoeffPerZone = wakeEffCoeffPerZone + pow((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke*(MU2/coeffScalingYaw)*(turbXd[turbI]-turbXd[turb])),2.0)*wakeOverlap[turb][turbI][1];
					wakeEffCoeffPerZone = wakeEffCoeffPerZone + pow((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke*(MU3/coeffScalingYaw)*(turbXd[turbI]-turbXd[turb])),2.0)*wakeOverlap[turb][turbI][2];
				};
				wakeEffCoeff = wakeEffCoeff+pow(axialInd[turb]*wakeEffCoeffPerZone,2.0);
			};
			effU[turbI] = effU[turbI] * (1 - 2 * sqrt(wakeEffCoeff));
			power[turbI] = pow(effU[turbI],3.0) * (0.5*rho*rotorA[turbI]*Cp[turbI]*pow(cos(yawD[turbI]),pP));
		};
	};

	for (turb = 0; turb<nTurb; turb++){
		for (turbI = 0; turbI<nTurb; turbI++){
			free(wakeDiameters[turb][turbI]);
			free(wakeOverlap[turb][turbI]);
		}
		free(wakeDiameters[turb]);
		free(wakeCenters[turb]);
		free(wakeOverlap[turb]);
	}
	free(wakeDiameters);
	free(wakeCenters);
	free(wakeOverlap);


	free(rotorA);
	free(turbXd);
	free(turbYd);
	free(yawD);
	free(front);
}

