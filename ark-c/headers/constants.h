#ifndef _CONSTANTS_H
#define _CONSTANTS_H

//	-----------------------------------------------------------------------------------
//	sound	-	sound speed
//	ro0_g	-	unperturbed density of the liquid
//	ro0_s	-	unperturbed density of the barriers material
//	dt		-	time step
//	VIS		-	kinematic viscosity
//	u10		-	initial speed along the axis x1
//	u20		-	initial speed along the axis X2
//	u30		-	initial speed along the axis X3
//	t0		-	initial temperature
//	TIME	-	current time
//	CFL		-	Courant number
//	pOutlet	-	pressure on the top border
//	u3Inlet	-	speed on the bottom border along the X3 axis
//	u2Inlet -	speed on the bottom border along the X2 axis
//	u1Inlet	-	speed on the bottom border along the x1 axis
//	tInlet	-	temperature on the bottom border

double sound, ro0_g, ro0_s, dt, VIS, u10, u20, u30, t0, TIME, CFL, pOutlet;
double u3Inlet, u2Inlet, u1Inlet, tInlet;

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

#endif
