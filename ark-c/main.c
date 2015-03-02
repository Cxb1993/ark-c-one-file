#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "headers/forces.h"
#include "headers/bulk.h"
#include "headers/constants.h"
#include "headers/stress.h"

void Input();
void InitializeData();
void TimeStepSize();
void Phase1();
void Phase2();
void StressTensor();
void UseForces();
void WriteDataTecplot();
void WriteDataParaView();
void FreeMemory();

double ***allocate3D(int, int, int);
void deallocate3D(double***, int, int);

double min3d(double, double, double);
double max3d(double, double, double);

int main (int argc, char** argv)
{
	Input();
	InitializeData();
	
	//first output
	WriteDataParaView();
	//WriteDataTecplot();
	
	nStep = 0;
	TIME = 0.0;
	
	do
	{
		TimeStepSize();
		nStep++;
		TIME += 0.5*dt;
		Phase1();
		StressTensor();
		UseForces();
		Phase2();
		Phase1();
		StressTensor();
		UseForces();
		TIME += 0.5*dt;

		if (nStep % nPrint == 0) 
		{
			WriteDataParaView();
			//WriteDataTecplot();
			printf("step: %d dt:%E\n", nStep, dt);
		}
	} while (nStep < nStop);
	FreeMemory();
    return 0;
}

void Input()
{
	// geometry index
	l = 1;

	// total number of grid nodes along the x1 axis
	n1_g = 16;
	// total number of grid nodes along the X2 axis
	n2_g = 16;
	// total number of grid nodes along the X3 axis
	n3_g = 32;

	// number of grid nodes along the x1 axis to 1 processor
	n1 = n1_g/1;
	// number of grid nodes along the X2 axis to 1 processor
	n2 = n2_g/1;
	// number of grid nodes along the X3 axis to 1 processor
	n3 = n3_g/1;

	// coordinates of west plane
	x1_w = 0.;
	// coordinates of east plane
	x1_e = 2.;

	// coordinates of south plane
	x2_s = 0.;
	// coordinates of north plane
	x2_n = 2. * M_PI;

	// coordinates of bottom plane
	x3_b = 0.;
	// coordinates of top plane
	x3_t = 4. * M_PI;

	// total number of steps
	nStop = 30000;
	// print interval
	nPrint = 100;

	// Courant number
	CFL = 0.2;

	// kinematic viscosity
	VIS = 0.5/50.;  // 0.5
	// initial temperature
	t0 = 1.;

	// initial velocity along the x1 axis
	u10 = 0.;
	// initial velocity along the X2 axis
	u20 = 0.;
	// initial velocity along the X3 axis
	u30 = 1.;
	// sound velocity
	sound = 10.;

	// pressure on the top plane
	pOutlet = 0;
	// velocity on the bottom plane along the X3 axis
	u3Inlet = u30;
	// velocity on the bottom plane along the X2 axis
	u2Inlet = u20;
	// velocity on the bottom plane along the x1 axis
	u1Inlet = u10;
	// temperature on the bottom plane
	tInlet = t0;
	// unperturbed density of the liquid
	ro0_g = 1.;
	// unperturbed density of the borders material
	ro0_s = 100000000.;

	// #####################################################
	// 				block of arrays allocation
	// #####################################################

	// coordinates of grid nodes along the all of the axises
	x1 = calloc(n1 + 2, sizeof(double));
	x2 = calloc(n2 + 2, sizeof(double));
	x3 = calloc(n3 + 2, sizeof(double));

	// variables on the current time step
	roCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	tCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	u1Con = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	u2Con = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	u3Con = allocate3D(n2 + 1, n1 + 1, n3 + 1);

	// variables on the next time step
	ronCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	tnCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	u1nCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	u2nCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	u3nCon = allocate3D(n2 + 1, n1 + 1, n3 + 1);

	// forces
	f1 = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	f2 = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	f3 = allocate3D(n2 + 1, n1 + 1, n3 + 1);
	
	// variables perpendicular to the axis x1
	ro1 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	t1 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u11 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u21 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u31 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	p1 = allocate3D(n2 + 2, n1 + 2, n3 + 2);

	// variables perpendicular to the axis X2
	ro2 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	t2 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u12 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u22 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u32 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	p2 = allocate3D(n2 + 2, n1 + 2, n3 + 2);

	// variables perpendicular to the axis X3
	ro3 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	t3 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u13 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u23 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	u33 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	p3 = allocate3D(n2 + 2, n1 + 2, n3 + 2);

	// NMAX = MAX(n1, n2, n3)
	int nmax = (int)max3d(n1, n2, n3);

	// additional buffers for phase 2
	rBuf	= calloc((nmax+2), sizeof(double));
	qBuf	= calloc((nmax+1), sizeof(double));
	tfBuf	= calloc((nmax+2), sizeof(double));
	tbBuf	= calloc((nmax+2), sizeof(double));
	u2fBuf	= calloc((nmax+2), sizeof(double));
	u2bBuf	= calloc((nmax+2), sizeof(double));
	u3fBuf	= calloc((nmax+2), sizeof(double));
	u3bBuf	= calloc((nmax+2), sizeof(double));

	// friction stress
    sigm11 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	sigm21 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	sigm31 = allocate3D(n2 + 2, n1 + 2, n3 + 2);

	sigm12 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	sigm22 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	sigm32 = allocate3D(n2 + 2, n1 + 2, n3 + 2);

	sigm13 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	sigm23 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
	sigm33 = allocate3D(n2 + 2, n1 + 2, n3 + 2);
}

void InitializeData()
{
//	double	ALFA0 = 0.204,
//			BETA = 0.3,
//			R0 = 0.05;

	// grid step along the x1 axis
	dx1=(x1_e-x1_w)/(n1-1);
	dx2=(x2_n-x2_s)/(n2-1);
	dx3=(x3_t-x3_b)/(n3-1);

	x1[0] = x1_w - dx1;
	x2[0] = x2_s - dx2;
	x3[0] = x3_b - dx3;

	// #####################################################
	//			block of arrays initialization
	// #####################################################

	// along the X1 axis
	for (int j = 0; j < n1 + 1 ; ++j) {
		x1[j + 1] = x1[j] + dx1;
	}

	// along the X2 axis
	for (int i = 0; i < n2 + 1; ++i) {
		x2[i + 1] = x2[i] + dx2;
	}

	// along the X3 axis
	for (int k = 0; k < n3 + 1; ++k) {
		x3[k + 1] = x3[k] + dx3;
	}

	for (int i = 0; i < n2 + 1; ++i) {
		for (int j = 0; j < n1 + 1; ++j) {
			for (int k = 0; k < n3 + 1; ++k) {
				u1Con[i][j][k] = u1nCon[i][j][k] = u10;
				u2Con[i][j][k] = u2nCon[i][j][k] = u20;
				u3Con[i][j][k] = u3nCon[i][j][k] = u30;
				roCon[i][j][k] = ronCon[i][j][k] = ro0_g;
				tCon[i][j][k] = tnCon[i][j][k] = t0;
			}
		}
	}

	for (int i = 0; i < n2 + 2; ++i) {
		for (int j = 0; j < n1 + 2; ++j) {
			for (int k = 0; k < n3 + 2; ++k) {
				p1[i][j][k] = p2[i][j][k] = p3[i][j][k] = 0.;
				ro1[i][j][k] = ro2[i][j][k] = ro3[i][j][k] = ro0_g;
				u11[i][j][k] = u12[i][j][k] = u13[i][j][k] = u10;
				u21[i][j][k] = u22[i][j][k] = u23[i][j][k] = u20;
				u31[i][j][k] = u32[i][j][k] = u33[i][j][k] = u30;
				t1[i][j][k] = t2[i][j][k] = t3[i][j][k] = t0;
			}
		}
	}
}

void TimeStepSize()
{
	double x1c, x1l, u1_c, u2_c, u3_c, dtu1, dtu2, dtu3, dtu, dtv1, dtv2, dtv3, dtv;

	dt = pow(10, 8);

	for (int i = 1; i < n2; ++i) {
		for (int j = 1; j < n1; ++j) {
			for (int k = 1; k < n3; ++k) {
				x1c = (x1[j+1] + x1[j])/2;
				x1l = 1+(l-1)*(x1c-1);

				u1_c = u1Con[i][j][k];
				u2_c = u2Con[i][j][k];
				u3_c = u3Con[i][j][k];

				dtu1 = CFL*dx1/(sound + fabs(u1_c));
				dtu2 = CFL*x1l*dx2/(sound + fabs(u2_c));
				dtu3 = CFL*dx3/(sound + fabs(u3_c));

				// DTU = MIN(DTU1, DTU2, DTU3)
				dtu = min3d(dtu1, dtu2, dtu3);

				if (VIS > pow(10, -16)) {
					dtv1 = CFL*dx1*dx1/(2.*VIS);
					dtv2 = CFL*(x1l*dx2)*(x1l*dx2)/(2.*VIS);
					dtv3 = CFL*dx3*dx3/(2.*VIS);

					// DTV = MIN (DTV1, DTV2, DTV3)
					dtv = min3d(dtv1, dtv2, dtv3);
				} else {
					dtv = pow(10, 16);
				}

				// DT = MIN(DT, DTU, DTV)
				dt = min3d(dt, dtu, dtv);
			}
		}
	}
}

void Phase1()
{
	// initialization
	for (int i = 1; i < n2; i++)
	{
		for (int j = 1; j < n1; j++)
		{
			for (int k = 1; k < n3; k++)
			{
				u1Con[i][j][k] = u1nCon[i][j][k];
				u2Con[i][j][k] = u2nCon[i][j][k];
				u3Con[i][j][k] = u3nCon[i][j][k];
				roCon[i][j][k] = ronCon[i][j][k];
				tCon[i][j][k] = tnCon[i][j][k];
			}
		}
	}

	// geometric characteristics of the computational cell
	double ds1, ds2, ds3, dvc, xle, xlw, xlt, xlb, alfa;

			// velocity, density, temperature and pressure on the eastern plane
	double	u1_e, u2_e, u3_e, ro_e, t_e, p_e,
			// velocity, density , temperature and pressure on the western plane
			u1_w, u2_w, u3_w, ro_w, t_w, p_w,
			// velocity, density , temperature and pressure on the northern plane
			u1_n, u2_n, u3_n, ro_n, t_n, p_n,
			// velocity, density , temperature and pressure on the southern plane
			u1_s, u2_s, u3_s, ro_s, t_s, p_s,
			// velocity, density , temperature and pressure on the top plane
			u1_t, u2_t, u3_t, ro_t, t_t, p_t,
			// velocity, density , temperature and pressure on the bottom plane
			u1_b, u2_b, u3_b, ro_b, t_b, p_b,
			// velocity, density and temperature in the cell center
			u1_c, u2_c, u3_c, ro_c, t_c,
			// velocity, density and temperature in the cell center on the next time step
			u1_cn, u2_cn, u3_cn, ro_cn, t_cn;

	// plane squares
	ds1 = dx2*dx3;
	ds2 = dx1*dx3; 
	ds3 = dx1*dx2;

	for (int i = 1; i < n2; i++)
	{
		for (int j = 1; j < n1; j++)
		{
			for (int k = 1; k < n3; k++)
			{
				// geometric factor of cylindricity
				xle = 1 + (l - 1)*(x1[j] - 1);
				// geometric factor of cylindricity
				xlw = 1 + (l - 1)*(x1[j - 1] - 1);
				// geometric factor of cylindricity
				xlt = 0.5*(xle + xlw);
				// geometric factor of cylindricity
				xlb = xlt;
				// cell volume
				dvc = 0.5*(xle + xlw)*dx1*dx2*dx3;
				// geometric parameter
				alfa = dt*(l - 1)*dx1*dx2*dx3 / (2 * dvc);

				// #########################################################
				//	get velocity, density , temperature and pressure values
				// #########################################################

				// east plane
				u1_e = u11[i][j + 1][k];
				u2_e = u21[i][j + 1][k];
				u3_e = u31[i][j + 1][k];
				ro_e = ro1[i][j + 1][k];
				t_e = t1[i][j + 1][k];
				p_e = p1[i][j + 1][k];
				// west plane
				u1_w = u11[i][j][k];
				u2_w = u21[i][j][k];
				u3_w = u31[i][j][k];
				ro_w = ro1[i][j][k];
				t_w = t1[i][j][k];
				p_w = p1[i][j][k];

				// north plane
				u1_n = u12[i + 1][j][k];
				u2_n = u22[i + 1][j][k];
				u3_n = u32[i + 1][j][k];
				ro_n = ro2[i + 1][j][k];
				t_n = t2[i + 1][j][k];
				p_n = p2[i + 1][j][k];
				// south plane
				u1_s = u12[i][j][k];
				u2_s = u22[i][j][k];
				u3_s = u32[i][j][k];
				ro_s = ro2[i][j][k];
				t_s = t2[i][j][k];
				p_s = p2[i][j][k];

				// top plane
				u1_t = u13[i][j][k + 1];
				u2_t = u23[i][j][k + 1];
				u3_t = u33[i][j][k + 1];
				ro_t = ro3[i][j][k + 1];
				t_t = t3[i][j][k + 1];
				p_t = p3[i][j][k + 1];
				// bottom plane
				u1_b = u13[i][j][k];
				u2_b = u23[i][j][k];
				u3_b = u33[i][j][k];
				ro_b = ro3[i][j][k];
				t_b = t3[i][j][k];
				p_b = p3[i][j][k];

				// cell center
				u1_c = u1Con[i][j][k];
				u2_c = u2Con[i][j][k];
				u3_c = u3Con[i][j][k];
				ro_c = roCon[i][j][k];
				t_c = tCon[i][j][k];
				
				// #####################################################
				//				new values evaluating
				// #####################################################

				// new density
				ro_cn = (ro_c*dvc - 0.5*dt*((xle*ro_e*u1_e - xlw*ro_w*u1_w)*ds1 +
					(ro_n*u2_n - ro_s*u2_s)*ds2 + (xlt*ro_t*u3_t - xlb*ro_b*u3_b)*ds3)) / dvc;
				
				// new conservative velocity along the x1 axis
				double u1_cp = (ro_c*u1_c*dvc - 0.5*dt*((ro_n*u2_n*u1_n - ro_s*u2_s*u1_s)*ds2 +
					(xlt*ro_t*u3_t*u1_t - xlb*ro_b*u3_b*u1_b)*ds3 +
					(xle*(ro_e*u1_e*u1_e + p_e) - xlw*(ro_w*u1_w*u1_w + p_w))*ds1
					- 0.5*(l - 1)*(p_e + p_w)*dx1*dx2*dx3)) / (dvc*ro_cn);
				// new conservative velocity along the X2 axis
				double u2_cp = (ro_c*u2_c*dvc - 0.5*dt*((xle*ro_e*u1_e*u2_e - xlw*ro_w*u1_w*u2_w)*ds1 +
					((ro_n*u2_n*u2_n + p_n) - (ro_s*u2_s*u2_s + p_s))*ds2 +
					(xlt*ro_t*u3_t*u2_t - xlb*ro_b*u3_b*u2_b)*ds3)) / (dvc*ro_cn);

				// take into account of centrifugal and Coriolis forces
				u1_cn = (u1_cp - alfa*u2_c*u1_cp) / (1 + (alfa*u2_c)*(alfa*u2_c));
				u2_cn = u2_cp - alfa*u2_c*u1_cn;

				// new conservative velocity along the X3 axis
				u3_cn = (ro_c*u3_c*dvc - 0.5*dt*((xle*ro_e*u1_e*u3_e - xlw*ro_w*u1_w*u3_w)*ds1 +
					(ro_n*u2_n*u3_n - ro_s*u2_s*u3_s)*ds2 +
					(xlt*(ro_t*u3_t*u3_t + p_t) - xlb*(ro_b*u3_b*u3_b + p_b))*ds3)) / (dvc*ro_cn);

				// new temperature
				t_cn = (ro_c*t_c*dvc - 0.5*dt*((xle*ro_e*t_e*u1_e - xlw*ro_w*t_w*u1_w)*ds1 +
					(ro_n*t_n*u2_n - ro_s*t_s*u2_s)*ds2 +
					(xlt*ro_t*t_t*u3_t - xlb*ro_b*t_b*u3_b)*ds3)) / (dvc*ro_cn);

				// finally
				u1nCon[i][j][k] = u1_cn;
				u2nCon[i][j][k] = u2_cn;
				u3nCon[i][j][k] = u3_cn;
				ronCon[i][j][k] = ro_cn;
				tnCon[i][j][k] = t_cn;
			}
		}
	}

	// #####################################################
	// 					boudary conditions
	// #####################################################

	// periodicity along the X3 axis
	for (int i = 1; i < n2; ++i)
	{
		for (int j = 1; j < n1; ++j)
		{
			// periodicity conditions on the north plane
			u1nCon[i][j][n3] = u1nCon[i][j][1];
			u2nCon[i][j][n3] = u2nCon[i][j][1];
			u3nCon[i][j][n3] = u3nCon[i][j][1];
			ronCon[i][j][n3] = ronCon[i][j][1];
			tnCon[i][j][n3] = tnCon[i][j][1];

			// periodicity conditions on the south plane
			u1nCon[i][j][0] = u1nCon[i][j][n3 - 1];
			u2nCon[i][j][0] = u2nCon[i][j][n3 - 1];
			u3nCon[i][j][0] = u3nCon[i][j][n3 - 1];
			ronCon[i][j][0] = ronCon[i][j][n3 - 1];
			tnCon[i][j][0] = tnCon[i][j][n3 - 1];
		}
	}

	// periodicity along the X2 axis
	for (int j = 1; j < n1; j++)
	{
		for (int k = 1; k < n3; k++)
		{
			// periodicity conditions on the north plane
			u1nCon[n2][j][k] = u1nCon[1][j][k];
			u2nCon[n2][j][k] = u2nCon[1][j][k];
			u3nCon[n2][j][k] = u3nCon[1][j][k];
			ronCon[n2][j][k] = ronCon[1][j][k];
			tnCon[n2][j][k] = tnCon[1][j][k];

			// periodicity conditions on the south plane
			u1nCon[0][j][k] = u1nCon[n2 - 1][j][k];
			u2nCon[0][j][k] = u2nCon[n2 - 1][j][k];
			u3nCon[0][j][k] = u3nCon[n2 - 1][j][k];
			ronCon[0][j][k] = ronCon[n2 - 1][j][k];
			tnCon[0][j][k] = tnCon[n2 - 1][j][k];
		}		
	}
	
	// no-slip conditions along the X1 axis
	for (int i = 1; i < n2; ++i)
	{
		for (int k = 1; k < n3; ++k)
		{
			// no-slip conditions on the west plane
			u1nCon[i][0][k] = 0.;
			u2nCon[i][0][k] = 0.;
			u3nCon[i][0][k] = 0.;
			ronCon[i][0][k] = ronCon[i][1][k];
			tnCon[i][0][k] = t0;

			// no-slip conditions on the east plane
			u1nCon[i][n1][k] = 0.;
			u2nCon[i][n1][k] = 0.;
			u3nCon[i][n1][k] = 0.;
			ronCon[i][n1][k] = ronCon[i][n1 - 1][k];
			tnCon[i][n1][k] = t0;
		}
	}
}

void StressTensor()
{
	// initialization of friction stress arrays
	for (int i = 0; i < n2 + 2; ++i) {
		for (int j = 0; j < n1 + 2; ++j) {
			for (int k = 0; k < n3 + 2; ++k) {
				sigm11[i][j][k] = 0.0;
				sigm21[i][j][k] = 0.0;
				sigm31[i][j][k] = 0.0;

				sigm12[i][j][k] = 0.0;
				sigm22[i][j][k] = 0.0;
				sigm32[i][j][k] = 0.0;

				sigm13[i][j][k] = 0.0;
				sigm23[i][j][k] = 0.0;
				sigm33[i][j][k] = 0.0;
			}
		}
	}


	// #####################################################
	// 					boudary conditions
	// #####################################################

	// periodicity along the X3 axis
	for (int i = 1; i < n2; ++i)
	{
		for (int j = 1; j < n1; ++j)
		{
			// periodicity conditions on the north plane
			u1Con[i][j][n3] = u1Con[i][j][1];
			u2Con[i][j][n3] = u2Con[i][j][1];
			u3Con[i][j][n3] = u3Con[i][j][1];
			roCon[i][j][n3] = roCon[i][j][1];
			tCon[i][j][n3] = tCon[i][j][1];

			// periodicity conditions on the south plane
			u1Con[i][j][0] = u1Con[i][j][n3 - 1];
			u2Con[i][j][0] = u2Con[i][j][n3 - 1];
			u3Con[i][j][0] = u3Con[i][j][n3 - 1];
			roCon[i][j][0] = roCon[i][j][n3 - 1];
			tCon[i][j][0] = tCon[i][j][n3 - 1];
		}
	}

	// periodicity along the X2 axis
	for (int j = 1; j < n1; j++)
	{
		for (int k = 1; k < n3; k++)
		{
			// periodicity conditions on the north plane
			u1Con[n2][j][k] = u1nCon[1][j][k];
			u2Con[n2][j][k] = u2nCon[1][j][k];
			u3Con[n2][j][k] = u3nCon[1][j][k];
			roCon[n2][j][k] = ronCon[1][j][k];
			tCon[n2][j][k] = tnCon[1][j][k];

			// periodicity conditions on the south plane
			u1Con[0][j][k] = u1Con[n2 - 1][j][k];
			u2Con[0][j][k] = u2Con[n2 - 1][j][k];
			u3Con[0][j][k] = u3Con[n2 - 1][j][k];
			roCon[0][j][k] = roCon[n2 - 1][j][k];
			tCon[0][j][k] = tCon[n2 - 1][j][k];
		}
	}

	// no-slip conditions along the X1 axis
	for (int i = 1; i < n2; ++i)
	{
		for (int k = 1; k < n3; ++k)
		{
			// no-slip conditions on the west plane
			u1Con[i][0][k] = 0.;
			u2Con[i][0][k] = 0.;
			u3Con[i][0][k] = 0.;
			roCon[i][0][k] = roCon[i][1][k];
			tCon[i][0][k] = t0;
			
			// no-slip conditions on the east plane
			u1Con[i][n1][k] = 0.;
			u2Con[i][n1][k] = 0.;
			u3Con[i][n1][k] = 0.;
			roCon[i][n1][k] = roCon[i][n1 - 1][k];
			tCon[i][n1][k] = t0;
		}
	}

	// #####################################################
	// 				bypassing along the faces
	// #####################################################

	double xle, xlw, xlt, xln, u1_c, u1_ce, u2_c, u2_ce, u3_c, u3_ce;

	// bypassing along the face perpendicular to X1
	for (int i = 1; i <= n2; ++i) {
		for (int k = 1; k <= n3; ++k) {
			for (int j = 1; j <= n1; ++j) {
				// geometric characteristics of the computational cell
				// geometric factor of cylindricity
				xle = 1 + (l-1)*(x1[j] - 1);

				// velocity components in cell centers
				u1_c = u1Con[i][j][k];
				u1_ce = u1Con[i][j - 1][k];

				u2_c = u2Con[i][j][k];
				u2_ce = u2Con[i][j - 1][k];

				u3_c = u3Con[i][j][k];
				u3_ce = u3Con[i][j - 1][k];

				// friction stress
				sigm11[i][j][k]=-VIS*xle*(u1_ce-u1_c)/dx1;
				sigm21[i][j][k]=-VIS*xle*(u2_ce-u2_c)/dx1;
				sigm31[i][j][k]=-VIS*xle*(u3_ce-u3_c)/dx1;
			}
		}
	}

	double u1_cn, u2_cn, u3_cn;

	// bypassing along the face perpenditcular to X2
	for (int j = 1; j <= n1; ++j) {
		for (int k = 1; k <= n3; ++k) {
			for (int i = 1; i <= n2; ++i) {
				// geometric characteristics of the computational cell
				// geometric factor of cylindricity
				xle = 1 + (l-1)*(x1[j] - 1);
				// geometric factor of cylindricity
				xlw = 1 + (l-1)*(x1[j-1] - 1);
				// geometric factor of cylindricity
				xlt = 0.5*(xle + xlw);
				// geometric factor of cylindricity
				xln = xlt;

				// velocity components in cell centers
				u1_c = u1Con[i][j][k];
				u1_cn = u1Con[i - 1][j][k];

				u2_c = u2Con[i][j][k];
				u2_cn = u2Con[i - 1][j][k];

				u3_c = u3Con[i][j][k];
				u3_cn = u3Con[i - 1][j][k];

				// friction stress
				sigm12[i][j][k]=-VIS*((u1_cn-u1_c)/dx2 -(l-1)*(u2_c+u2_cn))/xln;
				sigm22[i][j][k]=-VIS*((u2_cn-u2_c)/dx2 +(l-1)*(u1_c+u1_cn))/xln;
				sigm32[i][j][k]=-VIS*(u3_cn-u3_c)/dx2;
			}
		}
	}

	double u1_ct, u2_ct, u3_ct;

	// bypassing along the face perpenditcular to X3
	for (int i = 1; i <= n2; ++i) {
		for (int j = 1; j <= n1; ++j) {
			for (int k = 1; k <= n3; ++k) {
				// geometric characteristics of the computational cell
				// geometric factor of cylindricity
				xle = 1 + (l-1)*(x1[j] - 1);
				// geometric factor of cylindricity
				xlw = 1 + (l-1)*(x1[j-1] - 1);
				// geometric factor of cylindricity
				xlt = 0.5*(xle + xlw);

				// velocity components in the cell centers
				u1_c = u1Con[i][j][k];
				u1_ct = u1Con[i][j][k-1];

				u2_c = u2Con[i][j][k];
				u2_ct = u2Con[i][j][k-1];

				u3_c = u3Con[i][j][k];
				u3_ct = u3Con[i][j][k-1];

				// friction stress
				sigm13[i][j][k]=-VIS*xlt*(u1_ct-u1_c)/dx3;
				sigm23[i][j][k]=-VIS*xlt*(u2_ct-u2_c)/dx3;
				sigm33[i][j][k]=-VIS*xlt*(u3_ct-u3_c)/dx3;
			}
		}
	}

	// #####################################################
	// 				friction forces computation
	// #####################################################

	double ds1, ds2, ds3, sigm1c, sigm2c;

	// area of the face perpendicuar to x1
	ds1 = dx2*dx3;
	ds2 = dx1*dx3;
	ds3 = dx1*dx2;

	for (int i = 1; i < n2; ++i) {
		for (int j = 1; j < n1; ++j) {
			for (int k = 1; k < n3; ++k) {
				// geometric characteristics of the computational cell
				// geometric factor of cylindricity
				xle = 1 + (l-1)*(x1[j] - 1);
				// geometric factor of cylindricity
				xlw = 1 + (l-1)*(x1[j-1] - 1);
				// geometric factor of cylindricity
				xlt = 0.5*(xle + xlw);

				sigm1c = VIS*u1Con[i][j][k]/xlt;
				sigm2c = VIS*u2Con[i][j][k]/xlt;

				// friction forces
				f1[i][j][k] =
					(sigm11[i][j + 1][k] - sigm11[i][j][k]) * ds1 +
					(sigm12[i + 1][j][k] - sigm12[i][j][k]) * ds2 +
					(sigm13[i][j][k + 1] - sigm13[i][j][k]) * ds3 -
					(l - 1)*sigm1c*dx1*dx2*dx3;

				f2[i][j][k] =
					(sigm21[i][j + 1][k] - sigm21[i][j][k]) * ds1 +
					(sigm22[i + 1][j][k] - sigm22[i][j][k]) * ds2 +
					(sigm23[i][j][k + 1] - sigm23[i][j][k]) * ds3 -
					(l - 1)*sigm2c*dx1*dx2*dx3;

				f3[i][j][k] =
					(sigm31[i][j + 1][k] - sigm31[i][j][k]) * ds1 +
					(sigm32[i + 1][j][k] - sigm32[i][j][k]) * ds2 +
					(sigm33[i][j][k + 1] - sigm33[i][j][k]) * ds3;
			}
		}
	}
}

void UseForces()
{
	double xle, xlw, dvc, ro_c, ro_cn;

	for (int i = 1; i < n2; ++i) {
		for (int j = 1; j < n1; ++j) {
			for (int k = 1; k < n3; ++k) {
				// geometric factor of cylindricity
				xle = 1+(l-1)*(x1[j+1]-1);
				xlw = 1+(l-1)*(x1[j]-1);
				// cell volume
				dvc = 0.5*(xle+xlw)*dx1*dx2*dx3;

				ro_c = roCon[i][j][k];
				ro_cn = ronCon[i][j][k];

				u1nCon[i][j][k] = (ro_c*dvc*u1nCon[i][j][k] + 0.5*dt*f1[i][j][k])/(dvc*ro_cn);
				u2nCon[i][j][k] = (ro_c*dvc*u2nCon[i][j][k] + 0.5*dt*f2[i][j][k])/(dvc*ro_cn);
				u3nCon[i][j][k] = (ro_c*dvc*u3nCon[i][j][k] + 0.5*dt*f3[i][j][k])/(dvc*ro_cn);
			}
		}
	}
}

void Phase2()
{
	double u1_f, u1_b, u1_cn, u1_c,
		u2_f, u2_fn, u2_b, u2_bn, u2_cn, u2_c,
		u3_f, u3_fn, u3_b, u3_bn, u3_cn, u3_c,
		ro_f, ro_b, ro_cn, ro_c,
		t_f, t_fn, t_b, t_bn, t_cn, t_c,
		p_f, p_b, p_cn, p_c,
		r_f, r_fn, r_b, r_cn, r_c,
		q_f, q_b, q_bn, q_cn, q_c;

	double gr, gt, gu2, gu3, gq;

	double rmax, rmin, qmax, qmin, tmax, tmin, u2_max, u2_min, u3_max, u3_min;

	double ro0_b, ro0_f, qn, pn, rn, ro_n, tn, un, u2_n, u3_n, ucf, ucb;

	// first local invariants for the interior faces puts in the buffer arrays, bypass on the center of the cell
	// then by taking into account the boundary condition calculates extreme elements of the buffers
	// and only then calculates the flow variables

	// flow variables calculation on DS1 faces orthogonal X1 axis 

	// bypassing along the X1 axis

	// only interior faces !
	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n2; i++)
		{
			for (int j = 1; j < n1; j++)
			{
				u1_f = u11[i][j + 1][k];
				u1_b = u11[i][j][k];
				u1_cn = u1nCon[i][j][k];
				u1_c = u1Con[i][j][k];

				u2_f = u21[i][j + 1][k];
				u2_b = u21[i][j][k];
				u2_cn = u2nCon[i][j][k];
				u2_c = u2Con[i][j][k];

				u3_f = u31[i][j + 1][k];
				u3_b = u31[i][j][k];
				u3_cn = u3nCon[i][j][k];
				u3_c = u3Con[i][j][k];

				ro_f = ro1[i][j + 1][k];
				ro_b = ro1[i][j][k];
				ro_cn = ronCon[i][j][k];
				ro_c = roCon[i][j][k];

				t_f = t1[i][j + 1][k];
				t_b = t1[i][j][k];
				t_cn = tnCon[i][j][k];
				t_c = tCon[i][j][k];

				p_f = p1[i][j + 1][k];
				p_b = p1[i][j][k];
				p_cn = sound*sound*(ro_cn - ro0_g);
				p_c = sound*sound*(ro_c - ro0_g);

				// invariant calculation

				r_f = u1_f + p_f / (ro0_g * sound);
				r_b = u1_b + p_b / (ro0_g * sound);
				r_cn = u1_cn + p_cn / (ro0_g * sound);
				r_c = u1_c + p_c / (ro0_g * sound);

				r_fn = 2 * r_cn - r_b;

				q_f = u1_f - p_f / (ro0_g * sound);
				q_b = u1_b - p_b / (ro0_g * sound);
				q_cn = u1_cn - p_cn / (ro0_g * sound);
				q_c = u1_c - p_c / (ro0_g * sound);

				q_bn = 2 * q_cn - q_f;

				t_fn = 2 * t_cn - t_b;
				t_bn = 2 * t_cn - t_f;

				u2_fn = 2 * u2_cn - u2_b;
				u2_bn = 2 * u2_cn - u2_f;

				u3_fn = 2 * u3_cn - u3_b;
				u3_bn = 2 * u3_cn - u3_f;

				// the permissible range of changes
				gr = 2 * (r_cn - r_c) / dt + (u1_cn + sound)*(r_f - r_b) / dx1;
				gq = 2 * (r_cn - r_c) / dt + (u1_cn - sound)*(q_f - q_b) / dx1;

				gt = 2 * (t_cn - t_c) / dt + u1_cn*(t_f - t_b) / dx1;
				gu2 = 2 * (u2_cn - u2_c) / dt + u1_cn*(u2_f - u2_b) / dx1;
				gu3 = 2 * (u3_cn - u3_c) / dt + u1_cn*(u3_f - u3_b) / dx1;

				// RMAX=MAX(RF,RC,RB) +dt*GR
				rmax = max3d(r_f, r_c, r_b) + dt*gr;

				// RMIN=MIN(RF,RC,RB) +dt*GR
				rmin = min3d(r_f, r_c, r_b) + dt*gr;

				// QMAX=MAX(QF,QC,QB) +dt*GQ
				qmax = max3d(q_f, q_c, q_b) + dt*gq;

				// QMIN=MIN(QF,QC,QB) +dt*GQ
				qmin = min3d(q_f, q_c, q_b) + dt*gq;

				// TMAX=MAX(TF,TC,TB) +dt*GT
				tmax = max3d(t_f, t_c, t_b) + dt*gt;

				// TMIN=MIN(TF,TC,TB) +dt*GT
				tmin = min3d(t_f, t_c, t_b) + dt*gt;

				// U2MAX=MAX(U2F,U2C,U2B) +dt*GU2
				u2_max = max3d(u2_f, u2_c, u2_b) + dt*gu2;

				// U2MIN=MIN(U2F,U2C,U2B) +dt*GU2
				u2_min = min3d(u2_f, u2_c, u2_b) + dt*gu2;

				// U3MAX=MAX(U3F,U3C,U3B) +dt*GU3
				u3_max = max3d(u3_f, u3_c, u3_b) + dt*gu3;

				// U3MIN=MIN(U3F,U3C,U3B) +dt*GU3 
				u3_min = min3d(u3_f, u3_c, u3_b) + dt*gu3;

				// invariants correction
				if (r_fn > rmax) r_fn = rmax;
				if (r_fn < rmin) r_fn = rmin;

				if (q_bn > qmax) q_bn = qmax;
				if (q_bn < qmin) q_bn = qmin;

				if (t_fn > tmax) t_fn = tmax;
				if (t_fn < tmin) t_fn = tmin;

				if (t_bn > tmax) t_bn = tmax;
				if (t_bn < tmin) t_bn = tmin;

				if (u2_fn > u2_max) u2_fn = u2_max;
				if (u2_fn < u2_min) u2_fn = u2_min;

				if (u2_bn > u2_max) u2_bn = u2_max;
				if (u2_bn < u2_min) u2_bn = u2_min;

				if (u3_fn > u3_max) u3_fn = u3_max;
				if (u3_fn < u3_min) u3_fn = u3_min;

				if (u3_bn > u3_max) u3_bn = u3_max;
				if (u3_bn < u3_min) u3_bn = u3_min;

				// put invariants to buffers
				rBuf[j + 1] = r_fn;
				qBuf[j] = q_bn;
				tfBuf[j + 1] = t_fn;
				tbBuf[j] = t_bn;
				u2fBuf[j + 1] = u2_fn;
				u2bBuf[j] = u2_bn;
				u3fBuf[j + 1] = u3_fn;
				u3bBuf[j] = u3_bn;
			}

			// boundary conditions along the X1 axis
			// assignment of boundary invatiants and add them to the buffer arrays

			// periodicity conditions
			/*rBuf[1] = rBuf[n1];
			tfBuf[1] = tfBuf[n1];
			u2fBuf[1] = u2fBuf[n1];
			u3fBuf[1] = u3fBuf[n1];*/

			// periodicity conditions
			/*qBuf[n1] = qBuf[1];
			tbBuf[n1] = tbBuf[1];
			u2bBuf[n1] = u2bBuf[1];
			u3bBuf[n1] = u3bBuf[1];*/

			// no-slip conditions
			// j == 1
			ro0_b = ronCon[i][n1 - 1][k];
			ro0_f = ronCon[i][1][k];

			qn = qBuf[1];
			un = 0;
			pn = -qn*sound*ro0_g;
			ro_n = ro0_g + pn / (sound*sound);

			tn = t0;
			u2_n = 0;
			u3_n = 0;

			p1[i][1][k] = pn;
			u11[i][1][k] = un;
			ro1[i][1][k] = ro_n;
			t1[i][1][k] = tn;
			u21[i][1][k] = u2_n;
			u31[i][1][k] = u3_n;


			// j == n1
			rn = rBuf[n1];

			un = 0;
			pn = rn*sound*ro0_g;
			ro_n = ro0_g + pn / (sound*sound);

			tn = t0;
			u2_n = 0;
			u3_n = 0;

			p1[i][n1][k] = pn;
			u11[i][n1][k] = un;
			ro1[i][n1][k] = ro_n;
			t1[i][n1][k] = tn;
			u21[i][n1][k] = u2_n;
			u31[i][n1][k] = u3_n;

			// the flow variables calculations
			for (int j = 2; j < n1; j++)
			{
				ro0_b = ronCon[i][j - 1][k];
				ro0_f = ronCon[i][j][k];

				rn = rBuf[j];
				qn = qBuf[j];

				pn = (rn - qn)*sound*ro0_g / 2;
				un = (rn + qn) / 2;

				ro_n = ro0_g + pn / (sound*sound);

				ucf = u1nCon[i][j][k];
				ucb = u1nCon[i][j - 1][k];

				if (ucf >= 0 && ucb >= 0)
				{
					tn = tfBuf[j];
					u2_n = u2fBuf[j];
					u3_n = u3fBuf[j];
				}
				else if (ucf <= 0 && ucb <= 0)
				{
					tn = tbBuf[j];
					u2_n = u2bBuf[j];
					u3_n = u3bBuf[j];
				}
				else if (ucb >= 0 && ucf <= 0) 
				{
					if (ucb > -ucf) 
					{
						tn = tfBuf[j];
						u2_n = u2fBuf[j];
						u3_n = u3fBuf[j];
					}
					else 
					{
						tn = tbBuf[j];
						u2_n = u2bBuf[j];
						u3_n = u3bBuf[j];
					}
				}
				else if (ucb <= 0 && ucf >= 0)
				{
					tn = tnCon[i][j][k] + tnCon[i][j - 1][k] - t1[i][j][k];
					u2_n = u2nCon[i][j][k] + u2nCon[i][j - 1][k] - u21[i][j][k];
					u3_n = u3nCon[i][j][k] + u3nCon[i][j - 1][k] - u31[i][j][k];
				}

				p1[i][j][k] = pn;
				u11[i][j][k] = un;
				ro1[i][j][k] = ro_n;
				t1[i][j][k] = tn;
				u21[i][j][k] = u2_n;
				u31[i][j][k] = u3_n;
			}
		}
	}

	// flow variables calculation on DS2 faces orthogonal X2 axis 

	// bypassing along the X2 axis

	double xle, xlw, xls;

	double u1_fn, u1_bn;

	double gu1;

	double u1_max, u1_min, u1_n;

	for (int k = 1; k < n3; k++)
	{
		for (int j = 1; j < n1; j++)
		{
			xle = 1 + (l - 1)*(x1[j + 1] - 1);
			xlw = 1 + (l - 1)*(x1[j] - 1);
			xls = 0.5*(xle + xlw);
			
			for (int i = 1; i < n2; i++)
			{
				u2_f = u22[i + 1][j][k];
				u2_b = u22[i][j][k];
				u2_cn = u2nCon[i][j][k];
				u2_c = u2Con[i][j][k];

				u1_f = u12[i + 1][j][k];
				u1_b = u12[i][j][k];
				u1_cn = u1nCon[i][j][k];
				u1_c = u1Con[i][j][k];

				u3_f = u32[i + 1][j][k];
				u3_b = u32[i][j][k];
				u3_cn = u3nCon[i][j][k];
				u3_c = u3Con[i][j][k];

				ro_f = ro2[i + 1][j][k];
				ro_b = ro2[i][j][k];
				ro_cn = ronCon[i][j][k];
				ro_c = roCon[i][j][k];

				t_f = t2[i + 1][j][k];
				t_b = t2[i][j][k];
				t_cn = tnCon[i][j][k];
				t_c = tCon[i][j][k];

				p_f = p2[i + 1][j][k];
				p_b = p2[i][j][k];
				p_cn = sound*sound*(ro_cn - ro0_g);
				p_c = sound*sound*(ro_c - ro0_g);

				// invariant calculation
				r_f = u2_f + p_f / (ro0_g*sound);
				r_b = u2_b + p_b / (ro0_g*sound);
				r_cn = u2_cn + p_cn / (ro0_g*sound);
				r_c = u2_c + p_c / (ro0_g*sound);

				r_fn = 2 * r_cn - r_b;

				q_f = u2_f - p_f / (ro0_g*sound);
				q_b = u2_b - p_b / (ro0_g*sound);
				q_cn = u2_cn - p_cn / (ro0_g*sound);
				q_c = u2_c - p_c / (ro0_g*sound);

				q_bn = 2 * q_cn - q_f;

				t_fn = 2 * t_cn - t_b;
				t_bn = 2 * t_cn - t_f;

				u1_fn = 2 * u1_cn - u1_b;
				u1_bn = 2 * u1_cn - u1_f;

				u3_fn = 2 * u3_cn - u3_b;
				u3_bn = 2 * u3_cn - u3_f;

				// the permissible range of changes
				gr = 2 * (r_cn - r_c) / dt + (u2_cn + sound)*(r_f - r_b) / dx2 / xls;
				gq = 2 * (q_cn - q_c) / dt + (u2_cn - sound)*(q_f - q_b) / dx2 / xls;
				gt = 2 * (t_cn - t_c) / dt + u2_cn*(t_f - t_b) / dx2 / xls;
				gu1 = 2 * (u1_cn - u1_c) / dt + u2_cn*(u1_f - u1_b) / dx2 / xls;
				gu3 = 2 * (u3_cn - u3_c) / dt + u2_cn*(u3_f - u3_b) / dx2 / xls;


				// RMAX=MAX(RF,RC,RB) +dt*GR
				rmax = max3d(r_f, r_c, r_b) + dt*gr;

				// RMIN=MIN(RF,RC,RB) +dt*GR
				rmin = min3d(r_f, r_c, r_b) + dt*gr;

				// QMAX=MAX(QF,QC,QB) +dt*GQ
				qmax = max3d(q_f, q_c, q_b) + dt*gq;

				// QMIN=MIN(QF,QC,QB) +dt*GQ
				qmin = min3d(q_f, q_c, q_b) + dt*gq;

				// TMAX=MAX(TF,TC,TB) +dt*GT
				tmax = max3d(t_f, t_c, t_b) + dt*gt;

				// TMIN=MIN(TF,TC,TB) +dt*GT
				tmin = min3d(t_f, t_c, t_b) + dt*gt;

				// U1MAX=MAX(U1F,U1C,U1B) +dt*GU1
				u1_max = max3d(u1_f, u1_c, u1_b) + dt*gu1;

				// U1MIN=MIN(U1F,U1C,U1B) +dt*GU1
				u1_min = min3d(u1_f, u1_c, u1_b) + dt*gu1;

				// U3MAX=MAX(U3F,U3C,U3B) +dt*GU3
				u3_max = max3d(u3_f, u3_c, u3_b) + dt*gu3;

				// U3MIN=MIN(U3F,U3C,U3B) +dt*GU3 
				u3_min = min3d(u3_f, u3_c, u3_b) + dt*gu3;

				// invariants correction
				if (r_fn > rmax) r_fn = rmax;
				if (r_fn < rmin) r_fn = rmin;

				if (q_bn > qmax) q_bn = qmax;
				if (q_bn < qmin) q_bn = qmin;

				if (t_fn > tmax) t_fn = tmax;
				if (t_fn < tmin) t_fn = tmin;

				if (t_bn > tmax) t_bn = tmax;
				if (t_bn < tmin) t_bn = tmin;

				if (u1_fn > u1_max) u1_fn = u1_max;
				if (u1_fn < u1_min) u1_fn = u1_min;

				if (u1_bn > u1_max) u1_bn = u1_max;
				if (u1_bn < u1_min) u1_bn = u1_min;

				if (u3_fn > u3_max) u3_fn = u3_max;
				if (u3_fn < u3_min) u3_fn = u3_min;

				if (u3_bn > u3_max) u3_bn = u3_max;
				if (u3_bn < u3_min) u3_bn = u3_min;

				// put invariants to buffers
				// ==================================================
				// !!! IMPORTANT !!!
				// ==================================================
				// u2fBuf and u2bBuf are actially the u1fBuf and u1bBuf
				// It's not an error. We do it to save dynamic memory
				// ==================================================
				rBuf[i + 1] = r_fn;
				qBuf[i] = q_bn;
				tfBuf[i + 1] = t_fn;
				tbBuf[i] = t_bn;
				u2fBuf[i + 1] = u1_fn;
				u2bBuf[i] = u1_bn;
				u3fBuf[i + 1] = u3_fn;
				u3bBuf[i] = u3_bn;
			}

			// boundary conditions along the X2 axis
			// assignment of boundary invatiants and add them to the buffer arrays

			// periodicity conditions
			rBuf[1] = rBuf[n2];
			tfBuf[1] = tfBuf[n2];
			u2fBuf[1] = u2fBuf[n2];
			u3fBuf[1] = u3fBuf[n2];

			// periodicity conditions
			qBuf[n2] = qBuf[1];
			tbBuf[n2] = tbBuf[1];
			u2bBuf[n2] = u2bBuf[1];
			u3bBuf[n2] = u3bBuf[1];

			// the flow variables calculations
			for (int i = 1; i <= n2; i++)
			{
				ro0_b = ronCon[i - 1][j][k];
				ro0_f = ronCon[i][j][k];

				rn = rBuf[i];
				qn = qBuf[i];

				pn = (rn - qn)*sound*ro0_g / 2;
				un = (rn + qn) / 2;

				ro_n = ro0_g + pn / (sound*sound);

				ucf = u2nCon[i][j][k];
				ucb = u2nCon[i - 1][j][k];

				if (ucf >= 0 && ucb >= 0)
				{
					tn = tfBuf[i];
					u1_n = u2fBuf[i];
					u3_n = u3fBuf[i];
				}
				else if (ucf <= 0 && ucb <= 0)
				{
					tn = tbBuf[i];
					u1_n = u2bBuf[i];
					u3_n = u3bBuf[i];
				}
				else if (ucb >= 0 && ucf <= 0)
				{
					if (ucb > -ucf)
					{
						tn = tfBuf[i];
						u1_n = u2fBuf[i];
						u3_n = u3fBuf[i];
					}
					else
					{
						tn = tbBuf[i];
						u1_n = u2bBuf[i];
						u3_n = u3bBuf[i];
					}
				}
				else
				if (ucb <= 0 && ucf >= 0)
				{
					tn = tnCon[i][j][k] + tnCon[i - 1][j][k] - t2[i][j][k];
					u1_n = u1nCon[i][j][k] + u1nCon[i - 1][j][k] - u12[i][j][k];
					u3_n = u3nCon[i][j][k] + u3nCon[i - 1][j][k] - u32[i][j][k];
				}

				p2[i][j][k] = pn;
				u22[i][j][k] = un;
				ro2[i][j][k] = ro_n;
				t2[i][j][k] = tn;
				u12[i][j][k] = u1_n;
				u32[i][j][k] = u3_n;
			}
		}
	}

	// flow variables calculation on DS3 faces orthogonal X3 axis 

	// bypassing along the X3 axis

	for (int i = 1; i < n2; i++)
	{
		for (int j = 1; j < n1; j++)
		{
			for (int k = 1; k < n3; k++)
			{
				u3_f = u33[i][j][k + 1];
				u3_b = u33[i][j][k];
				u3_cn = u3nCon[i][j][k];
				u3_c = u3Con[i][j][k];

				u1_f = u13[i][j][k + 1];
				u1_b = u13[i][j][k];
				u1_cn = u1nCon[i][j][k];
				u1_c = u1Con[i][j][k];

				u2_f = u23[i][j][k + 1];
				u2_b = u23[i][j][k];
				u2_cn = u2nCon[i][j][k];
				u2_c = u2Con[i][j][k];

				ro_f = ro3[i][j][k + 1];
				ro_b = ro3[i][j][k];
				ro_cn = ronCon[i][j][k];
				ro_c = roCon[i][j][k];

				t_f = t3[i][j][k + 1];
				t_b = t3[i][j][k];
				t_cn = tnCon[i][j][k];
				t_c = tCon[i][j][k];

				p_f = p3[i][j][k + 1];
				p_b = p3[i][j][k];
				p_cn = sound*sound*(ro_cn - ro0_g);
				p_c = sound*sound*(ro_c - ro0_g);

				// invariant calculation
				r_f = u3_f + p_f / (ro0_g*sound);
				r_b = u3_b + p_b / (ro0_g*sound);
				r_cn = u3_cn + p_cn / (ro0_g*sound);
				r_c = u3_c + p_c / (ro0_g*sound);

				r_fn = 2 * r_cn - r_b;

				q_f = u3_f - p_f / (ro0_g*sound);
				q_b = u3_b - p_b / (ro0_g*sound);
				q_cn = u3_cn - p_cn / (ro0_g*sound);
				q_c = u3_c - p_c / (ro0_g*sound);

				q_bn = 2 * q_cn - q_f;

				t_fn = 2 * t_cn - t_b;
				t_bn = 2 * t_cn - t_f;

				u2_fn = 2 * u2_cn - u2_b;
				u2_bn = 2 * u2_cn - u2_f;

				u1_fn = 2 * u1_cn - u1_b;
				u1_bn = 2 * u1_cn - u1_f;

				// the permissible range of changes
				gr = 2 * (r_cn - r_c) / dt + (u3_cn + sound)*(r_f - r_b) / dx3;
				gq = 2 * (r_cn - r_c) / dt + (u3_cn - sound)*(q_f - q_b) / dx3;

				gt = 2 * (t_cn - t_c) / dt + u3_cn*(t_f - t_b) / dx3;
				gu1 = 2 * (u1_cn - u1_c) / dt + u3_cn*(u1_f - u1_b) / dx3;
				gu2 = 2 * (u2_cn - u2_c) / dt + u3_cn*(u2_f - u2_b) / dx3;

				// RMAX=MAX(RF,RC,RB) +dt*GR
				rmax = max3d(r_f, r_c, r_b) + dt*gr;

				// RMIN=MIN(RF,RC,RB) +dt*GR
				rmin = min3d(r_f, r_c, r_b) + dt*gr;

				// QMAX=MAX(QF,QC,QB) +dt*GQ
				qmax = max3d(q_f, q_c, q_b) + dt*gq;

				// QMIN=MIN(QF,QC,QB) +dt*GQ
				qmin = min3d(q_f, q_c, q_b) + dt*gq;

				// TMAX=MAX(TF,TC,TB) +dt*GT
				tmax = max3d(t_f, t_c, t_b) + dt*gt;

				// TMIN=MIN(TF,TC,TB) +dt*GT
				tmin = min3d(t_f, t_c, t_b) + dt*gt;

				// U1MAX=MAX(U1F,U1C,U1B) +dt*GU1
				u1_max = max3d(u1_f, u1_c, u1_b) + dt*gu1;

				// U1MIN=MIN(U1F,U1C,U1B) +dt*GU1
				u1_min = min3d(u1_f, u1_c, u1_b) + dt*gu1;

				// U2MAX=MAX(U2F,U2C,U2B) +dt*GU2
				u2_max = max3d(u2_f, u2_c, u2_b) + dt*gu2;

				// U2MIN=MIN(U2F,U2C,U2B) +dt*GU2
				u2_min = min3d(u2_f, u2_c, u2_b) + dt*gu2;

				// invariants correction
				if (r_fn > rmax) r_fn = rmax;
				if (r_fn < rmin) r_fn = rmin;

				if (q_bn > qmax) q_bn = qmax;
				if (q_bn < qmin) q_bn = qmin;

				if (t_fn > tmax) t_fn = tmax;
				if (t_fn < tmin) t_fn = tmin;

				if (t_bn > tmax) t_bn = tmax;
				if (t_bn < tmin) t_bn = tmin;

				if (u1_fn > u1_max) u1_fn = u1_max;
				if (u1_fn < u1_min) u1_fn = u1_min;

				if (u1_bn > u1_max) u1_bn = u1_max;
				if (u1_bn < u1_min) u1_bn = u1_min;

				if (u2_fn > u2_max) u2_fn = u2_max;
				if (u2_fn < u2_min) u2_fn = u2_min;

				if (u2_bn > u2_max) u2_bn = u2_max;
				if (u2_bn < u2_min) u2_bn = u2_min;

				// put invariants to buffers
				// ====================================================
				// !!! IMPORTANT !!!
				// ====================================================
				// u2fBuf and u2bBuf are actially the u1fBuf and u1bBuf
				// u3fBuf and u3bBuf are actially the u2fBuf and u2bBuf
				// It's not an error. We do it to save dynamic memory
				// ====================================================
				rBuf[k + 1] = r_fn;
				qBuf[k] = q_bn;
				tfBuf[k + 1] = t_fn;
				tbBuf[k] = t_bn;
				u2fBuf[k + 1] = u1_fn;
				u2bBuf[k] = u1_bn;
				u3fBuf[k + 1] = u2_fn;
				u3bBuf[k] = u2_bn;
			}

			// boundary conditions along the X3 axis
			// assignment of boundary invatiants and add them to the buffer arrays

			// periodicity conditions
			rBuf[1] = rBuf[n3];
			tfBuf[1] = tfBuf[n3];
			u2fBuf[1] = u2fBuf[n3];
			u3fBuf[1] = u3fBuf[n3];

			// periodicity conditions
			qBuf[n3] = qBuf[1];
			tbBuf[n3] = tbBuf[1];
			u2bBuf[n3] = u2bBuf[1];
			u3bBuf[n3] = u3bBuf[1];

			// the flow variables calculations
			for (int k = 1; k <= n3; k++)
			{
				ro0_b = ronCon[i][j][k - 1];
				ro0_f = ronCon[i][j][k];

				rn = rBuf[k];
				qn = qBuf[k];

				pn = (rn - qn)*sound*ro0_g / 2;
				un = (rn + qn) / 2;

				ro_n = ro0_g + pn / (sound*sound);

				ucf = u3nCon[i][j][k];
				ucb = u3nCon[i][j][k - 1];

				if (ucf >= 0 && ucb >= 0)
				{
					tn = tfBuf[j];
					u1_n = u2fBuf[j];
					u2_n = u3fBuf[j];
				}
				else if (ucf <= 0 && ucb <= 0)
				{
					tn = tbBuf[j];
					u1_n = u2bBuf[j];
					u2_n = u3bBuf[j];
				}
				else if (ucb >= 0 && ucf <= 0)
				{
					if (ucb > -ucf)
					{
						tn = tfBuf[j];
						u1_n = u2fBuf[j];
						u2_n = u3fBuf[j];
					}
					else
					{
						tn = tbBuf[j];
						u1_n = u2bBuf[j];
						u2_n = u3bBuf[j];
					}
				}
				else
				if (ucb <= 0 && ucf >= 0)
				{
					tn = tnCon[i][j][k] + tnCon[i][j][k - 1] - t3[i][j][k];
					u1_n = u1nCon[i][j][k] + u1nCon[i][j][k - 1] - u13[i][j][k];
					u2_n = u2nCon[i][j][k] + u2nCon[i][j][k - 1] - u23[i][j][k];
				}

				p3[i][j][k] = pn;
				u33[i][j][k] = un;
				ro3[i][j][k] = ro_n;
				t3[i][j][k] = tn;
				u13[i][j][k] = u1_n;
				u23[i][j][k] = u2_n;
			}

			// inlet conditions
			/*qn = qBuf[1];
			rn = u3Inlet + (ronCon[i][j][1] - ro0_g)*sound / ro0_g;
			un = (rn + qn) / 2;
			pn = (rn - qn)*sound*ro0_g / 2;
			ro_n = ro0_g + pn / sound / sound;
			u2_n = u2Inlet;
			u1_n = u1Inlet;
			tn = tInlet;
			p3[i][j][1] = pn;
			u33[i][j][1] = un;
			ro3[i][j][1] = ro_n;
			t3[i][j][1] = tn;
			u23[i][j][1] = u2_n;
			u13[i][j][1] = u1_n;*/

			// outlet conditions
			/*rn = rBuf[n3];
			pn = pOutlet;
			un = rn - pn / ro0_g / sound;
			tn = tfBuf[n3];
			u2_n = u2fBuf[n3];
			u1_n = u3fBuf[n3];
			p3[i][j][n3] = pn;
			u33[i][j][n3] = un;
			ro3[i][j][n3] = ro_n;
			t3[i][j][n3] = tn;
			u23[i][j][n3] = u2_n;
			u13[i][j][n3] = u1_n;*/
		}
	}
}

void WriteDataTecplot()
{
	char filename[50];

	//sprintf(filename, "out_%d.tec", nStep);
	sprintf_s(filename, 50, "out_%d.dat", nStep);

	//FILE *fd = fopen(filename, "w");
	FILE *fd;
	fopen_s(&fd, filename, "w");

	fprintf(fd, "TITLE=\"OUT\"\nVARIABLES=\"X\",\"Y\",\"Z\",\"U1\",\"U2\",\"U3\",\"PC\",\"TC\"\nZONE T=\"%10.4f\", I=%4d, J=%4d, K=%4d\n", TIME, n2 - 1, n1 - 1, n3 - 1);

	for (int i = 1; i < n2; ++i)
	{
		for (int j = 1; j < n1; ++j)
		{
			for (int k = 1; k < n3; ++k)
			{
				fprintf(fd, "%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n", x1[j], x2[i], x3[k], u1Con[i][j][k], u2Con[i][j][k], u3Con[i][j][k], ronCon[i][j][k], tnCon[i][j][k]);
			}
		}
	}

	fclose(fd);
}

void WriteDataParaView()
{
	char filename[50];

	//sprintf(filename, "out_%d.vtk", nStep);
	sprintf_s(filename, 50, "out_%d.vtk", nStep);

	//FILE *fd = fopen(filename, "w");
	FILE *fd;
	fopen_s(&fd, filename, "w");

	fprintf(fd, "# vtk DataFile Version 3.0\nvtk output\nASCII\n");
	fprintf(fd, "DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d", n1, n2, n3);

	fprintf(fd, "\nX_COORDINATES %d float\n", n1);
	for (int i = 1; i <= n1; i++)
	{
		fprintf(fd, "%f ", x1[i]);
	}

	fprintf(fd, "\nY_COORDINATES %d float\n", n2);
	for (int i = 1; i <= n2; i++)
	{
		fprintf(fd, "%f ", x2[i]);
	}

	fprintf(fd, "\nZ_COORDINATES %d float\n", n3);
	for (int i = 1; i <= n3; i++)
	{
		fprintf(fd, "%f ", x3[i]);
	}

	fprintf(fd, "\nCELL_DATA %d\nscalars U1 float\nLOOKUP_TABLE default\n", (n1-1)*(n2-1)*(n3-1));
	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n2; i++)
		{
			for (int j = 1; j < n1; j++)
			{
				fprintf(fd, "%f ", u1nCon[i][j][k]);
			}
		}
	}

	fprintf(fd, "\nscalars U2 float\nLOOKUP_TABLE default\n");
	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n2; i++)
		{
			for (int j = 1; j < n1; j++)
			{
				fprintf(fd, "%f ", u2nCon[i][j][k]);
			}
		}
	}

	fprintf(fd, "\nscalars U3 float\nLOOKUP_TABLE default\n");
	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n2; i++)
		{
			for (int j = 1; j < n1; j++)
			{
				fprintf(fd, "%f ", u3nCon[i][j][k]);
			}
		}
	}

	fprintf(fd, "\nscalars PC float\nLOOKUP_TABLE default\n");
	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n2; i++)
		{
			for (int j = 1; j < n1; j++)
			{
				fprintf(fd, "%f ", ronCon[i][j][k]);
			}
		}
	}

	fprintf(fd, "\nscalars TC float\nLOOKUP_TABLE default\n");
	for (int k = 1; k < n3; k++)
	{
		for (int i = 1; i < n2; i++)
		{
			for (int j = 1; j < n1; j++)
			{
				fprintf(fd, "%f ", tnCon[i][j][k]);
			}
		}
	}

	fclose(fd);
}

double ***allocate3D(int dim1, int dim2, int dim3)
{
	double ***arr;

	arr = calloc(dim1, sizeof(double**));

	for (int i = 0; i < dim1; ++i) {
		arr[i] = calloc(dim2, sizeof(double*));
		for (int j = 0; j < dim2; ++j) {
			arr[i][j] = calloc(dim3, sizeof(double));
		}
	}

	return arr;
}

void deallocate3D(double*** arr, int dim1, int dim2)
{
	for (int i = 0; i < dim1; ++i) {
		for (int j = 0; j < dim2; ++j) {
			free(arr[i][j]);
		}
		free(arr[i]);
	}
	free(arr);
	arr = NULL;
}

void FreeMemory()
{
	free(x1);
	free(x2);
	free(x3);

	deallocate3D(roCon, n2 + 1, n1 + 1);
	deallocate3D(u1Con, n2 + 1, n1 + 1);
	deallocate3D(u2Con, n2 + 1, n1 + 1);
	deallocate3D(u3Con, n2 + 1, n1 + 1);
	deallocate3D(tCon, n2 + 1, n1 + 1);

	deallocate3D(ronCon, n2 + 1, n1 + 1);
	deallocate3D(u1nCon, n2 + 1, n1 + 1);
	deallocate3D(u2nCon, n2 + 1, n1 + 1);
	deallocate3D(u3nCon, n2 + 1, n1 + 1);
	deallocate3D(tnCon, n2 + 1, n1 + 1);

	deallocate3D(ro1, n2 + 2, n1 + 2);
	deallocate3D(t1, n2 + 2, n1 + 2);
	deallocate3D(u11, n2 + 2, n1 + 2);
	deallocate3D(u21, n2 + 2, n1 + 2);
	deallocate3D(u31, n2 + 2, n1 + 2);
	deallocate3D(p1, n2 + 2, n1 + 2);

	deallocate3D(ro2, n2 + 2, n1 + 2);
	deallocate3D(t2, n2 + 2, n1 + 2);
	deallocate3D(u12, n2 + 2, n1 + 2);
	deallocate3D(u22, n2 + 2, n1 + 2);
	deallocate3D(u32, n2 + 2, n1 + 2);
	deallocate3D(p2, n2 + 2, n1 + 2);

	deallocate3D(ro3, n2 + 2, n1 + 2);
	deallocate3D(t3, n2 + 2, n1 + 2);
	deallocate3D(u13, n2 + 2, n1 + 2);
	deallocate3D(u23, n2 + 2, n1 + 2);
	deallocate3D(u33, n2 + 2, n1 + 2);
	deallocate3D(p3, n2 + 2, n1 + 2);

    deallocate3D(f1, n2 + 1, n1 + 1);
	deallocate3D(f2, n2 + 1, n1 + 1);
	deallocate3D(f3, n2 + 1, n1 + 1);
	
    deallocate3D(sigm11, n2 + 2, n1 + 2);
	deallocate3D(sigm21, n2 + 2, n1 + 2);
	deallocate3D(sigm31, n2 + 2, n1 + 2);

	deallocate3D(sigm12, n2 + 2, n1 + 2);
	deallocate3D(sigm22, n2 + 2, n1 + 2);
	deallocate3D(sigm32, n2 + 2, n1 + 2);

	deallocate3D(sigm13, n2 + 2, n1 + 2);
	deallocate3D(sigm23, n2 + 2, n1 + 2);
	deallocate3D(sigm33, n2 + 2, n1 + 2);

	free(rBuf);
	free(qBuf);
	free(tfBuf);
	free(tbBuf);
	free(u2fBuf);
	free(u2bBuf);
	free(u3fBuf);
	free(u3bBuf);
}

double min3d(double x1, double x2, double x3)
{
	x1 = x1 > x2 ? x2 : x1;
	x1 = x1 > x3 ? x3 : x1;
	return x1;
}

double max3d(double x1, double x2, double x3)
{
	x1 = x1 < x2 ? x2 : x1;
	x1 = x1 < x3 ? x3 : x1;
	return x1;
}
