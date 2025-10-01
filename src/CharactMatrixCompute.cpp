
#ifndef D1_Charact_CPP
#define D1_Charact_CPP
#include "D1EulerSolver.h"

void D1Euler::CalculateCharactMatrices(int I, char flag)
{

	//--------------------------
	// Simple denotations of WorkingVaribales For reconstruction----Begin
	double *Rho = Phy.VarWorking.dens;
	double *u = Phy.VarWorking.u;
	double *p = Phy.VarWorking.p;

	double *E = Phy.VarWorking.eng;
	double *RhoU = Phy.VarWorking.momt;
	double *SqrtRho = Phy.VarWorking.SqrtRho;
	double *a = Phy.VarWorking.a;
	double *H = Phy.VarWorking.H;
	// Simple denotations of WorkingVaribales For reconstruction----End
	//-------------------------

	// Roe averaged variables--Begin
	double *AvRho = Phy.RoeAv.AvRho;
	double *AvH = Phy.RoeAv.AvH;
	double *AvU = Phy.RoeAv.AvU;
	double *AvA = Phy.RoeAv.AvA;
	double *AvP = Phy.RoeAv.AvP;
	// Roe averaged variables--End

	// Set the exception index ----the left-most---point----Begin
	if (Phy.Ghost - 1 == I)
	{
		SqrtRho[I] = sqrt((Rho[I]));
		H[I] = (E[I] + p[I]) / Rho[I];
		a[I] = sqrt(Phy.gamma * p[I] / Rho[I]);
	}
	// Set the exception index ----the left-most---point----End

	int iP1 = I + 1;

	{ // Roe average
		//-------------------------------------------------------
		//---------------Calculate Roe Averaged variables--Begin
		SqrtRho[iP1] = sqrt(Rho[iP1]);
		H[iP1] = (E[iP1] + p[iP1]) / Rho[iP1];

		double SqrtRhoLPlusSqrtRhoR = SqrtRho[I] + SqrtRho[iP1];
		double D1 = SqrtRho[I] / SqrtRhoLPlusSqrtRhoR;
		double D2 = SqrtRho[iP1] / SqrtRhoLPlusSqrtRhoR;

		AvRho[I] = (SqrtRho[I] + SqrtRho[iP1]) * (SqrtRho[I] + SqrtRho[iP1]) * 0.25;
		// AvRho[IChP] = (0.5 * SqrtRho[IChP] + 0.5 * SqrtRho[iP1]) * (0.5 * SqrtRho[IChP] + 0.5 * SqrtRho[iP1]);
		AvH[I] = D1 * H[I] + D2 * H[iP1];
		AvU[I] = D1 * u[I] + D2 * u[iP1];
		double uuAV = Phy.RoeAv.AvU[I] * Phy.RoeAv.AvU[I];
		AvA[I] = sqrt((Phy.RoeAv.AvH[I] - 0.5 * uuAV) * (Phy.gamma - 1.0));
		//---------------Calculate Roe Averaged variables--End
		//-------------------------------------------------------
	}

	double OneOverA = 1.0 / AvA[I];

	double XP = 0.5 * (Phy.gamma - 1) * OneOverA * OneOverA;
	double XPTimesUU = XP * AvU[I] * AvU[I];
	double XPTImesU = XP * AvU[I];
	double UOverA = AvU[I] / AvA[I];
	double ATimesU = AvA[I] * AvU[I];

	//---Jacobian matrix---Begin
	if ('n' == flag || 'f' == flag)
	{
		Phy.L[1] = (0.5 * XPTimesUU + 0.5 * UOverA);
		Phy.L[2] = ((-XPTImesU - 0.5 * OneOverA));
		Phy.L[3] = (XP);
		// Phy.L[4] = (1.0 - XPTimesUU);
		// Phy.L[5] = (2 * XPTImesU);
		// Phy.L[6] = (-2 * XP);
		Phy.L[7] = (0.5 * XPTimesUU - 0.5 * UOverA);
		Phy.L[8] = (-XPTImesU + 0.5 * OneOverA);
		Phy.L[9] = (XP);

		Phy.R[1] = 1.0;
		// Phy.R[2] = 1.0;
		Phy.R[3] = 1.0;
		Phy.R[4] = (AvU[I] - AvA[I]);
		// Phy.R[5] = AvU[I];
		Phy.R[6] = (AvU[I] + AvA[I]);
		Phy.R[7] = (AvH[I] - ATimesU);
		// Phy.R[8] = 0.5 * AvU[I] * AvU[I];
		Phy.R[9] = (AvH[I] + ATimesU);
	}

	if ('l' == flag || 'f' == flag)
	{
		Phy.L[4] = (1.0 - XPTimesUU);
		Phy.L[5] = (2 * XPTImesU);
		Phy.L[6] = (-2 * XP);
	}

	if ('f' == flag)
	{
		Phy.R[2] = 1.0;
		Phy.R[5] = AvU[I];
		Phy.R[8] = 0.5 * AvU[I] * AvU[I];
	}

	//---Jacobian matrix---End
}

#endif