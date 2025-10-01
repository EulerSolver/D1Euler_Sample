#ifndef D1_RiemannSolver_CPP
#define D1_RiemannSolver_CPP
#include "D1EulerSolver.h"

void D1Euler::GLF_FVS()
{
	double *Rho = Phy.VarWorking.dens;
	double *u = Phy.VarWorking.u;
	double *p = Phy.VarWorking.p;
	double *a = Phy.VarWorking.a;

	double *E = Phy.VarWorking.eng;
	double *RhoU = Phy.VarWorking.momt;
	

	for (int i = Phy.Ghost; i < Phy.NxAug - Phy.Ghost; i++)
	{
		Phy.infs[i] = max(abs(u[i]), abs(u[i] + a[i]));
		Phy.infs[i] = max(Phy.infs[i], abs(u[i] - a[i]));
	}

	double MaxAlpha = *max_element(Phy.infs + Phy.Ghost, Phy.infs + Phy.NxAug - Phy.Ghost);

	for (int i = 0; i < Phy.NxAug; i++)
	{
		double f1 = RhoU[i];
		double f2 = RhoU[i] * u[i] + p[i];
		double f3 = u[i] * (E[i] + p[i]);

		Phy.FP.dens[i] = 0.5 * (f1 + MaxAlpha * Rho[i]);
		Phy.FP.momt[i] = 0.5 * (f2 + MaxAlpha * RhoU[i]);
		Phy.FP.eng[i] = 0.5 * (f3 + MaxAlpha * E[i]);

		// Phy.FM.dens[i] = 0.5 * (f1 - MaxAlpha * Rho[i]);
		// Phy.FM.momt[i] = 0.5 * (f2 - MaxAlpha * RhoU[i]);
		// Phy.FM.eng[i] = 0.5 * (f3 - MaxAlpha * E[i]);

		Phy.FM.dens[i] = (f1 - Phy.FP.dens[i]);
		Phy.FM.momt[i] = (f2 - Phy.FP.momt[i]);
		Phy.FM.eng[i] = (f3 - Phy.FP.eng[i]);
	}
}

#endif