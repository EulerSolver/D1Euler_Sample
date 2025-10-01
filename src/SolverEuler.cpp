
#ifndef D1_CPP
#define D1_CPP
#include "D1EulerSolver.h"


D1Euler::D1Euler(std::string whatflow, std::string ReconstructType,
				 double tf, double CFL, int nx,
				 double eps)
{
	// std::cout<<"Im Here"<<std::endl;
	Phy.Nx = nx;
	Phy.Ghost = 3;
	Phy.NxAug = nx + Phy.Ghost * 2;
	//-------------
	allocateAll(Phy); // Calloc physical quantities according to the grid size;
	//-------------
	Phy.whatflow = whatflow;
	Phy.ReconstructType = ReconstructType;
	Phy.CFL = CFL;
	Phy.Nx = nx;
	Phy.eps = eps;
	Phy.VarWorking = Phy.Var;
}




void D1Euler::boundaryConditionsCall()
{

	{

		if (Phy.boundaryCD == "zeroGrad")
		{
			zeroGrad(Phy.VarWorking.dens, Phy.Ghost, Phy.Nx);
			zeroGrad(Phy.VarWorking.momt, Phy.Ghost, Phy.Nx);
			zeroGrad(Phy.VarWorking.eng, Phy.Ghost, Phy.Nx);
		}
		else if (Phy.boundaryCD == "reflect")
		{
			symmetry(Phy.VarWorking.dens, Phy.Ghost, Phy.Nx);
			reflect(Phy.VarWorking.momt, Phy.Ghost, Phy.Nx);
			symmetry(Phy.VarWorking.eng, Phy.Ghost, Phy.Nx);
		}
	}
}

void D1Euler::Dt_Calculate()
{

	//--------------------------
	double *Rho = Phy.Var.dens;
	double *u = Phy.Var.u;
	double *p = Phy.Var.p;
	double *a = Phy.Var.a;
	double *E = Phy.Var.eng;
	double *RhoU = Phy.Var.momt;
	//-------------------------

	for (int i = Phy.Ghost; i < Phy.NxAug - Phy.Ghost; i++)
	{
		u[i] = RhoU[i] / Rho[i];
		p[i] = (E[i] - 0.5 * Rho[i] * u[i] * u[i]) * (Phy.gamma - 1);
		a[i] = sqrt(Phy.gamma * p[i] / Rho[i]);
		Phy.infs[i] = abs(u[i]) + a[i];
	}
	Phy.maxinfs = (*max_element(Phy.infs + Phy.Ghost, Phy.infs + Phy.NxAug - Phy.Ghost));
	Phy.dt = Phy.CFL * Phy.dx / Phy.maxinfs;
}

void D1Euler::get_u_p_a()
{

	double *Rho = Phy.VarWorking.dens;
	double *u = Phy.VarWorking.u;
	double *p = Phy.VarWorking.p;
	double *a = Phy.VarWorking.a;

	double *E = Phy.VarWorking.eng;
	double *RhoU = Phy.VarWorking.momt;

	for (int i = 0; i < Phy.NxAug; i++)
	{
		u[i] = RhoU[i] / Rho[i];
		p[i] = (E[i] - 0.5 * Rho[i] * u[i] * u[i]) * (Phy.gamma - 1);
		a[i] = sqrt(Phy.gamma * p[i] / Rho[i]);
	}
}

D1Euler::~D1Euler()
{
	freeAll(Phy);
}

#endif
