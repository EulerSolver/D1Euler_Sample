#ifndef INI_CPP
#define INI_CPP

#include "InitialCDs.h"
// #include "PhysicsValueAndControl.h"
#include <iostream>
// #pragma warning(disable:26495)
#include <cmath>
#include <iostream>

InitailCDs::InitailCDs(PhysicsValueAndParameters &Phy)
{
	double *rho = Phy.Var.dens;
	double *u = Phy.Var.u;
	double *p = Phy.Var.p;
	double *rhoU = Phy.Var.momt;
	double *eng = Phy.Var.eng;

	if (Phy.whatflow == "Sod")
	{
		Phy.boundaryCD = "ZeroGrad";
		Phy.tf = 0.2;
		Phy.dx = (1.0 - (0.0)) / (Phy.Nx - 1);
		Phy.gamma = 1.4;
		for (int i = 0; i < Phy.NxAug; i++)
		{
			Phy.x[i] = 0 + (i - Phy.Ghost) * Phy.dx;
			if (Phy.x[i] < 0.5)
			{
				rho[i] = 1.0;
				u[i] = 0.0;
				p[i] = 1.0;
				
			}
			else
			{
				rho[i] = 0.125;
				u[i] = 0.0;
				p[i] = 0.1;
			}
			rhoU[i] = rho[i] * u[i];
			eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
		}

		zeroGrad(rho, Phy.Ghost, Phy.Nx);
		zeroGrad(u, Phy.Ghost, Phy.Nx);
		zeroGrad(p, Phy.Ghost, Phy.Nx);
		zeroGrad(rhoU, Phy.Ghost, Phy.Nx);
		zeroGrad(eng, Phy.Ghost, Phy.Nx);
	}
    else if (Phy.whatflow == "Blast")
	{
		Phy.boundaryCD = "reflect";
		Phy.tf = 0.038;
		Phy.dx = (1.0 - (0.0)) / (double)(Phy.Nx);
		Phy.gamma = 1.4;
		for (int i = 0; i < Phy.NxAug; i++)
		{
			Phy.x[i] = Phy.dx * (i + 0.5 - Phy.Ghost);
			if (Phy.x[i] <= 0.1)
			{
				rho[i] = 1.0;
				u[i] = 0.0;
				p[i] = 1000.0;
			}
			else if (Phy.x[i] >= 0.1 && Phy.x[i] <= 0.9)
			{
				rho[i] = 1.0;
				u[i] = 0.0;
				p[i] = 0.01;
			
			}
			else
			{
				rho[i] = 1.0;
				u[i] = 0.0;
				p[i] = 100;
			}
			rhoU[i] = rho[i] * u[i];
			eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
			symmetry(rho, Phy.Ghost, Phy.Nx);
			reflect(u, Phy.Ghost, Phy.Nx);
			symmetry(p, Phy.Ghost, Phy.Nx);
			reflect(rhoU, Phy.Ghost, Phy.Nx);
			symmetry(eng, Phy.Ghost, Phy.Nx);
		}
	}

	else if (Phy.whatflow == "ShuOsher")
	{
		Phy.boundaryCD = "ZeroGrad";
		Phy.tf = 1.8;
		Phy.dx = (5.0 - (-5.0)) / (Phy.Nx - 1);
		Phy.gamma = 1.4;
		double k = 5.0;
		for (int i = 0; i < Phy.NxAug; i++)
		{
			Phy.x[i] = -5 + (-Phy.Ghost + i) * Phy.dx;
			if (Phy.x[i] < -4.0)
			{
				rho[i] = 3.857143;
				u[i] = 2.629369;
				p[i] = 31.0 / 3.0;
				rhoU[i] = rho[i] * u[i];
				eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
			}
			else
			{
				rho[i] = 1.0 + 0.2 * sin(k * Phy.x[i]);
				u[i] = 0.0;
				p[i] = 1.0;
				rhoU[i] = rho[i] * u[i];
				eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
			}
			// cout << pv->dens[i]<<" ";
		}
		zeroGrad(rho, Phy.Ghost, Phy.Nx);
		zeroGrad(u, Phy.Ghost, Phy.Nx);
		zeroGrad(p, Phy.Ghost, Phy.Nx);
		zeroGrad(rhoU, Phy.Ghost, Phy.Nx);
		zeroGrad(eng, Phy.Ghost, Phy.Nx);
	}
	else if (Phy.whatflow == "Lax")
	{
		{
			Phy.boundaryCD = "ZeroGrad";
			Phy.tf = 0.26;
			Phy.dx = (1.0 - (-1.0)) / (Phy.Nx - 1.0);
			Phy.gamma = 1.4;

			for (int i = 0; i < Phy.NxAug; i++)
			{
				Phy.x[i] = -1.0 + (i - Phy.Ghost) * Phy.dx;
				if (Phy.x[i] <= 0.0)
				{
					rho[i] = 0.445;
					u[i] = 0.698;
					p[i] = 3.528;
					rhoU[i] = rho[i] * u[i];
					eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
				}
				else
				{
					rho[i] = 0.5;
					u[i] = 0.0;
					p[i] = 0.571;
					// Z[i] = 0.0;
					rhoU[i] = rho[i] * u[i];
					eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
					// RhoZ[i] = dens[i] * Z[i];
				}
			}

			zeroGrad(rho, Phy.Ghost, Phy.Nx);
			zeroGrad(u, Phy.Ghost, Phy.Nx);
			zeroGrad(p, Phy.Ghost, Phy.Nx);
			zeroGrad(rhoU, Phy.Ghost, Phy.Nx);
			zeroGrad(eng, Phy.Ghost, Phy.Nx);
		}
	}

	else if (Phy.whatflow == "LeBlanc")
	{
		double domainL = 0.0, domainMid = 10, domainR = 20;
		double rhoL = 2, uL = 0.0, pL = 1e9;
		double rhoR = 1e-3, uR = 0.0, pR = 1.0;
		double tt = 1e-4;
		double gamma = 1.4;

		// Phy.boundary = "ref";
		Phy.boundaryCD = "fix";
		Phy.tf = tt;
		Phy.dx = (domainR - domainL) / (Phy.Nx - 1);
		Phy.gamma = gamma;

		for (int i = 0; i < Phy.NxAug; i++)
		{
			// Phy.x[i] =  (i + 0.5 - Phy.Ghost) * Phy.dx;
			Phy.x[i] = domainL + (i - Phy.Ghost) * Phy.dx;

			if (Phy.x[i] <= domainMid)
			{
				rho[i] = rhoL;
				u[i] = uL;
				p[i] = pL;
				rhoU[i] = rho[i] * u[i];
				eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
			}
			else
			{
				rho[i] = rhoR;
				u[i] = uR;
				p[i] = pR;
				rhoU[i] = rho[i] * u[i];
				eng[i] = 0.5 * rho[i] * u[i] * u[i] + p[i] / (Phy.gamma - 1);
			}
			zeroGrad(rho, Phy.Ghost, Phy.Nx);
			zeroGrad(u, Phy.Ghost, Phy.Nx);
			zeroGrad(p, Phy.Ghost, Phy.Nx);
			// period(Z, Phy.Ghost, Phy.Nx);
			zeroGrad(rhoU, Phy.Ghost, Phy.Nx);
			zeroGrad(eng, Phy.Ghost, Phy.Nx);
		}
	}

	else
	{
		std::cout << "No such flow-name. Please Check the initial file (in.in)!" << std::endl;
	}
}
#endif
