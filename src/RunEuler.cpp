#ifndef D1_RUN_CPP
#define D1_RUN_CPP
#include "D1EulerSolver.h"

void D1Euler::Run()
{
	Phy.t = 0.0;
	int count = 0;
	while (Phy.t < Phy.tf)
	{
		Dt_Calculate();

		if (abs(Phy.t - Phy.tf) < Phy.dt)
			Phy.dt = Phy.tf - Phy.t;

		// if (count++ % 200 == 0)
		// 	std::cout << "Dt is " << Phy.dt << std::endl;

		TVD_RK3(1);
		Phy.t += Phy.dt;
	}
	get_u_p_a();
}

#endif