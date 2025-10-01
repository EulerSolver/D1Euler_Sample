#ifndef BOUNDARYCDS_CPP
#define BOUNDARYCDS_CPP
#include "boundaryCDs.h"

void zeroGrad(double *Quantity, int Ghost, int Nx)
{
	int NxAug = Nx + 2 * Ghost;
	if (Ghost == 3)
	{
		Quantity[0] = Quantity[3];
		Quantity[1] = Quantity[3];
		Quantity[2] = Quantity[3];

		Quantity[NxAug - 1] = Quantity[NxAug - 4];
		Quantity[NxAug - 2] = Quantity[NxAug - 4];
		Quantity[NxAug - 3] = Quantity[NxAug - 4];
	}

}

void reflect(double *Quantity, int Ghost, int Nx)
{
	int NxAug = Nx + 2 * Ghost;
	if (Ghost == 3)
	{

		Quantity[0] = -Quantity[5];
		Quantity[1] = -Quantity[4];
		Quantity[2] = -Quantity[3];

		Quantity[NxAug - 1] = -Quantity[NxAug - 6];
		Quantity[NxAug - 2] = -Quantity[NxAug - 5];
		Quantity[NxAug - 3] = -Quantity[NxAug - 4];

		// Quantity[0] = -Quantity[6];
		// Quantity[1] = -Quantity[5];
		// Quantity[2] = -Quantity[4];

		// Quantity[NxAug - 1] = -Quantity[NxAug - 7];
		// Quantity[NxAug - 2] = -Quantity[NxAug - 6];
		// Quantity[NxAug - 3] = -Quantity[NxAug - 5];
	}
	else if (4 == Ghost)
	{
		// std::cout<<"yes ,the ref boundary is called"<<std::endl;
		Quantity[0] = -Quantity[7];
		Quantity[1] = -Quantity[6];
		Quantity[2] = -Quantity[5];
		Quantity[3] = -Quantity[4];

		Quantity[NxAug - 1] = -Quantity[NxAug - 8];
		Quantity[NxAug - 2] = -Quantity[NxAug - 7];
		Quantity[NxAug - 3] = -Quantity[NxAug - 6];
		Quantity[NxAug - 4] = -Quantity[NxAug - 5];

		// Quantity[0] = -Quantity[8];
		// Quantity[1] = -Quantity[7];
		// Quantity[2] = -Quantity[6];
		// Quantity[3] = -Quantity[5];

		// Quantity[NxAug - 1] = -Quantity[NxAug - 9];
		// Quantity[NxAug - 2] = -Quantity[NxAug - 8];
		// Quantity[NxAug - 3] = -Quantity[NxAug - 7];
		// Quantity[NxAug - 4] = -Quantity[NxAug - 6];
	}
}

void symmetry(double *Quantity, int Ghost, int Nx)
{
	int NxAug = Nx + 2 * Ghost;
	if (Ghost == 3)
	{
		Quantity[0] = Quantity[5];
		Quantity[1] = Quantity[4];
		Quantity[2] = Quantity[3];

		Quantity[NxAug - 1] = Quantity[NxAug - 6];
		Quantity[NxAug - 2] = Quantity[NxAug - 5];
		Quantity[NxAug - 3] = Quantity[NxAug - 4];

		// Quantity[0] = Quantity[6];
		// Quantity[1] = Quantity[5];
		// Quantity[2] = Quantity[4];

		// Quantity[NxAug - 1] = Quantity[NxAug - 7];
		// Quantity[NxAug - 2] = Quantity[NxAug - 6];
		// Quantity[NxAug - 3] = Quantity[NxAug - 5];
	}
	else if (4 == Ghost)
	{
		Quantity[0] = Quantity[7];
		Quantity[1] = Quantity[6];
		Quantity[2] = Quantity[5];
		Quantity[3] = Quantity[4];

		Quantity[NxAug - 1] = Quantity[NxAug - 8];
		Quantity[NxAug - 2] = Quantity[NxAug - 7];
		Quantity[NxAug - 3] = Quantity[NxAug - 6];
		Quantity[NxAug - 4] = Quantity[NxAug - 5];

		// Quantity[0] = Quantity[8];
		// Quantity[1] = Quantity[7];
		// Quantity[2] = Quantity[6];
		// Quantity[3] = Quantity[5];

		// Quantity[NxAug - 1] = Quantity[NxAug - 9];
		// Quantity[NxAug - 2] = Quantity[NxAug - 8];
		// Quantity[NxAug - 3] = Quantity[NxAug - 7];
		// Quantity[NxAug - 4] = Quantity[NxAug - 6];
	}
	
}

#endif // !BOUNDARY_CPP

