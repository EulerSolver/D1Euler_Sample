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

	}
}

#endif // !BOUNDARY_CPP
