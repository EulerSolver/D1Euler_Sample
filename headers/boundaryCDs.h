
#ifndef BOUNDARYCDS_H
#define BOUNDARYCDS_H


void zeroGrad(double* Quantity, int Ghost, int Nx);
void reflect(double* Quantity, int Ghost, int Nx);
void symmetry(double* Quantity, int Ghost, int Nx);
#endif