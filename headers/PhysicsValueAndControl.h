#pragma once
#ifndef PHYSICSVALUE_H
#define PHYSICSVALUE_H

#include <iostream>
#include "time.h"
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <stack>
#include <unordered_map>
#include <chrono>

std::istream &ReadIOFile_2(std::istream &ios, std::map<std::string, std::string> &Model);
void readParaFrom_in_in_file();

struct VarsUsed_Euler
{
	double *dens = nullptr;
	double *u = nullptr;
	double *p = nullptr;
	double *momt = nullptr;
	double *eng = nullptr;

	double *InternalE = nullptr;

	double *a = nullptr;
	double *H = nullptr;
	double *Entropy = nullptr;
	double *SqrtRho = nullptr;
};

struct VarOfRoeAv
{
	double *AvRho = nullptr;
	double *AvU = nullptr;
	double *AvP = nullptr;
	double *AvH = nullptr;
	double *AvA = nullptr;
};

struct Inlet
{
	double rho = 0.0;
	double u = 0.0;
	double p = 0.0;
	double mom = 0.0;
	double eng = 0.0;
};

void allocateVarsMemory(VarsUsed_Euler &Var, int size);

void allocateVarsOfRoeAV(VarOfRoeAv &AV, int size);

void FreeVarsAllocations(VarsUsed_Euler &Var, int size);

void freeVarsOfRoeAV(VarOfRoeAv &AV, int size);

struct PhysicsValueAndParameters
{
	VarsUsed_Euler VarWorking;
	VarsUsed_Euler Var;
	VarsUsed_Euler VarRK1;
	VarsUsed_Euler VarRK2;
	VarsUsed_Euler RHS;

	VarsUsed_Euler Flux_interface;
	VarsUsed_Euler FP;
	VarsUsed_Euler FM;

	VarOfRoeAv RoeAv;

	Inlet Inlet;

	std::string whatscheme;
	std::string ReconstructType;
	std::string whatflow;
	std::string boundaryCD;

	double *L = nullptr;
	double *R = nullptr;

	int Nx = 0.0;
	int Ghost = 3;
	int NxAug = 0.0;
	double CFL = 0.0;
	double gamma = 0;

	double eps = 0;

	double dx = 0;
	double dt = 0;
	double tf = 0;
	double maxinfs;

	double t = 0;

	double *x = nullptr;
	double *infs = nullptr;
};

void allocateAll(PhysicsValueAndParameters &phy);
void freeAll(PhysicsValueAndParameters &phy);

// void SetPara();

//-------------------

#endif
