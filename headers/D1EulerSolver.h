#ifndef D1Euler_H
#define D1Euler_H

#include "PhysicsValueAndControl.h"
#include "boundaryCDs.h"
using namespace std;
class InitailCDs;
class D1Euler
{
public:
	friend class InitailCDs;
	void Run();
	PhysicsValueAndParameters Phy;

	D1Euler(std::string whatflow, std::string ReconstructType,
			double tf, double CFL, int nx,
			double eps);
	~D1Euler();
private:
	

	void (D1Euler::*spacialScheme)(double *const Quantiy, int IndexFromSorce, int sign, double &RecontReceipt, double[], double[]);

	void boundaryConditionsCall();
	void Dt_Calculate();
	void get_u_p_a();
	void CalculateCharactMatrices(int iChP, char flag);

	void TVD_RK3(int order);

	void ChooseReconstructionMethod_AndGo();

	void Component_Wise();
	void Characteristic_Wise();
	void ImprovedCo();
	void PressureEntropyReplacement();

	void weno5(double *__restrict__ Quantiy, int IndexFromSorce, int sign, double &RecontReceipt, double* __restrict__ Q, double*__restrict__ W, const char flag);

	void GLF_FVS();
};
#endif
