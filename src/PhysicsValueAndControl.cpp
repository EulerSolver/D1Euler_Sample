#ifndef PHYSICSVALUE_CPP
#define PHYSICSVALUE_CPP
#include "PhysicsValueAndControl.h"
#include <unordered_map>

std::map<std::string, std::string> Para = {
    {"nx", ""},
    {"epsilon", ""},
    {"cfl", ""},
    {"whatflow", ""},
    {"ReconstructType", ""},
    {"Time", ""}
};

	// void SetPara()
	// {
	// Para["nx"] = " ";
	// Para["epsilon"] = " ";
	// Para["cfl"] = " ";
	// Para["whatflow"] = " ";
	// Para["ReconstructType"] = " ";
	// Para["Time"] = " ";
	// }

std::istream& ReadIOFile_2(std::istream& ios, std::map<std::string, std::string>& Model) {

	std::string word0;
	std::string word1;
	ios >> word0 >> word1;
	if (Model.find(word0) != Model.end()) {
		Model[word0] = word1;
	}
	else if (word0.size()) {
		std::cout << "make sure MAP discipine" << std::endl;
	}
	return ios;
}

void readParaFrom_in_in_file() {
    std::ifstream ff("in.in");
    auto ParaSize = Para.size();
    // std::cout<< " ParaSize is " << ParaSize<<std::endl;
    if (ff)
    {
        while (ReadIOFile_2(ff, Para))
            ParaSize--;
        if (ParaSize != 0)
        {
            std::cout << "Went Wrong, Check Matchs of ParaMap above and in.in file" << std::endl;
        }
    }
    else
    {
        std::cout << "in.in file don't exist, create one" << std::endl;
    }
    ff.close();
    for (const auto &pair : Para)
    {
        std::cout << pair.first << "= " << pair.second << std::endl;
    }
}


void allocateVarsMemory(VarsUsed_Euler &Var, int size) {
	Var.dens = (double*)calloc(size, sizeof(double));
	Var.u = (double*)calloc(size, sizeof(double));
	Var.p = (double*)calloc(size, sizeof(double));
	Var.momt = (double*)calloc(size, sizeof(double));
	Var.eng = (double*)calloc(size, sizeof(double));
	//Var.engG = (double*)calloc(size, sizeof(double));  //this special variable is for Gamma-independent characteristic-space;

	Var.InternalE = (double*)calloc(size, sizeof(double));

	Var.a = (double*)calloc(size, sizeof(double));
	Var.H = (double*)calloc(size, sizeof(double));
	Var.Entropy = (double*)calloc(size, sizeof(double));
	Var.SqrtRho = (double*)calloc(size, sizeof(double));
	
}

void FreeVarsAllocations(VarsUsed_Euler &Var, int size) {
	
	free(Var.dens );
	free(Var.u );
	free(Var.p );
	free(Var.momt );
	free(Var.eng );
	free(Var.InternalE );
	free(Var.a );
	free(Var.H );
	free(Var.Entropy);
	free(Var.SqrtRho);
	
}

void allocateVarsOfRoeAV(VarOfRoeAv& AV, int size) {
	AV.AvRho = (double*)calloc(size, sizeof(double));
	AV.AvU = (double*)calloc(size, sizeof(double));
	AV.AvP = (double*)calloc(size, sizeof(double));
	AV.AvH = (double*)calloc(size, sizeof(double));
	AV.AvA = (double*)calloc(size, sizeof(double));
	
}

void freeVarsOfRoeAV(VarOfRoeAv& AV, int size) {
	free(AV.AvRho );
	free(AV.AvU );
	free(AV.AvP);
	free(AV.AvH );
	free(AV.AvA );
	
}

void allocateAll(PhysicsValueAndParameters &Phy)
{
	// iniVarMultiCompEq4(Phy.VarWorking, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.Var, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.VarRK1, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.VarRK2, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.RHS, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.Flux_interface, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.FP, Phy.Nx + Phy.Ghost * 2);
	allocateVarsMemory(Phy.FM, Phy.Nx + Phy.Ghost * 2);
	
	allocateVarsOfRoeAV(Phy.RoeAv, Phy.NxAug);			   

	Phy.x = (double *)calloc(Phy.Nx + Phy.Ghost * 2, sizeof(double));
	Phy.infs = (double *)calloc(Phy.Nx + Phy.Ghost * 2, sizeof(double));

	//-----For the characteristic structure
	Phy.L = (double *)calloc(100, sizeof(double));
	Phy.R = (double *)calloc(100, sizeof(double));
	//-----For the characteristic structure
}

void freeAll(PhysicsValueAndParameters &Phy)
{
	
	FreeVarsAllocations(Phy.Var, Phy.Nx + Phy.Ghost * 2);
	FreeVarsAllocations(Phy.VarRK1, Phy.Nx + Phy.Ghost * 2);
	FreeVarsAllocations(Phy.VarRK2, Phy.Nx + Phy.Ghost * 2);
	FreeVarsAllocations(Phy.RHS, Phy.Nx + Phy.Ghost * 2);
	FreeVarsAllocations(Phy.Flux_interface, Phy.Nx + Phy.Ghost * 2);
	FreeVarsAllocations(Phy.FP, Phy.Nx + Phy.Ghost * 2);
	FreeVarsAllocations(Phy.FM, Phy.Nx + Phy.Ghost * 2);
	
	freeVarsOfRoeAV(Phy.RoeAv, Phy.NxAug);			   

	free(Phy.x);
	free(Phy.L);
	free(Phy.R);
	
}

#endif // !PHYSICSVALUE_CPP
