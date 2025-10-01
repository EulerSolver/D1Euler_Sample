#ifndef WTP_CPP
#define WTP_CPP
#include "outputFile.h"
#include <fstream>
using namespace std;

outputFile::outputFile(PhysicsValueAndParameters phy)
{
	ofstream outfile;
	string NX = to_string(phy.Nx);
	string Eps = to_string(phy.eps);
	string OutputName = "./result/" + phy.whatflow + phy.whatscheme + NX + phy.ReconstructType + " - " + ".dat";

	outfile.open(OutputName, ios::out | ios::trunc);
	if (!outfile.is_open())
	{
		std::cerr << "Error: Path does not exist." << std::endl;
	}
	outfile << "Title=\"D1\"" << endl;
	outfile << "variables=\"x\",\"rho\",\"u\",\"p\",\"a\",\"M\"" << endl;
	outfile << "zone t=\"Box\",F=POINT" << endl;
	for (int i = phy.Ghost; i < phy.NxAug - phy.Ghost; i++)
	{
		outfile << phy.x[i] << "\t" << phy.Var.dens[i] << "\t" << phy.Var.u[i] << "\t" << phy.Var.p[i]
				<< "\t" << sqrt((phy.gamma) * phy.Var.p[i] / phy.Var.dens[i])
				<< "\t" << phy.Var.u[i] / sqrt((phy.gamma) * phy.Var.p[i] / phy.Var.dens[i]) << "\n";
	}
	outfile.close();
}


#endif