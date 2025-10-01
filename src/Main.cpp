#ifndef RunMain_H
#define RunMain_H
#include "D1EulerSolver.h"
#include "InitialCDs.h"
#include "PhysicsValueAndControl.h"
#include "outputFile.h"
using namespace std;

int main()
{

    extern std::map<std::string, std::string> Para;
    
    readParaFrom_in_in_file();
   
    D1Euler D1_Euler_Solver(Para["whatflow"], 
                 Para["ReconstructType"],
               std::stod(Para["Time"]), std::stod(Para["cfl"]),
               std::stoi(Para["nx"]), std::stod(Para["epsilon"])
    );
    InitailCDs initialize(D1_Euler_Solver.Phy);
    if (std::stod(Para["Time"]) != 0)
    {
        D1_Euler_Solver.Phy.tf = std::stod(Para["Time"]);
    }

    auto t_s = std::chrono::high_resolution_clock::now();

    D1_Euler_Solver.Run();

    auto t_e = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(t_e - t_s).count();


    cout << "  ----------- \n \n \n   " << endl;
    cout << "CPU time of " << Para["whatflow"] << " at nx " << Para["nx"] << "= " << duration << "s" << endl; // ���ʱ�䣨��λ����
    cout << "  \n \n \n -----------   " << endl;

    outputFile *tp = new outputFile(D1_Euler_Solver.Phy);
    delete tp;    //
    tp = nullptr; //

    
    // for (const auto &pair : Para)
    // {
    //     std::cout << pair.first << "= " << pair.second << std::endl;
    // }
    
}
#endif // !RunMain_H
