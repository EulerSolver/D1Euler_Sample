#ifndef D1_TVDRK_CPP
#define D1_TVDRK_CPP
#include "D1EulerSolver.h"

void D1Euler::TVD_RK3(int Level)
{

    boundaryConditionsCall();

    get_u_p_a();

    GLF_FVS();  //Global L-F splitting

    //Calculate the fluxes at each interface
    //Four choices: 1.Component-Wise, 2:Characteristic-wise, 3. Pressure and  entropy replacement, 4.Improved common-weights.
    ChooseReconstructionMethod_AndGo();  
    //Calculate the fluxes at each interface

    {  
        int i_Minus_1 = 0;
        int boundL = Phy.Ghost, boundR = Phy.NxAug - Phy.Ghost;
        for (int i = boundL; i < boundR; i++)
        {
            i_Minus_1 = i - 1;

            Phy.RHS.dens[i] = -(Phy.Flux_interface.dens[i] - Phy.Flux_interface.dens[i_Minus_1]) / Phy.dx;
            Phy.RHS.momt[i] = -(Phy.Flux_interface.momt[i] - Phy.Flux_interface.momt[i_Minus_1]) / Phy.dx;
            Phy.RHS.eng[i] = -(Phy.Flux_interface.eng[i] - Phy.Flux_interface.eng[i_Minus_1]) / Phy.dx;
        }
    }

    if (Level == 1)
    {
        for (int i = 0; i < Phy.NxAug; i++)
        {
            Phy.VarRK1.dens[i] = Phy.Var.dens[i] + Phy.dt * Phy.RHS.dens[i];
            Phy.VarRK1.momt[i] = Phy.Var.momt[i] + Phy.dt * Phy.RHS.momt[i];
            Phy.VarRK1.eng[i] = Phy.Var.eng[i] + Phy.dt * Phy.RHS.eng[i];
        }
        //------------------------------
        Phy.VarWorking = Phy.VarRK1;
        //------------------------------
        TVD_RK3(2);
    }
    else if (Level == 2)
    {
        for (int i = 0; i < Phy.NxAug; i++)
        {
            Phy.VarRK2.dens[i] = 3.0 / 4.0 * Phy.Var.dens[i] + 1.0 / 4.0 * Phy.VarRK1.dens[i] + 1.0 / 4.0 * Phy.dt * Phy.RHS.dens[i];
            Phy.VarRK2.momt[i] = 3.0 / 4.0 * Phy.Var.momt[i] + 1.0 / 4.0 * Phy.VarRK1.momt[i] + 1.0 / 4.0 * Phy.dt * Phy.RHS.momt[i];
            Phy.VarRK2.eng[i] = 3.0 / 4.0 * Phy.Var.eng[i] + 1.0 / 4.0 * Phy.VarRK1.eng[i] + 1.0 / 4.0 * Phy.dt * Phy.RHS.eng[i];
        }
        Phy.VarWorking = Phy.VarRK2;
        TVD_RK3(3);
    }
    else if (Level == 3)
    {
        for (int i = 0; i < Phy.NxAug; i++)
        {
            Phy.Var.dens[i] = 1.0 / 3.0 * Phy.Var.dens[i] + 2.0 / 3.0 * Phy.VarRK2.dens[i] + 2.0 / 3.0 * Phy.dt * Phy.RHS.dens[i];
            Phy.Var.momt[i] = 1.0 / 3.0 * Phy.Var.momt[i] + 2.0 / 3.0 * Phy.VarRK2.momt[i] + 2.0 / 3.0 * Phy.dt * Phy.RHS.momt[i];
            Phy.Var.eng[i] = 1.0 / 3.0 * Phy.Var.eng[i] + 2.0 / 3.0 * Phy.VarRK2.eng[i] + 2.0 / 3.0 * Phy.dt * Phy.RHS.eng[i];
        }
        Phy.VarWorking = Phy.Var;
    }
}

#endif