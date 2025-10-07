#ifndef D1_Aug_CPP
#define D1_Aug_CPP
#include "D1EulerSolver.h"

void D1Euler::ChooseReconstructionMethod_AndGo()
{

	if (Phy.ReconstructType == "Comp")
	{
		Component_Wise();
	}
	else if (Phy.ReconstructType == "Char")
	{

		Characteristic_Wise();
	}
	else if (Phy.ReconstructType == "PS")
	{
		PressureEntropyReplacement();
	}
	else if (Phy.ReconstructType == "Co")
	{
		ImprovedCo();
	}
	else {
		std::cout<<"\n Please choose a right reconstruction-type name in \' in.in \' file, or define a new function using this name. \n"<<std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void D1Euler::Component_Wise()
{
	// std::cout<<" Im here in the flux reconstruction methods :  Component-Wise" << std::endl;

	double Reconstructed_FP_dens = 0, Reconstructed_FP_momt = 0, Reconstructed_FP_eng = 0;
	double Reconstructed_FM_dens = 0, Reconstructed_FM_momt = 0, Reconstructed_FM_eng = 0;
	double Q[10], W[10]; // Q and W are just placeholders to avoid passing nullptr. // They are not used in the later part of this function.

	int boundL = Phy.Ghost - 1, boundR = Phy.NxAug - Phy.Ghost;
	for (int i = boundL; i < boundR; i++)
	{
		//--------------UpWind
		// f^+_{i+1/2} = WENO5(f_{i-2,...i+2})
		weno5(Phy.FP.dens, i, 1, Reconstructed_FP_dens, Q, W, 'r');
		weno5(Phy.FP.momt, i, 1, Reconstructed_FP_momt, Q, W, 'r');
		weno5(Phy.FP.eng, i, 1, Reconstructed_FP_eng, Q, W, 'r');
		//--------------DownWind
		// f^-_{i+1/2} = WENO5(f_{i-1,...i+3})
		weno5(Phy.FM.dens, i + 1, -1, Reconstructed_FM_dens, Q, W, 'r');
		weno5(Phy.FM.momt, i + 1, -1, Reconstructed_FM_momt, Q, W, 'r');
		weno5(Phy.FM.eng, i + 1, -1, Reconstructed_FM_eng, Q, W, 'r');

		// Calculate Flux_interface-Begin
		Phy.Flux_interface.dens[i] = Reconstructed_FP_dens + Reconstructed_FM_dens;
		Phy.Flux_interface.momt[i] = Reconstructed_FP_momt + Reconstructed_FM_momt;
		Phy.Flux_interface.eng[i] = Reconstructed_FP_eng + Reconstructed_FM_eng;
		// Calculate Flux_interface-End
	}
}

void D1Euler::Characteristic_Wise()
{

	double Ch1P[10] = {0}, Ch2P[10] = {0}, Ch3P[10] = {0}; // Charcteristic variables for F^+
	double Ch1M[10] = {0}, Ch2M[10] = {0}, Ch3M[10] = {0}; // For F^-

	double Ch1P_reconstructed = 0, Ch2P_reconstructed = 0, Ch3P_reconstructed = 0;
	double Ch1M_reconstructed = 0, Ch2M_reconstructed = 0, Ch3M_reconstructed = 0;

	int boundL = Phy.Ghost - 1, boundR = Phy.NxAug - Phy.Ghost;
	for (int i = boundL; i < boundR; i++)
	{

		CalculateCharactMatrices(i, 'f');

		int c_L = Phy.Ghost - 1, c_R = Phy.Ghost;
		// Local characteristic projections for L \dot F^+_{i-2,...i+2}
		for (int c = -c_L; c < c_R; c++)
		{
			int i_c = i + c;
			Ch1P[2 + c] = Phy.FP.dens[i_c] * Phy.L[1] + Phy.FP.momt[i_c] * Phy.L[2] + Phy.FP.eng[i_c] * Phy.L[3];
			Ch2P[2 + c] = Phy.FP.dens[i_c] * Phy.L[4] + Phy.FP.momt[i_c] * Phy.L[5] + Phy.FP.eng[i_c] * Phy.L[6];
			Ch3P[2 + c] = Phy.FP.dens[i_c] * Phy.L[7] + Phy.FP.momt[i_c] * Phy.L[8] + Phy.FP.eng[i_c] * Phy.L[9];
		}
		// Local characteristic projections for L \dot F^-_{i-1,...i+3}
		for (int c = -c_L; c < c_R; c++)
		{
			int i_c = i + c + 1;
			Ch1M[2 + c] = Phy.FM.dens[i_c] * Phy.L[1] + Phy.FM.momt[i_c] * Phy.L[2] + Phy.FM.eng[i_c] * Phy.L[3];
			Ch2M[2 + c] = Phy.FM.dens[i_c] * Phy.L[4] + Phy.FM.momt[i_c] * Phy.L[5] + Phy.FM.eng[i_c] * Phy.L[6];
			Ch3M[2 + c] = Phy.FM.dens[i_c] * Phy.L[7] + Phy.FM.momt[i_c] * Phy.L[8] + Phy.FM.eng[i_c] * Phy.L[9];
		}
		//--------------
		double Q[10], W[10]; // Q and W are just placeholders to avoid passing nullptr. // They are not used in the later part of this function.
		weno5(Ch1P, 2, 1, Ch1P_reconstructed, Q, W, 'r');
		weno5(Ch2P, 2, 1, Ch2P_reconstructed, Q, W, 'r');
		weno5(Ch3P, 2, 1, Ch3P_reconstructed, Q, W, 'r');
		//--------------
		weno5(Ch1M, 2, -1, Ch1M_reconstructed, Q, W, 'r');
		weno5(Ch2M, 2, -1, Ch2M_reconstructed, Q, W, 'r');
		weno5(Ch3M, 2, -1, Ch3M_reconstructed, Q, W, 'r');

		double temp1 = Ch1P_reconstructed + Ch1M_reconstructed;
		double temp2 = Ch2P_reconstructed + Ch2M_reconstructed;
		double temp3 = Ch3P_reconstructed + Ch3M_reconstructed;

		Phy.Flux_interface.dens[i] = temp1 * Phy.R[1] + temp2 * Phy.R[2] + temp3 * Phy.R[3];
		Phy.Flux_interface.momt[i] = temp1 * Phy.R[4] + temp2 * Phy.R[5] + temp3 * Phy.R[6];
		Phy.Flux_interface.eng[i] = temp1 * Phy.R[7] + temp2 * Phy.R[8] + temp3 * Phy.R[9];
	}
}

void D1Euler::ImprovedCo()
{

	double *AvU = Phy.RoeAv.AvU;

	double *u = Phy.VarWorking.u;
	double *p = Phy.VarWorking.p;
	double *rho = Phy.VarWorking.dens;

	double phi_P[10], phi_M[10];
	double *Psi = p;

	int boundL = Phy.Ghost - 1, boundR = Phy.NxAug - Phy.Ghost;
	for (int i = boundL; i < boundR; i++)
	{

		CalculateCharactMatrices(i, 'l');

		// Calculate the characteristic variable associated with the linearly degenerated field.--Begin
		int c_L = Phy.Ghost - 1, c_R = Phy.Ghost;
		for (int c = -c_L; c < c_R; c++)
		{ // for F^+_{i-2,...i+2}
			int iPc = i + c;
			phi_P[2 + c] = Phy.FP.dens[iPc] * Phy.L[4] + Phy.FP.momt[iPc] * Phy.L[5] + Phy.FP.eng[iPc] * Phy.L[6];
		}
		for (int c = -c_L; c < c_R; c++)
		{ // for  F^-_{i-1,...i+3}
			int iPc = i + 1 + c;
			phi_M[2 + c] = Phy.FM.dens[iPc] * Phy.L[4] + Phy.FM.momt[iPc] * Phy.L[5] + Phy.FM.eng[iPc] * Phy.L[6];
		}
		// Calculate the characteristic variable associated with the linearly degenerated field.--End

		double Q_F1P[5], Q_F2P[5], Q_F3P[5];
		double Q_F1M[5], Q_F2M[5], Q_F3M[5];
		double Q_phiP[5];
		double Q_phiM[5];

		{ // Get the substencil linear reconstructions
			double rc;
			double W[10];
			weno5(Phy.FP.dens, i, 1, rc, Q_F1P, W, 'l');
			weno5(Phy.FP.momt, i, 1, rc, Q_F2P, W, 'l');
			weno5(Phy.FP.eng, i, 1, rc, Q_F3P, W, 'l');

			weno5(phi_P, 2, 1, rc, Q_phiP, W, 'l');

			weno5(Phy.FM.dens, i + 1, -1, rc, Q_F1M, W, 'l');
			weno5(Phy.FM.momt, i + 1, -1, rc, Q_F2M, W, 'l');
			weno5(Phy.FM.eng, i + 1, -1, rc, Q_F3M, W, 'l');

			weno5(phi_M, 2, -1, rc, Q_phiM, W, 'l');
		}

		double w_Psi_P[5], w_Psi_M[5];
		double w_phi_P[5], w_phi_M[5];
		{ // Get the weights used (\omega(\psi)_k and \omega(\phi)_k)
			double placeholder2[10], placeholder1;
			weno5(Psi, i, 1, placeholder1, placeholder2, w_Psi_P, 'w');
			weno5(Psi, i + 1, -1, placeholder1, placeholder2, w_Psi_M, 'w');

			weno5(phi_P, 2, 1, placeholder1, placeholder2, w_phi_P, 'w');
			weno5(phi_M, 2, -1, placeholder1, placeholder2, w_phi_M, 'w');
		}

		// Calculte the kept term
		double avu = AvU[i];
		double halfAvUU = 0.5 * avu * avu;
		double Kept_Temrm = 0.0;
		for (int j = 0; j < 3; ++j)
		{ //                    Kept term for    F^+             and   kept term for  // F^-
			Kept_Temrm += (Q_phiP[j] * (w_phi_P[j] - w_Psi_P[j])) + (Q_phiM[j] * (w_phi_M[j] - w_Psi_M[j]));
		}
		// Calculte the kept term

		// Calculate the first term Comonmon-Weighted fluxes
		double Reconstructed_F_interface_dens = 0, Reconstructed_F_interface_momt = 0, Reconstructed_F_interface_eng = 0;
		for (int j = 0; j < 3; ++j)
		{ //                                  First term for  F^+   and   first term for  F^-
			Reconstructed_F_interface_dens += (Q_F1P[j] * w_Psi_P[j]) + (Q_F1M[j] * w_Psi_M[j]);
			Reconstructed_F_interface_momt += (Q_F2P[j] * w_Psi_P[j]) + (Q_F2M[j] * w_Psi_M[j]);
			Reconstructed_F_interface_eng += (Q_F3P[j] * w_Psi_P[j]) + (Q_F3M[j] * w_Psi_M[j]);
		}
		// Calculate the first term

		// Add up the two terms
		Reconstructed_F_interface_dens += Kept_Temrm;
		Reconstructed_F_interface_momt += avu * Kept_Temrm;
		Reconstructed_F_interface_eng += halfAvUU * Kept_Temrm;
		// Add up the two terms

		Phy.Flux_interface.dens[i] = Reconstructed_F_interface_dens;
		Phy.Flux_interface.momt[i] = Reconstructed_F_interface_momt;
		Phy.Flux_interface.eng[i] = Reconstructed_F_interface_eng;
	}
}

void D1Euler::PressureEntropyReplacement()
{
	// Efficient Implementation of Weighted ENO Schemes.G.S.Jiang & C.W.Shu.1996. Equation(4.6);

	double *AvU = Phy.RoeAv.AvU;
	double *AvA = Phy.RoeAv.AvA;

	double *u = Phy.VarWorking.u;
	double *p = Phy.VarWorking.p;
	double *rho = Phy.VarWorking.dens;
	double *E = Phy.VarWorking.eng;

	// double Entropy[4000] = {0};
	double* Entropy = Phy.VarWorking.Entropy;

	for (int i = 0; i < Phy.NxAug; i++)
	{
		// Entropy[i] = log(p[i] / pow(Phy.VarWorking.dens[i], Phy.gamma));
		Entropy[i] = (p[i] / pow(rho[i], Phy.gamma)); // This way seems more robust
	}

	int boundL = Phy.Ghost - 1, boundR = Phy.NxAug - Phy.Ghost;
	for (int i = boundL; i < boundR; i++)
	{

		CalculateCharactMatrices(i, 'n');

		double L_F1P[5], L_F2P[5], L_F3P[5]; // Linear reconstruction of the split fluxes on each sub-stencil.--UpWinding
		double L_F1M[5], L_F2M[5], L_F3M[5]; // --DownWinding

		{ // Get the linear reconstructions
			double rc, W[10];
			weno5(Phy.FP.dens, i, 1, rc, L_F1P, W, 'l');
			weno5(Phy.FP.momt, i, 1, rc, L_F2P, W, 'l');
			weno5(Phy.FP.eng, i, 1, rc, L_F3P, W, 'l');

			weno5(Phy.FM.dens, i + 1, -1, rc, L_F1M, W, 'l');
			weno5(Phy.FM.momt, i + 1, -1, rc, L_F2M, W, 'l');
			weno5(Phy.FM.eng, i + 1, -1, rc, L_F3M, W, 'l');
		}

		double wP_P[5], wpP_M[5];
		double wEtr_P[5], wEtr_M[5];

		{ // Get the weights used
			double rc, Q[5];
			weno5(p, i, 1, rc, Q, wP_P, 'w');
			weno5(p, i + 1, -1, rc, Q, wpP_M, 'w');

			weno5(Entropy, i, 1, rc, Q, wEtr_P, 'w');
			weno5(Entropy, i + 1, -1, rc, Q, wEtr_M, 'w');
		}

		double FtP[5]; // vassel temporary reconstruct f1 f2 f3
		double FtM[5];
		int loopTemp = 3;

		// Calculate the first and second terms of Eq.(4.6), projections and inverse projections
		for (size_t h = 0; h < loopTemp; h++)
		{
			// (F^+_{j+1/2,1/m}-F^+_{j+1/2,2})
			double wt_P = (wP_P[h] - wEtr_P[h]);
			FtP[1] += L_F1P[h] * wt_P;
			FtP[2] += L_F2P[h] * wt_P;
			FtP[3] += L_F3P[h] * wt_P;
			// (F^-_{j+1/2,1/m}-F^-_{j+1/2,2})
			double wt_M = (wpP_M[h] - wEtr_M[h]);
			FtM[1] += L_F1M[h] * wt_M;
			FtM[2] += L_F2M[h] * wt_M;
			FtM[3] += L_F3M[h] * wt_M;
		}
		//---
		// ForFirstTerm = l_1 \cdot (F_{j+1/2,1/m}-F_{j+1/2,2})
		double ForFirstTerm = (FtP[1] + FtM[1]) * Phy.L[1] + (FtP[2] + FtM[2]) * Phy.L[2] + (FtP[3] + FtM[3]) * Phy.L[3];
		// ForSecondTerm = l_m \cdot (F_{j+1/2,1/m}-F_{j+1/2,2})
		double ForSecondTerm = (FtP[1] + FtM[1]) * Phy.L[7] + (FtP[2] + FtM[2]) * Phy.L[8] + (FtP[3] + FtM[3]) * Phy.L[9];

		
		double Reconstructed_F_interface_dens = 0, Reconstructed_F_interface_momt = 0, Reconstructed_F_interface_eng = 0;
      //                                  ForFirstTerm R_1       and  ForSecondTerm R_m  
		Reconstructed_F_interface_dens += ForFirstTerm * Phy.R[1] + ForSecondTerm * Phy.R[3];
		Reconstructed_F_interface_momt += ForFirstTerm * Phy.R[4] + ForSecondTerm * Phy.R[6];
		Reconstructed_F_interface_eng += ForFirstTerm * Phy.R[7] + ForSecondTerm * Phy.R[9];
		// Calculate the first and second terms of Eq.(4.6), projections and inverse projections

		// Adding up the last term of Eq.(4.6)
		for (size_t h = 0; h <= loopTemp; h++)
		{
			Reconstructed_F_interface_dens += (L_F1P[h] * wEtr_P[h]) + (L_F1M[h] * wEtr_M[h]);
			Reconstructed_F_interface_momt += (L_F2P[h] * wEtr_P[h]) + (L_F2M[h] * wEtr_M[h]);
			Reconstructed_F_interface_eng += (L_F3P[h] * wEtr_P[h]) + (L_F3M[h] * wEtr_M[h]);
		}
		// Adding up the last term of Eq.(4.6)

		Phy.Flux_interface.dens[i] = Reconstructed_F_interface_dens;
		Phy.Flux_interface.momt[i] = Reconstructed_F_interface_momt;
		Phy.Flux_interface.eng[i] = Reconstructed_F_interface_eng;
	}
}
#endif
