#ifndef D1_SpacialScheme_CPP
#define D1_SpacialScheme_CPP
#include "D1EulerSolver.h"

void D1Euler::weno5(double *__restrict__ splitF, int i, int sign, double &ReconstructedValue, double *__restrict__ Q, double *__restrict__ w, const char flag)
{
	// flag-> 'l':only the linear reconstruction on each substencil
	//        'w':only the weights
	//        'r':the whole WENOJS procedure
	int OneSign = sign;
	int TwoSign = 2 * sign;

	double f_im2 = splitF[i - TwoSign];
	double f_im1 = splitF[i - OneSign];
	double f_i = splitF[i];
	double f_ip1 = splitF[i + OneSign];
	double f_ip2 = splitF[i + TwoSign];

	if ('l' == flag || 'r' == flag)
	{
		Q[0] = 1.0 / 3.0 * f_im2 - 7.0 / 6.0 * f_im1 + 11.0 / 6.0 * f_i;
		Q[1] = -1.0 / 6.0 * f_im1 + 5.0 / 6.0 * f_i + 1.0 / 3.0 * f_ip1;
		Q[2] = 1.0 / 3.0 * f_i + 5.0 / 6.0 * f_ip1 - 1.0 / 6.0 * f_ip2;
	}

	if ('w' == flag || 'r' == flag)
	{
		double temp1 = f_im2 - 4.0 * f_im1 + 3.0 * f_i;
		double temp2 = f_im2 - 2.0 * f_im1 + f_i;
		double beta1 = 0.25 * temp1 * temp1 + 13.0 / 12.0 * temp2 * temp2;

		temp1 = f_im1 - f_ip1;
		temp2 = f_im1 - 2.0 * f_i + f_ip1;
		double beta2 = 0.25 * temp1 * temp1 + 13.0 / 12.0 * temp2 * temp2;

		temp1 = 3.0 * f_i - 4.0 * f_ip1 + f_ip2;
		temp2 = f_i - 2.0 * f_ip1 + f_ip2;
		double beta3 = 0.25 * temp1 * temp1 + 13.0 / 12.0 * temp2 * temp2;

		double beta1PlusEps = (beta1 + Phy.eps);
		double beta2PlusEps = (beta2 + Phy.eps);
		double beta3PlusEps = (beta3 + Phy.eps);
		double alpha1 = 0.1 / (beta1PlusEps * beta1PlusEps);
		double alpha2 = 0.6 / (beta2PlusEps * beta2PlusEps);
		double alpha3 = 0.3 / (beta3PlusEps * beta3PlusEps);

		double sum = alpha1 + alpha2 + alpha3;

		w[0] = alpha1 / sum;
		w[1] = alpha2 / sum;
		w[2] = alpha3 / sum;
	}
	if ('r' == flag)
		ReconstructedValue = w[0] * Q[0] + w[1] * Q[1] + w[2] * Q[2];
}

#endif
