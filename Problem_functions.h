#include <math.h>
#include <cmath>
#include <iostream>

void initial_value(double * out, int K, double x);
void flux_val(double * out, int K, const double * in);
void Roe_flux(double * out, double * u_left, double * u_right, int K);
void wall(double * q, int K);
void copy_flux(double *out, double* in, int K);

using namespace std;

void initial_value(double* out, int K,double x){
	
	out[0] =  1. + 0.1*exp(-pow((x-5.),2));
	out[1] =  0.;

}

void flux_val(double * out, int K, const double * in){
	double q1 = in[0];
	double q2 = in[1];

	out[0] = q2;
	out[1] = pow(q2,2)/q1 + 9.8/2.*pow(q1,2);

} 

void Roe_flux(double * out, double * q_left, double * q_right, int K){

	int k;
	double FL[K]; double FR[K]; double q_av[K]; double x[K];
	
	flux_val(FL,K,q_left);
	flux_val(FR,K,q_right);

	for (k = 0; k < K; k++){
		q_av[k] = 0.5*(q_left[k] + q_right[k]);
	}
	
	double h_av = q_av[0]; double u_av = q_av[1]/q_av[0]; 
	double L1 = u_av + sqrt(9.8*h_av); double L2 = u_av - sqrt(9.8*h_av);

	double a = 2.*sqrt(9.8*h_av);

	for (k = 0; k < K; k++){
		x[k] = (q_right[k] - q_left[k]);
	}

	out[0] = 0.5*(FL[0] + FR[0]) - 0.5 / a * ( abs(L1) * (-L2*x[0] + x[1]) + 
							abs(L2) * (L1*x[0] - x[1]) );

	out[1] = 0.5*(FL[1] + FR[1]) - 0.5 / a * ( abs(L1) * (-L2*L1*x[0] + L1*x[1]) +
							abs(L2) * (L2*L1*x[0] - L2*x[1]) );

}

void wall(double * q, int K){
	/* Functions that defines the boundary value*/
	q[1] = 0;
}

void copy_vec(double *out, double* in, int N){
	int n;
	for (n = 0; n < N; n++){
		out[n] = in[n];
	}
}
