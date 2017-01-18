#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <thread>

#include "DG_1D.h"
#include "discpp.h"
#include "anim.h"

using namespace std;

/* --- Main program --- */

double maxi(double * vec, int N){
	return *std::max_element(vec,vec+N);
};

double mini(double * vec, int N){
	return *std::min_element(vec,vec+N);
};

int main(int argc, char **argv){

	double p_sec = 0.005; /* Pause time*/ 
	double dt = 0.002;
	double T = 4.;	
	double Nt = int (T/dt) + 1;
	
	const int K = 2,  N_E = 32;

	DG_1D<K,0> problem_0(N_E);
	DG_1D<K,1> problem_1(N_E);
	DG_1D<K,2> problem_2(N_E);
	DG_1D<K,3> problem_3(N_E);

	int npts_0 = N_E*(0 + 1);
	int npts_1 = N_E*(1 + 1);
	int npts_2 = N_E*(2 + 1);
	int npts_3 = N_E*(3 + 1);

	problem_0.read_global_sol();
	problem_1.read_global_sol();
	problem_2.read_global_sol();
	problem_3.read_global_sol();
	
	/*for (int i = 0; i < npts; i++){
		X[i]   = problem.domain[i];
	}*/

	double out_0[npts_0];
	double out_1[npts_1];
	double out_2[npts_2];
	double out_3[npts_3];

	//problem.write_sol(out,0);
	
	// Plotting with DISLIN

	Dislin g;

	plot_init_setup(g); 
	plot_axis_setup(g);	
	
	double min, max;

	/*g.graf(X[0],X[npts-1],0.,5,0.8,1.2,0.8,0.2);
	g.title();
	g.curve(X,out,npts);
	g.endgrf();*/

	for (int i = 0; i < Nt; i++){

		problem_0.advance_RK(dt,4);
		problem_1.advance_RK(dt,4);
		problem_2.advance_RK(dt,4);
		problem_3.advance_RK(dt,4);

		problem_0.read_global_sol();
		problem_1.read_global_sol();
		problem_2.read_global_sol();
		problem_3.read_global_sol();
		
		problem_0.write_sol(out_0,0);	
		problem_1.write_sol(out_1,0);
		problem_2.write_sol(out_2,0);
		problem_3.write_sol(out_3,0);
	
		pause(p_sec);
		g.erase();		
		
		string time =  "Time t = ";
		time += std::to_string((float) problem_0.t);	
		
		set_title(g,time);
				
		//min = mini(out,npts);
		//max = maxi(out,npts);

		//cout << min << " " << max << endl;	

		g.graf(0.,10., 0. , 5  , 0.95, 1.15, 1.0, 0.05); 
		g.title();
	
		g.curve(problem_0.domain,out_0,npts_0);
		g.curve(problem_1.domain,out_1,npts_1);
		g.curve(problem_2.domain,out_2,npts_2);
		g.curve(problem_3.domain,out_3,npts_3);

		g.endgrf();
	}

	g.disfin(); 

	return 0;
}
