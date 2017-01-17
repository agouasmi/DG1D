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

	double p_sec = 0.002; /* Pause time*/ 
	double dt = 0.001;	
	double Nt = 5000;
	
	const int K = 2, order = 2 , N_E = 100;
	DG_1D<K,order> problem(N_E);

	int npts = N_E*(order + 1);

	double out[npts];
	double X[npts];	

	problem.read_global_sol();
	
	for (int i = 0; i < npts; i++){
		X[i]   = problem.domain[i];
	}

	problem.write_sol(out,0);
	
	// Plotting with DISLIN

	Dislin g;

	plot_init_setup(g); 
	plot_axis_setup(g);	
	
	double min, max;

	g.graf(X[0],X[npts-1],0.,5,0.8,1.2,0.8,0.2);
	g.title();
	g.curve(X,out,npts);
	g.endgrf();

	for (int i = 0; i < Nt; i++){

		problem.advance_RK(dt,3);
		problem.read_global_sol();
		problem.write_sol(out,0);	
	
		pause(p_sec);
		g.erase();		
		
		string time =  "Time t = ";
		time += std::to_string((float) problem.t);	
		
		set_title(g,time);
				
		//min = mini(out,npts);
		//max = maxi(out,npts);

		//cout << min << " " << max << endl;	

		g.graf(X[0],X[npts-1], 0. , 5  , 0.95, 1.15, 1.0, 0.05); 
		g.title();
	
		g.curve(X,out,npts);
		
		g.endgrf();
	}

	g.disfin(); 

	return 0;
}
