/*      -------- 1D Finite: Element Class ----------- */
/*      -------- Everything put in the header file -- */

#include <iostream>
#include <stdio.h>

#include "Problem_functions.h"

using namespace std;

template <int K, int order>
class Element{
	// class containing basic info on what's happening inside an 1D element
        // K is the number of physical variables

	public:
		double center, size;
		int  dofs_number;

		double **dofs_backup;   // Temporary var for low storage RK4

		double **dofs;        //  dofs_number * (K physical vars) - unknowns
		double **dofs_fluxes; //  corresponding flux values

		double *flux_face[2]; // Left and Right flux interface values
		double *dofs_face[2]; // Left and Right dofs interface values 	
	               	
		unsigned int  quad_rule ; // number of quadrature points
		double * quad_weights;
		double * quad_points;

		// Element auxiliary methods
		double x_at_node(unsigned int node);
                double xi_at_node(unsigned int node);	
                double basis_function(unsigned int node, double xi);
                double basis_gradient(unsigned int node, double xi);
		void   save_dofs();
	
		// Test methods
		void print_elem();  
	
		// Main methods
		void setup(double x, double h, unsigned int  quad_rule, 
				double * quad_weights, double * quad_points);
		void compute_face_values();
		void compute_dof_fluxes();
		void step(double dt);
		

		Element();
		~Element();
	
};

template <int K, int order>
Element<K,order>::Element(){
	
	int i,k;

	dofs_number = order + 1;

        dofs        = new double*[dofs_number];
	dofs_fluxes = new double*[dofs_number];
	dofs_backup = new double*[dofs_number];

	for (i = 0; i < dofs_number; i++){
		dofs[i]        = new double[K];
		dofs_fluxes[i] = new double[K];
		dofs_backup[i]   = new double[K];
	}

	flux_face[0] = new double[K]; flux_face[1] = new double[K];
	dofs_face[0] = new double[K]; dofs_face[1] = new double[K];
 
	quad_points  = NULL;
	quad_weights = NULL;

};

template <int K, int order>
Element<K,order>::~Element(){

	delete [] quad_points;
	delete [] quad_weights;

	int i;
	for (i = 0; i < dofs_number; i++){
		delete [] dofs[i];
		delete [] dofs_fluxes[i];
		delete [] dofs_backup[i];		
	}

	delete [] flux_face[0]; delete [] flux_face[1];
};

template <int K, int order>
void Element<K,order>::setup(double x, double h,
		 unsigned int  quad_r, double * quad_pts, double * quad_wts){
	
	center = x;
	size = h;

	quad_rule        = quad_r;

	quad_points      = quad_pts;
	quad_weights     = quad_wts;
 
};

template <int K, int order>
double Element<K,order>::x_at_node(unsigned int node){
	
	return center + size/2.*xi_at_node(node);

}	

template <int K, int order>
double Element<K, order>::xi_at_node(unsigned int node){

	return quad_points[node]; // here the degrees of freedom are exactly the Gauss points
}


template <int K, int order>
void Element<K,order>::print_elem(){

	cout << " Center at x =  " << center << endl;

	int node;
	for (node = 0; node < dofs_number; node++){
		cout << " Node " << node << " at xi = " << xi_at_node(node) 
                << ",  x = " << x_at_node(node) << " -- quad_weight " <<
		quad_weights[node] << endl; 

	}
	cout << endl << " ----------------- " << endl;
};



template <int K, int order>
double Element<K,order>::basis_function(unsigned int node, double xi){
	/* Lagrange basis functions */	

	int i; double out = 1.;
	double xi_ref = xi_at_node(node);
	double xi_temp;

	for (i = 0; i < dofs_number; i++ ){
		if (i != node){
		xi_temp = xi_at_node(i);
		out *= (xi - xi_temp) / (xi_ref - xi_temp);
		}
	}
		
	return out;
}

template <int K, int order>
double Element<K,order>::basis_gradient(unsigned int node, double xi){
	/* Gradient of Lagrange Basis functions */ 

	int i,j; double denom = 1.; double temp = 1.; double out = 0.;
	double xi_ref = xi_at_node(node);
	double xi_temp;

	// Compute the common denominator first
        for (i = 0; i < dofs_number; i++){
		if (i != node){
		xi_temp = xi_at_node(i);
		denom *= 1 / (xi_ref - xi_temp);
		}
	}
	
	// Do the rest	
        for (j = 0; j < dofs_number; j++){
		if (j != node){
			for (i = 0; i < dofs_number; i++){
				if ((i != node)&&(i != j)){	
				xi_temp = xi_at_node(i);
				temp *= (xi - xi_temp);
				}
			}
		out += temp;
		temp = 1.;
		}
	}
	
	return out*denom;
}


template<int K, int order>
void Element<K,order>::compute_dof_fluxes(){
	int node;
	for (node = 0; node < dofs_number; node++){
		flux_val(dofs_fluxes[node],K,dofs[node]);
	}
}

template<int K, int order>
void Element<K,order>::compute_face_values(){
	/* Compute the interface values in preparation for the Riemann solver */
	
	int node, k;
	
	for (k = 0; k < K; k++){

		dofs_face[0][k] = 0.;
		dofs_face[1][k] = 0.;
		
		for (node = 0; node < dofs_number; node++){
			dofs_face[0][k] += dofs[node][k]*basis_function(node,-1.);
			dofs_face[1][k] += dofs[node][k]*basis_function(node,1.);
		}
	}

	//cout << "LEFT: h = " << dofs_face[0][0] << ", hu = " << dofs_face[0][1] << endl;	
	//cout << "RIGHT: h = " << dofs_face[1][0] << ", hu = " << dofs_face[1][1] << endl;
}


template<int K, int order>
void Element<K,order>::step(double dt){
	/* Euler explicit */
	
	compute_dof_fluxes();	

	int node, k, p;
	
	for (node = 0; node < dofs_number; node++){
		for (k = 0; k < K; k++){

		dofs[node][k] = dofs_backup[node][k]  -	 2*dt/(quad_weights[node]*size) * 
					( (flux_face[1][k])*(basis_function(node,1)) - 
					 (flux_face[0][k])*(basis_function(node,-1)) );

			for (p = 0; p < quad_rule; p++){
                                 dofs[node][k] += 2*dt/(quad_weights[node]*size)
					*(quad_weights[p])*(dofs_fluxes[p][k]) *
							(basis_gradient(node,quad_points[p]));

			}

		}
	} 
}

template <int K, int order>
void Element<K,order>::save_dofs(){
	
	int i, k;
	for (i = 0; i < dofs_number; i++){
		for (k = 0; k < K; k++){
			dofs_backup[i][k] = dofs[i][k];
		}
	}

}
