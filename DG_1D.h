/*       -------- 1D Discontinuous Galerkin Class  --------- */

#include <iostream>
#include <stdio.h>
#include <math.h>

#include "Element.h"

using namespace std;

template <int K, int order>
class DG_1D{
	
	public:
		// Constructor and destructor
                DG_1D(int num_elems);
                ~DG_1D();
	
		double t;
		
		int    N_E;                    // number of elements
		int    global_dof_num;         // total number of dofs = (N_E * quad_rule)

		double LEFT_BC, RIGHT_BC;      // Boundaries
					
		unsigned int quad_rule;        // number of (Gauss-Legendre) quadrature points
		double *  quad_weights;        // quadrature weights
		double *   quad_points;        // quadrature points
		
		Element<K,order> * elements;   // 1D Finite Elements
	
		double **  solution;	       // Discontinuous solution - (global_dof_num) * K
		double *     domain;           // Physical domain
		
		// Solution steps
		void compute_fluxes();
		void advance(double dt);
		void advance_RK(double dt, int RK);

		// Auxiliary functions
                void print_elements();
		void read_global_sol();
		void write_sol(double * dest, int k);		
		
	private:	
		void set_quad_rule();
		void setup_system();           // set quad points, quad weights, initialize elements
		void set_initial_conditions();		
};


template <int K, int order>
DG_1D<K, order>::DG_1D(int num_elements){	
	
	N_E = num_elements;
	t = 0;
	
	quad_rule = order + 1;

	global_dof_num = N_E*quad_rule;

	solution = new double*[global_dof_num];

	quad_weights = new double[quad_rule];
	quad_points  = new double[quad_rule];	
	
	LEFT_BC = 0.; RIGHT_BC = 10.;
	
	for (int i = 0; i < global_dof_num; i++){
		solution[i] = new double[K];
	}

	domain = new double[quad_rule*N_E];	

	elements = new Element<K,order>[N_E];	

	setup_system();
	set_initial_conditions();	
};

template <int K, int order>
DG_1D<K, order>::~DG_1D(){
	
	delete [] quad_weights;
	delete [] quad_points;
	delete [] domain;
	for (int i = 0; i < global_dof_num; i++){
		delete [] solution[i];
	}
	delete [] solution;
	
	delete [] elements;	

};

template <int K, int order>
void DG_1D<K, order>::print_elements(){

	int ELEM;
	for (ELEM = 0; ELEM < N_E; ELEM++){

		cout << " Element " << ELEM << endl;
		elements[ELEM].print_elem();
	} 
};

template <int K, int order>
void DG_1D<K,order>::set_quad_rule(){
	
	switch (quad_rule){
		case 1: 
		    	quad_points[0] = 0.; quad_weights[0] = 2.; 
			break;
		case 2:
			quad_points[0] = -sqrt(1./3.); quad_weights[0] = 1;
			quad_points[1] =  sqrt(1./3.); quad_weights[1] = 1;
			break;
		case 3:
			quad_points[0] = -sqrt(3./5.); quad_weights[0] = 5./9.;
			quad_points[1] =           0.; quad_weights[1] = 8./9.;
			quad_points[2] =  sqrt(3./5.); quad_weights[2] = 5./9.;
		        break;
		case 4:
			quad_points[0]  = -sqrt(3./7. + 2./7.*sqrt(6./5.));
			quad_weights[0] = (18. - sqrt(30.))/36.;
			quad_points[1]  = -sqrt(3./7. - 2./7.*sqrt(6./5.));
			quad_weights[1] =  (18. + sqrt(30.))/36.;

			quad_points[2] = -quad_points[1];
		        quad_weights[2] = quad_weights[1];
			quad_points[3] = -quad_points[0];
			quad_weights[3] = quad_weights[0];
			break;
		case 5:
			quad_points[0] = -1./3.*sqrt(5. + 2.*sqrt(10./7.));
			quad_points[1] = -1./3.*sqrt(5. - 2.*sqrt(10./7.));
			quad_points[2] = 0;
			quad_points[3] = -quad_points[1];
			quad_points[4] = -quad_points[0];

			quad_weights[0] = (322. - 13.*sqrt(70.))/900.;
			quad_weights[1] = (322. + 13.*sqrt(70.))/900.;
			quad_weights[2] = (128./225.);			
			quad_weights[3] = quad_weights[1];
			quad_weights[4] = quad_weights[0];
	}

};

template <int K, int order>
void DG_1D<K, order>::setup_system(){
	
	// Define the quadrature rule
	set_quad_rule();

	// Initialize the elements (center, size, dofs through the quad points
	int i, q;
	double h = (RIGHT_BC - LEFT_BC)/double(N_E);
	double x_c;
	for (i = 0; i < N_E; i++){
		x_c = LEFT_BC + h*i + h/2.;
		elements[i].setup(x_c, h, quad_rule,
					quad_points, quad_weights);

		for (q = 0; q < quad_rule; q++){
			domain[i*quad_rule + q] = x_c + h/2.*quad_points[q];
		}
	}

};

template <int K, int order>
void DG_1D<K, order>:: set_initial_conditions(){

	int ELEM, node, k;
	for (ELEM = 0; ELEM < N_E; ELEM++){
		for (node = 0; node < elements[ELEM].dofs_number; node++){
			for(k = 0; k < K; k++){

			initial_value(elements[ELEM].dofs[node],K, 
					elements[ELEM].x_at_node(node));
			}
		}
	}

};

template <int K, int order>
void DG_1D<K, order>::read_global_sol(){
		
	int ELEM, q;
	int k;	
	for (ELEM = 0; ELEM < N_E; ELEM++){
		for (k = 0; k < K; k++){
			for (q = 0; q < (elements[ELEM]).dofs_number; q++){
			solution[ELEM*quad_rule + q][k] = elements[ELEM].dofs[q][k];
			}
		}
	}

};

template <int K, int order>
void DG_1D<K, order>::compute_fluxes(){
	/* Compute the fluxes at the interfaces */

	int k, ELEM;
 
	/* Compute the dof values at the interfaces */  
	for (ELEM = 0; ELEM < N_E; ELEM++){
		elements[ELEM].compute_face_values();
	}
	
	/* Left and right boundary interfaces */
	wall(elements[0].dofs_face[0],K); 
	wall(elements[N_E-1].dofs_face[1],K);

	flux_val((elements[0]).flux_face[0],    K, (elements[0]).dofs_face[0]);
	flux_val((elements[N_E-1]).flux_face[1],K, (elements[N_E-1]).dofs_face[1]);

	/* Interior interfaces */
	for (ELEM = 0; ELEM < N_E-1; ELEM++){

		Roe_flux(elements[ELEM].flux_face[1],
					elements[ELEM].dofs_face[1],
						elements[ELEM+1].dofs_face[0], K);

		copy_vec(elements[ELEM+1].flux_face[0],
				elements[ELEM].flux_face[1], K);
	}

};

template <int K, int order>
void DG_1D<K, order>:: write_sol(double * dest, int k){
	
	int i;
	for (i = 0; i < global_dof_num; i++){
		dest[i] = solution[i][k];
	}

}

template <int K, int order>
void DG_1D<K, order>::advance(double dt){
		
	int ELEM;
	compute_fluxes();
	for (ELEM = 0; ELEM < N_E; ELEM++){
		elements[ELEM].save_dofs();
		elements[ELEM].step(dt);
	}

	t += dt;
}

template <int K, int order>
void DG_1D<K, order>::advance_RK(double dt, int RK){

	int ELEM, level;

	for (ELEM = 0; ELEM < N_E; ELEM++){
		elements[ELEM].save_dofs();
	}
		
	for (level = RK; level > 0; level--){
		compute_fluxes();
	
		for (ELEM = 0; ELEM < N_E; ELEM++){
			elements[ELEM].step(dt/(double) level);
		}
	}

	t+=dt;
}
