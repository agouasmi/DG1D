#include "DG_1D.h"
#include "Problem_functions.h"

template <int K, int order>
void DG_1D<K, order>::print_elements(){

	int ELEM;
	for (ELEM = 0; ELEM < N_E; ELEM++){

		cout << " Element " << ELEM << endl;
		(this->elements[ELEM]).print_elem();
	} 

}

template <int K, int order>
void DG_1D<K,order>::set_quad_rule(){
	
	switch (quad_rule){
		case 1: 
		    	this->quad_points[0] = 0.; this->quad_weights[0] = 2.; 
			break;
		case 2:
			this->quad_points[0] = -sqrt(1./3.); this->quad_weigths[0] = 1;
			this->quad_points[1] =  sqrt(1./3.); this->quad_weights[1] = 1;
			break;
		case 3:
			this->quad_points[0] = -sqrt(3./5.); this->quad_weigths[0] = 5./9.;
			this->quad_points[1] =           0.; this->quad_weigths[1] = 8./9.;
			this->quad_points[2] =  sqrt(3./5.); this->quad_weigths[2] = 5./9.;
		        break;
		case 4:
			this->quad_points[0]  = -sqrt(3./5. + 2./7.*sqrt(6./5.));
			this->quad_weigths[0] = (18. - sqrt(30.))/36.;
			this->quad_points[1]  = -sqrt(3./5. - 2./7.*sqrt(6./5.));
			this->quad_weights[1] =  (18. + sqrt(30.))/36.;

			this->quad_points[2] = -this->quad_points[1];
		        this->quad_weights[2] = this->quad_weights[1];
			this->quad_points[3] = -this->quad_points[0];
			this->quad_weights[3] = this->quad_weights[0];
			break;			
	}

}

template <int K, int order>
void DG_1D<K, order>::setup_system(){
	
	// Define the quadrature rule
	this->set_quad_rule();

	// Initialize the elements (center, size, dofs through the quad points
	int i;
	double h = 1./double(N_E);
	double x_c;
	for (i = 0; i < N_E; i++){
		x_c = h*i + h/2.;
		(this->elements[i]).setup(x_c, h, &(this->quad_rule),
					this->quad_points, this->quad_weights);
	}

}

template <int K, int order>
void DG_1D<K, order>:: set_initial_conditions(){

	int i, node;
	for (i = 0; i < N_E; i++){
		for (node = 0; node < (this->element[i]).dofs_number; i++){
		initial_value((this->elements[i]).dofs[node], K, (this->element[i]).x_at_node(node));
		}
	}

}
