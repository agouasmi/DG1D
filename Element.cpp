#include "Element.h"

template <int K, int order>
double Element<K,order>::x_at_node(unsigned int node){
	
	return center + size/2.*xi_at_node(node);

}	

template <int K, int order>
double Element<K, order>::xi_at_node(unsigned int node){

	return quad_points[node]; // here the degrees of freedom are exactly the Lobatto points
}

template <int K, int order>
double Element<K,order>::basis_function(unsigned int node, double xi){
	/* Lagrange basis functions */	

	int i; double out = 1.;
	double xi_ref = xi_at_node(node);
	double xi_temp;

	for (i = 0; (i < dofs_number) && (i != node); i++ ){
		xi_temp = xi_at_node(i);
		out *= (xi - xi_temp) / (xi_ref - xi_temp);
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
        for (i = 0; (i < dofs_number) && (i != node); i++ ){
		xi_temp = xi_at_node(i);
		denom *= 1 / (xi_ref - xi_temp);
	}

	// Do the rest	
        for (j = 0; (j < dofs_number) && (j != node); j++){
		
		temp = 1.;
		for (i = 0; (i < dofs_number) && (i != node) && (i != j); i++ ){
			xi_temp = xi_at_node(i);
			temp *= (xi - xi_temp);
		}
		out += temp;
	}
	
	return out*denom;
}

