
#include "closed_form.hpp"
#include "Resolve.hpp"
#include "mesh_spot.hpp"
#include "vol.hpp"
#include "bound_conditions.hpp"
#include "payoff.hpp"
#include "Greeks.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>

////
int main(int argc, char* argv[])
{
	// Create the option parameters
	double S = 100;
	double K = 100;  // Strike price
	double r = 0.05;   // Risk-fr ee rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry
	double theta_ = 0.5;
	
	// mesh discretisation parameters
	long nb_step_spot =300;    // Spot goes from [0.0, 1.0]
	long nb_step_time =200; 
	
	
	// Create the PayOff object (strike as parameter)
	project::PayOff* option = new project::PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	project::mesh grille(S,T,v,nb_step_time,nb_step_spot,option);
	
	//Create PDE object (to solve the problem, we provide to the class, the grid previoulsy defined in the class PDE) 
	//PDE m(option,grille.getdx(),grille.getdt(),grille.Getvector_time(),grille.Getvector_stock());
	
	//test of the boundaries : 
	long s = grille.Getvector_stock().size();
	long _t_ = grille.Getvector_time().size();
	std::vector<double> sigma(s,v);
	std::vector<double> rate(s,r);
	project::Derichtlet c(grille, rate);
	project::Neumann c2(grille, theta_, sigma, rate);
	
	double tm = 0.1;
	double sco = 0.5;
	
	//project::vol_surface vo(v, grille, tm, sco);
	
	std::vector<std::vector<double>> test = project::transform_matrix(c.get_cond(),s);
	
	std::vector<std::vector<double>> vol_mat;
	std::vector<std::vector<double>> rate_mat;
	
	 
	for(long i=0; i<_t_;i++){
		
		vol_mat.push_back(sigma);
		rate_mat.push_back(rate);
	};
	
	//std::cout << vol_mat.size() << std::endl;
	//std::cout << vol_mat[0].size() << std::endl;
	
	std::vector<double> init_f(grille.get_init_vector());
	
	std::vector<std::vector<double>> res;
	
	//project::solver sol(grille, res);
		 
	std::vector<double> cond_test;
	
	for(int i =0; i <test.size();i++){
		
		
		cond_test.push_back(test[i].back());
	}

	project::solver sol(grille, theta_,c.get_cond(), vol_mat, rate_mat);
	
	//std::vector<double> ccc = c2.get_coef_neumann();
	
	std::vector<std::vector<double>> price = sol.get_vector_price();
	
/* 	for(int i=0; i<2;i++){
		
		std::cout << "price at time " << i << std::endl;
		project::print(price[i]);
	} */
	
	
	for(int i=1; i<grille.Getvector_stock().size()-1;i++){
		
		 if (grille.Getvector_stock()[i] == log(S)){
			 
			 double indice = i;
			 
			 double p = sol.get_price(i);
			 
			 double bs_ = project::bs_price(S, K, v, T, true);
			 
			 std::cout << "price is " << p  << " and index in vector is "<< indice << std::endl;
			 std::cout << "price BS is " << bs_ << std::endl;
		
		
	};
	
	}

	
	project::Greeks g(grille, sol);
	
	std::vector<double> delta = g.get_delta();
	std::vector<double> gamma = g.get_gamma();
	std::vector<double> theta = g.get_theta();
	
/* 	project::print(delta);
	project::print(gamma);
	project::print(theta);
		
	std::cout<< "fonction calcul cond init:" << std::endl;
	std::cout << grille.init_cond(S) << std::endl;

	std::cout<< "Vecteur de prix (log):" << std::endl;
	project::print(grille.Getvector_stock());
	std::cout<< "Size vecteur stock:" << std::endl;
	std::cout<< grille.Getvector_stock().size() <<std::endl;
	
	std::cout<< "Vecteur de temps:" << std::endl;
	project::print(grille.Getvector_time());
	std::cout<< "Size vecteur temps:" << std::endl;
	std::cout<< grille.Getvector_time().size() <<std::endl;
	
	std::cout<< "Vecteur cond init (not log):" << std::endl;
	project::print(grille.get_init_vector());
	std::cout<< "Size vecteur cond ini:" << std::endl;
	std::cout<< grille.get_init_vector().size() <<std::endl;
	
	std::cout<< "dx:" << std::endl;
	std::cout << grille.getdx() << std::endl;
	
	std::cout<< "dt:" << std::endl;
	std::cout << grille.getdt() << std::endl;
	
	std::cout<< "vecteur dirichlet:" << std::endl;
	project::print(c2.get_cond());
	std::cout<< "Size vecteur dirichlet:" << std::endl;
	std::cout<< c2.get_cond().size() <<std::endl;
	
	
	std::cout<< "vecteur Neumann:" << std::endl;
	project::print(c2.get_cond()); */

	delete option;



	return 0;
	
}
