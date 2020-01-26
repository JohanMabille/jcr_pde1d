
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
	// Spot goes from [0.0, 1.0]
	long nb_step_spot =258;    
	long nb_step_time =252; 
	
	
	// Create the PayOff object (strike as parameter)
	project::PayOff* option = new project::PayOffCall(K);
	// Create the mesh object (parameters (S,Vol, Nb time step, nb stock step)
	project::mesh grille(S,T,v,nb_step_time,nb_step_spot,option);
	

	//Boundaries : 
	std::size_t s = grille.Getvector_stock().size();
	std::size_t  _t_ = grille.Getvector_time().size();
	std::vector<double> sigma(s,v);
	std::vector<double> rate(s,r);
	project::Derichtlet c(grille, rate);
	project::Neumann c2(grille, theta_, sigma, rate);
	
	double coef_spot_v = 0.1;
	double coef_time_v = 0.01;
	
	//Create the non constante volatility object and rate object
	
	project::volatility vol_obj(v,grille);
	
	project::vol_surface vol_obj_2(v,grille,coef_spot_v,coef_time_v);
	
	project::volatility vol_obj_r(r,grille);
	
	project::vol_surface vol_obj_r_2(v,grille,coef_spot_v,coef_time_v);
	

	//Create the constante volatility and rate 
	std::vector<std::vector<double>> vol_mat;
	std::vector<std::vector<double>> rate_mat;
	 
	for(long i=0; i<_t_;i++){
		
		vol_mat.push_back(sigma);
		rate_mat.push_back(rate);
	};
	
	//create the initiale vector T
	std::vector<double> init_f(grille.get_init_vector());
	
	std::vector<std::vector<double>> res;
	
	//project::solver sol(grille, res);
	// Uncomment this part to run the solver with non constant vol and non constant rate 
	//project::solver sol_2(grille, theta_,c2.get_cond(), vol_obj_2.vector_vol(), vol_obj_r_2.vector_vol());
	// std::vector<std::vector<double>> price2 = sol2.get_vector_price();
	
	//solver with constant rate and constant vol
	project::solver sol(grille, theta_,c2.get_cond(), vol_mat, rate_mat);


	std::vector<std::vector<double>> price = sol.get_vector_price();

// Uncomment this to print all the vector of price through time 
/* 	for(int i=0; i<price.size();i++){
		
		std::cout << "price at time " << i << std::endl;
		project::print(price[i]);
	} */
	
	//this print both the closed formula price and the FDM one 
	for(int i=1; i<grille.Getvector_stock().size()-1;i++)
	{
		
		 if (grille.Getvector_stock()[i] == log(S))
		{
			 
			 size_t indice = i;
			 
			 double p = sol.get_price(i);
			 
			 double bs_ = project::bs_price(S, K, v, T, true);
			 
			 std::cout << "price is " << p  << " and index in vector is "<< indice << std::endl;
			 std::cout << "price BS is " << bs_ << std::endl;
		
		
		};
	
	}
	
	//Nous avons essayé différentes choses pour débuger nos prix, corriger des erreurs de coefficients, rendre plus intelligible l'algorithme de Thomas
	//Nous avons printer les tailles et les valeurs des différents vecteurs utilisés et aucun problème ne semble venir de ceci
	//Nous sommes conscients que le prix devrait être bon avec 1000 steps pour le spot néanmoins nous trouvant le bon résultat à 10e-5 lorsque nous utilisons 258 points de spot 
	//Ce résultat est valide pour r = 0.05, v= 0.2 et vol/rate constants et nous trouvons le même résultat avec Neumann et Derichtlet 
	//Lorsque nous prenons des matrices de vol et de taux non constants (fonction du spot et du temps) nous trouvons un résultat relativement prochain avec 1000 points de spot et 252 de temps
	//Néanmoins il semble il y avoir un problème plus important comme l'indique les valeurs assez folles prisent par les grecques.

	//this creates the Greek object from the solver object 
	project::Greeks g(grille, sol);
	
	std::vector<double> delta = g.get_delta();
	std::vector<double> gamma = g.get_gamma();
	std::vector<double> theta = g.get_theta();
	
	std::cout << "Delta " << std::endl;
	project::print(delta);
	std::cout << "Gamma " << std::endl;
	project::print(gamma);
	std::cout << "Theta " << std::endl;
	project::print(theta);
	
	delete option;
		
/* 	std::cout<< "fonction calcul cond init:" << std::endl;
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
	project::print(c.get_cond());
	std::cout<< "Size vecteur dirichlet:" << std::endl;
	std::cout<< c.get_cond().size() <<std::endl;
	
	std::cout<< "vecteur Neumann:" << std::endl;
	project::print(c2.get_cond()); */

	



	return 0;
	
}
