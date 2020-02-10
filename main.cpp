
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
	double r = 0.0;   // Risk-fr ee rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.0/365.;    // One year until expiry
	double theta_ = 0.5;
	
	// mesh discretisation parameters
	// Spot goes from [0.0, 1.0]
	long nb_step_spot =500;    
	long nb_step_time =1; 
	
	
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
	
	// Uncomment this part to run the solver with non constant vol and non constant rate 
	//project::solver sol_2(grille, theta_,c2.get_cond(), vol_obj_2.vector_vol(), vol_obj_r_2.vector_vol());
	// std::vector<std::vector<double>> price2 = sol2.get_vector_price();
	
	// Uncomment this part to run the solver with Derichtlet conditions 
	//project::solver sol_d(grille, theta_,c.get_cond(), vol_mat, rate_mat);
	// std::vector<std::vector<double>> price_d = sol_d.get_vector_price();
	
	//solver with constant rate and constant vol
	//project::solver sol(grille, theta_,c2.get_cond(), vol_mat, rate_mat);
	project::solver sol(grille, theta_,c.get_cond(), vol_mat, rate_mat);


	std::vector<std::vector<double>> price = sol.get_vector_price();

// Uncomment this to print all the vector of price through time 
/* 	for(int i=0; i<price.size();i++){
		
		std::cout << "price at time " << i << std::endl;
		project::print(price[i]);
	} */
	
	//this print both the closed formula price and the FDM one 
        std::vector<double> spots = grille.Getvector_stock();
        for(size_t i = 1; i < spots.size() - 1; ++i)
        {
            double s = std::exp(spots[i]);
            double df = std::exp(-r * T);
            
            // As stated on the gitter channel, the name is misleading. This
            // is actually the Black formula (on forward), not the BS formula
            // (on spot)
            //double bs = project::bs_price(s, K, v, T, true);
            double bs = project::bs_price(s / df, K, v, T, true) * df;

            std::cout << "i = " << i
                      << ", S = " << s
                      << ", BS = " << bs
                      << ", price = " << sol.get_price(i)
                      << ", diff = " << sol.get_price(i) - bs
                      << std::endl;
        }
        /*for(int i=1; i<grille.Getvector_stock().size()-1;i++)
	{

		
		 if (grille.Getvector_stock()[i] == log(S))
		{
			 
			 size_t indice = i;
			 
			 double p = sol.get_price(i);
			 
			 double bs_ = project::bs_price(S, K, v, T, true);
			 
			 std::cout << "price is " << p  << " and index in vector is "<< indice << std::endl;
			 std::cout << "price BS is " << bs_ << std::endl;
		
		
		};
	
	}*/
	
	//this creates the Greek object from the solver object 
	
	//project::Greeks g2(grille, sol); uncomment this to have 
	
	/*std::vector<double> coef = c2.get_coef_neumann();
	project::Greeks g(grille, sol, coef);
	
	std::vector<double> delta = g.get_delta();
	std::vector<double> gamma = g.get_gamma();
	std::vector<double> theta = g.get_theta();
	
	std::cout << "Delta " << std::endl;
	project::print(delta);
	std::cout << "Gamma " << std::endl;
	project::print(gamma);
	std::cout << "Theta " << std::endl;
	project::print(theta);*/
	
	//Le prix converge seulement pour 258 spot et 252 steps de temps pour les paramètres actuellement retenus (vol de 20%, rate de 5% et theta 1/2)
	//La convergence actuelle est achevée uniquement avec Neumann (différence d'environ 10 centimes avec Derichtlet)
	//APrès avoir vérifié à nouveau les coefficients ainsi que tous les outputs intermédiaires, il semble certain que nous passons à côté de quelque chose qui nous empêche de converger
	//Mais nous ne sommes pas certains de savoir quoi ni comment, pour débeuger nous avons printer différents outputs intermédiaires (conditions initiales, vector solution avant/après le passage dans l'inversion de matrice)
	//Mais rien ne nous semble péché de ce côté là...
	//L'absence de convergence totale et sans accroc en prix nous empêche aussi d'avoir des grecques corrects comme vous pourrez le constater.
	//Nous restons à votre disposition si vous aviez le temps après correction de nous expliquer ce que nous avons rater
	
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
