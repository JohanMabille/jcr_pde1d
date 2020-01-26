
#ifndef BOUND_CONDITIONS_HPP
#define BOUND_CONDITIONS_HPP
#include "mesh_spot.hpp"
#include<cmath>
#include<vector>
#include<algorithm>
#include <limits>
#include <iostream>

// we create a boundaries conditions class that is purely virtual and 2 other classes neumann and derichtlet

// boundaries class
	

namespace project{
	
class bound_conditions 
	{
		public:
		
			bound_conditions(const mesh& _grid);
			
			//~bound_conditions();
			
			//std::vector<double> cond() const;
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 //std::vector<double>  cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
		
		protected:
		
			mesh m_grille;
			//Parameters m_param;
	
		//private:
		
		//std::vector<std::vector<double>> Matrix_conditions;
			
	};
	
	
	class Neumann: public bound_conditions 
	{
		
		public:
		
		 
		 Neumann(const mesh& m_grid, const double& theta, const std::vector<double>& sigma, const std::vector<double>& rate);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;
			 std::vector<double> get_coef_neumann() const;

			std::vector<double> matrix_neumann;
			std::vector<double> coef_neumann;
			
	 //~Neumann();
			
		//std::vector<double> K_neuman; 
		
	};

	class Derichtlet: public bound_conditions 
	{
		
		public:
		
			Derichtlet(const mesh& m_grid, const std::vector<double>& rate);
			//virtual std::vector<std::vector<double>>  operator() (PDE _payoff, mesh grid, Parameters param, PayOff* option,std::vector<double>& K_neuman);
			 //std::vector<double> cond(PDE _payoff, mesh grid, Parameters param, PayOff* option);
			 std::vector<double> get_cond() const;

			std::vector<double> matrix_derichtlet;
			
		
	};
}
	
#endif 