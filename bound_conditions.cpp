
#include "bound_conditions.hpp"
#include<cmath>
#include<vector>
#include<algorithm>
#include <limits>
#include <iostream>


namespace project{
bound_conditions::bound_conditions(const mesh& _grid)
	:m_grille(_grid){};
	 
	 //bound_conditions::~bound_conditions(){};
	 
	
	Derichtlet::Derichtlet(const mesh& m_grid,const std::vector<double>& rate)
	:
	 bound_conditions(m_grid)
	{
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double T = m_grid.Getvector_time().back();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();
		
		//From PDE object
		double spot_max = m_grid.Getvector_stock().back();
		double spot_min = m_grid.Getvector_stock()[0];
		//double f_0_T = _init_cond[0];
		//double f_N_T = _init_cond.back();
		
		matrix_derichtlet.resize(size_vec*2); 
		//std::vector<double> lower_conditions(size_vec);

		
		for(size_t i = 0; i < size_vec; ++i)
        {
			double df = exp(-rate[i]*(T-dt*i));//exp(-rate[i]*(T-dt*i));
			matrix_derichtlet[i] = m_grid.init_cond(exp(spot_min),df); // for here
			matrix_derichtlet[size_vec+i] = m_grid.init_cond(exp(spot_max),df);
        }
	};
	
	std::vector<double> Derichtlet::get_cond() const
	{
		return matrix_derichtlet;
	};
	
	Neumann::Neumann(const mesh& m_grid, const double& theta, const std::vector<double>& sigma, const std::vector<double>& rate)
	:bound_conditions(m_grid){
		//From mesh object
		double dt = m_grid.getdt();
		double dx = m_grid.getdx(); //need the stock step 
		size_t size_vec = m_grid.Getvector_time().size();
		double maturity = m_grid.Getvector_time().back();
		double S0 = m_grid.get_Spot();
		size_t size_spot = m_grid.Getvector_stock().size();	
		double S_min = m_grid.Getvector_stock()[0]; 
		double S_max = m_grid.Getvector_stock().back(); 
		
		double h =0.0001;
		double df_der = 1.0;
		
		double K1 = (m_grid.init_cond(exp(S_min +h),df_der) - m_grid.init_cond(exp(S_min -h),df_der))/(2*h);
		
		//std::cout<< K1 <<std::endl;
		
		double K3 = (m_grid.init_cond(exp(S_max +h),df_der) - m_grid.init_cond(exp(S_max -h),df_der))/(2*h);
		
		//std::cout<< K3 <<std::endl;
		
		double K2 = (m_grid.init_cond(exp(S_min+2*h),df_der) - m_grid.init_cond(exp(S_min),df_der) + m_grid.init_cond(exp(S_min-2*h),df_der))/(4*h*h);
		
		//std::cout<< K2 <<std::endl;
		
		
		double K4 = (m_grid.init_cond(exp(S_max+2*h),df_der) - m_grid.init_cond(exp(S_max),df_der) + m_grid.init_cond(exp(S_max-2*h),df_der))/(4*h*h);
		
		//std::cout<< K4 <<std::endl;
		
		//std::vector<double> coef;
		
		
		coef_neumann.push_back(K1);
		coef_neumann.push_back(K2);
		coef_neumann.push_back(K3);
		coef_neumann.push_back(K4);
		
		
		
		//print(coef);
		//print(coef_neumann);
		
		//From PDE object
		std::vector<double>  _init_cond = m_grid.get_init_vector(); //get the terminal condition vector to get f(S0,T) and f(Smax,T)
		double spot_min = m_grid.Getvector_stock()[0];
		double df = exp(-rate.back()*maturity);
		double spot_max = m_grid.Getvector_stock().back();
		double f_0_T = m_grid.init_cond(exp(spot_min),df);
		double f_N_T = m_grid.init_cond(exp(spot_max),df);
		   
		double coef_; 
		double coef_K1_K2; 
		double coef_K3_K4; 
		
		matrix_neumann.resize(size_vec);
		matrix_neumann.back() =  f_0_T*exp(-maturity*rate[0]);
		matrix_neumann.resize(size_vec*2);
		matrix_neumann.back() = f_N_T*exp(-maturity*rate.back());
		
		for (size_t it = size_vec-1; it > 0; it--)
		{
		      //reverse iterator to fill the vector from the end to the beginning
			  
			  //std::cout << it << std::endl;
			  
			  size_t tt = size_vec + it -1;
			  
			  //std::cout << tt << std::endl;
			  
			  		coef_ =(1 - dt*(1-theta)*rate[it])/(1+ dt*theta*rate[it]);
					coef_K1_K2 = -dt*((-std::pow(sigma[it],2)*K1)/2 + (std::pow(sigma[it],2)/2 - rate[it])*K2)/(1+ dt*theta*rate[it]);
					coef_K3_K4 = -dt*((-std::pow(sigma[it],2)*K3)/2 + (std::pow(sigma[it],2)/2 - rate[it])*K4)/(1+ dt*theta*rate[it]);
					
					
			  matrix_neumann[it]  = (coef_*matrix_neumann[it-1] + coef_K1_K2);
			
			  matrix_neumann[tt] = (coef_*matrix_neumann[tt-1] + coef_K3_K4);
		 }
	};
	
	
	
    std::vector<double> Neumann::get_cond() const
	{
		return matrix_neumann;
	};
	
	std::vector<double> Neumann::get_coef_neumann() const
	{
		return coef_neumann;
	};
	
	//Neumann::~Neumann(){};
			
	
}	

	
	