#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <vector>
#include "Greeks.hpp"


namespace project {
Greeks::Greeks(const mesh& grid, solver sol, const std::vector<double>& init_values){
	
//price matrix is T row and nb_spots -2 columns;
double dt = grid.getdt(); //get the time step
double dx = grid.getdx(); //get the stock step
long T = grid.Getvector_time().size(); //get the number of time steps
double nb_spot_values = grid.Getvector_stock().size(); //get the number of spot values
std::vector<std::vector<double>> price_matrix = sol.get_vector_price(); //get the matrix of option prices
std::vector<double> stock = grid.Getvector_stock();
//price_matrix.pop_back(); //get rid of the last row that is the price at T
 
m_delta.resize(nb_spot_values -2,0.0);
m_gamma.resize(nb_spot_values -2,0.0);
m_theta.resize(nb_spot_values -2,0.0);

// Boundaries conditions
// Assuming that if the user imputs a Dirichlet condition, then the Neumann conditions (K1...) are null
m_delta[0] = init_values[0]/252.;
m_gamma[0] = init_values[1];
m_theta[0] = (price_matrix[1][0] - price_matrix[0][0])/(365.*dt);

//std::cout << m_delta[0] << std::endl;
//std::cout << m_gamma[0] << std::endl;
//std::cout << m_theta[0] << std::endl;

m_delta.back() = init_values[2]/252.;
m_gamma.back() = init_values.back();
m_theta.back() = (price_matrix[1].back() - price_matrix[0].back())/(365.*dt);
    
//std::cout << m_delta.back() << std::endl;
//std::cout << m_gamma.back() << std::endl;
//std::cout << m_theta.back() << std::endl;
  
 for(size_t i = 1; i < m_delta.size()-1; ++i){
 
	 double ds_left = exp(stock[i]) - exp(stock[i] - i*dx);
	 double ds_right = exp(stock[i] + i*dx) - exp(stock[i]);
	 double delta_right = (price_matrix[0][i+1] - price_matrix[0][i])/(252.*ds_right);
	 double delta_left = (price_matrix[0][i] - price_matrix[0][i-1])/(252.*ds_left);
     m_delta[i] = (delta_left + delta_right)/2.;
     m_gamma[i] = (m_delta[i+1] - m_delta[i-1])/(ds_left*ds_right);
     m_theta[i] = (price_matrix[1][i] - price_matrix[0][i])/365.*dt;
     
 };
 
};

std::vector<double> Greeks::get_delta(){
    
    return m_delta;
};

std::vector<double> Greeks::get_gamma(){
    
    return m_gamma;
};

std::vector<double> Greeks::get_theta(){
    
    return m_theta;
};

Greeks::~Greeks(){};

}