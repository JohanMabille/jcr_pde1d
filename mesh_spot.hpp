

#ifndef mesh_spot_HPP
#define mesh_spot_HPP

#include "payoff.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

// Mesh Class	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mesh Class
	//This class builds the mesh, giving dt, dS, and the axes of the mesh as output


namespace project {

class mesh
	{
	
	public: 
	
		mesh(const double& spot, const double& maturity, const double& volatility,const long& time_step,const long& steps, PayOff* _option);
		~mesh();
		PayOff* option;
		
		std::vector<double> Getvector_time() const;
		std::vector<double> Getvector_stock() const;
		double getdx() const;
		double getdt() const;
		double get_Spot() const;
		double init_cond(const double& x,const double& df = 0) const;
		std::vector<double> get_init_vector() const;
	
	
	private: 
	
		std::vector<double> m_vector_time;
		std::vector<double> m_vector_stock;
		double m_dx;
		double m_dt;
		double m_spot;
		long m_steps;
		long m_time_step;
		std::vector<double> m_init_vector;
		
	};
	
//Method to print vector content
void print(const std::vector<double>& v);
//Method to transform the vector boundaries into matrix for resolution;
std::vector<std::vector<double>> transform_matrix(const std::vector<double>& vector_init, const double& nb_rows);
}
#endif