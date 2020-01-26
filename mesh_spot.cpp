
#include "mesh_spot.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

namespace project{
	mesh::mesh(const double& spot, const double& maturity,const double& volatility, const long& time_step, const long& steps,PayOff* _option)
	:
	 option(_option),
	 m_dt(maturity/time_step),
	 m_spot(spot),
	 m_steps(steps),
	 m_time_step(time_step)
	{
		//The mesh is centered around log(S0)
		double S0 = log(m_spot);
		double high_bound = S0 + 5*volatility*sqrt(maturity);
		double low_bound = S0 - 5*volatility*sqrt(maturity);
		
		//dS
		m_dx = (high_bound - low_bound)/steps;
		
		//Axe of spot price
		for (long i = 0; i < m_steps+1 ; ++i)
		{
			
			m_vector_stock.push_back(high_bound + (i - m_steps)*m_dx);
		}
		//Axe of time
		for (std::size_t j = 0; j < m_time_step+1 ; ++j)
		{
			m_vector_time.push_back(j*m_dt);	
		}
		
		size_t nb_step = m_vector_stock.size();
		//std::vector<double> m_init_vector(nb_step);
		
		for (std::size_t i = 0; i < nb_step ; ++i)
		{
			m_init_vector.push_back(init_cond(exp(m_vector_stock[i])));
		}
	};

		
	//Multiple Get methods to return private variables of the class


	//useful to get the vector of time from the mesh 
	std::vector<double> mesh::Getvector_time()const
	{
		return m_vector_time;
	}; 
	//useful to get the vector of stock path from the mesh 
	std::vector<double> mesh::Getvector_stock() const
	{
		return m_vector_stock;
	}; 
	double mesh::init_cond(const double& x,const double& df) const 
	{
	  return option->operator()(x,df);
	};
	
	//we will need to get the dx for CN algo	
	double mesh::getdx() const
	{
		return m_dx;
	}; 
	std::vector<double> mesh::get_init_vector() const //const forbid to modify the state of my object
    {
        return m_init_vector;
    };
	
	// we will need the dt for CN algo 
	double mesh::getdt() const
	{
		return m_dt;
	}; 
	
	double mesh::get_Spot() const
	{
		return m_spot;
	};
	
	//Destructor
	mesh::~mesh() {};
    void print(const std::vector<double>& v)
    {
        for(size_t i = 0; i < v.size(); ++i)
        {
            std::cout << v[i] << ","; 
        }
        std::cout << std::endl;
    };

	std::vector<std::vector<double>> transform_matrix(const std::vector<double>& vector_init, const double& nb_rows){
	
	double endv; 
	//std::vector<double> upper_bound(vector_init);
	//std::vector<double> lower_bound(vector_init);
	
	
	 
	 if (vector_init.size()%2==0){endv=vector_init.size()/2;}
	 else {endv=std::floor(vector_init.size()/2) -1;}; //case useless !
	 
	 std::vector<double> upper_bound(endv);
	 std::vector<double> lower_bound(endv);
	 
	 std::copy(vector_init.begin(), vector_init.begin() + endv , upper_bound.begin());
	 
	 std::copy(vector_init.begin() + endv, vector_init.end(), lower_bound.begin());
	 
	 
	std::vector<double> row_0(upper_bound.size(), 0.0);
	 
	std::vector<std::vector<double>> matrix;
	
	matrix.resize(nb_rows, std::vector<double>(upper_bound.size()));
	
	matrix.front() = upper_bound;
	
	//std::cout << "row 0 size " << row_0.size() <<  "up size " << upper_bound.size() <<  "low size " << lower_bound.size() << std::endl;
	
	for(int i=1; i<nb_rows-1; i++){
		
		
		matrix[i] = row_0;
	};
	
	matrix.back() = lower_bound;
	 
	return matrix;		
	};
}
	
