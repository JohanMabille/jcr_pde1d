

#ifndef Resolve_HPP
#define Resolve_HPP

#include "mesh_spot.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

namespace project{	
	
	class solver
	{
	public:
		
            // Implementation: consder passing scalar types (bool, double, etc) by value
		solver(const mesh& grid,
                       const double& theta,
                       const std::vector<double>& boundaries,
                       const std::vector<std::vector<double>>& vol_mat,
                       const std::vector<std::vector<double>>& rate_mat); //constructor of the solver object
		
		std::vector<double> Mid_diag_coeff(const mesh& grid,
                                                   const bool& A,
                                                   const double& theta,
                                                   const std::vector<double>& sigma, 
                                                   const std::vector<double>& rate);

		std::vector<double> Upper_diag_coeff(const mesh& grid,
                                                     const bool& A,
                                                     const double& theta,
                                                     const std::vector<double>& sigma,
                                                     const std::vector<double>& rate);

		std::vector<double> Lower_diag_coeff(const mesh& grid,
                                                     const bool& A, 
                                                     const double& theta,
                                                     const std::vector<double>& sigma,
                                                     const std::vector<double>& rate);


		void thomas_algorithm(const std::vector<double>& upper_diag,
                                      const std::vector<double>& mid_diag,
                                      const std::vector<double>& lower_diag,
                                      const std::vector<double>& f_n1,
		                      std::vector<double>& f_sol);
		
		std::vector<std::vector<double>> get_vector_price();
		const double get_price(const size_t& n);
		
		std::vector<double> BX_vector(const std::vector<double>& upper,
                                              const std::vector<double>& mid,
                                              const std::vector<double>& low,
                                              const std::vector<double>& bound_diff,
                                              const std::vector<double>& Fn1);
		//this function will be used to compute at each time step the BX vector 
		
	private:
		
		mesh m_mesh;
		std::vector<std::vector<double>> m_results;
		

	};
}
#endif
