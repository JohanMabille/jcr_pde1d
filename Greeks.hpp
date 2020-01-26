#ifndef Greeks_hpp
#define Greeks_hpp
#include <iostream>
#include <vector>
#include <algorithm>
#include "mesh_spot.hpp"
#include "bound_conditions.hpp"
#include "Resolve.hpp"

namespace project{
class Greeks
{
    
public:

    Greeks(const mesh& grid, solver sol, const std::vector<double>& init_values ={0.0,0.0,0.0,0.0});
    
    ~Greeks();
	
	std::vector<double> Greeks::get_delta();
	std::vector<double> Greeks::get_gamma();
	std::vector<double> Greeks::get_theta();
	
	private:
	
	std::vector<double> m_delta;
	std::vector<double> m_gamma;
	std::vector<double> m_theta;
    
    
};

}


#endif /* Greeks_hpp */
