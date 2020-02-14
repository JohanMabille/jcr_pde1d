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

    // Implementation: consider passing sol by const ref to avoid many copies
    Greeks(const mesh& grid, solver sol, const std::vector<double>& init_values ={0.0,0.0,0.0,0.0});
    
    ~Greeks();
	
    // Implementation: returns by const reference
    // Implementation: these methods should be const
	std::vector<double> get_delta();
	std::vector<double> get_gamma();
	std::vector<double> get_theta();
	
	private:
	
	std::vector<double> m_delta;
	std::vector<double> m_gamma;
	std::vector<double> m_theta;
    
    
};

}


#endif /* Greeks_hpp */
