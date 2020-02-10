//Partie hpp	
	//Classe de base, permet aussi le cas le plus simple ou la vol et rate sont constants
	
#ifndef VOL_HPP
#define VOL_HPP
#include <iostream>
#include <vector>
#include <algorithm>

#include "mesh_spot.hpp"

namespace project {

        // Design: no virtual method, so you won't be able to replace volatility with an object from an inheriting
        // class in the client code.
	class volatility //is used to build non_constant rates also
	{
	
	public: 
		//Constructor
		volatility(const double& v,const mesh& grid); //Constructor, base class only takes initial value as input
		
		//Method to compute a vol (vol const donc les deux coeffs valent zéro)
		double compute_vol(const double& S = 0.,const double& t = 0.,const double& x = 0.,const double& y = 0.); //To compute a vol
		//Method to compute vector or vol (for one step of time)
                // Implementation: consider returning by constant reference
		std::vector<std::vector<double>> vector_vol();
		
	private:
	
		std::vector<std::vector<double>> m_vol_const;

	
	protected:
		
		// std::vector<double> m_vector_time;
		// std::vector<double> m_vector_stock;
		// std::size_t m_nb_time;
		// std::size_t m_nb_spot;	
		double m_init_vol;
		mesh m_grid;
		
	};
	 //Classe ou on va dire que la volatility est fonction du temps et du spot
	//On pourra créer ainsi d'autres classes avec différents comportement pour la vol (e.g. seulement une fonction du temps)
        // Design: entity semantic + virtual method
	class vol_surface : public volatility
	{
	public:
	
		//Niveau paramètre, la classe prends en plus les deux coeff pour définit la fonction
		vol_surface(const double& v,mesh grid, const double& coeff_tps, const double& coeff_spot);
		
		//On reprends les mêmes méthodes (du coup la les coeffs ne sont plus constant zéro comme au dessus)
		double compute_vol(const double& S,const double& t,const double& x,const double& y); 
		std::vector<std::vector<double>> vector_vol();
		
		
	private:
	
		double m_coeff_time;
		double m_coeff_spot;
		std::vector<std::vector<double>> m_vol_matrix;
		
	};
        // Design: having a compute_vol method in a class name rate_surface should ring a bell
	class rate_surface : public volatility
	{
	public:
		//Niveau paramètre, la classe prends en plus les deux coeff pour définit la fonction
		rate_surface(const double& r,mesh grid, const double& coeff_tps, const double& coeff_spot);
		
		//On reprends les mêmes méthodes (du coup la les coeffs ne sont plus constant zéro comme au dessus)
		double compute_vol(const double& S,const double& t,const double& x,const double& y); 
		std::vector<std::vector<double>> vector_vol();
		
		
	private:
	
		double m_coeff_time;
		double m_coeff_spot;
		std::vector<std::vector<double>> m_rate_matrix;
		
	};
}
#endif
