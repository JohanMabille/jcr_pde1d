#include "vol.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

//Classe de base
	//Constructor
	
namespace project {
	//Constructor
	volatility::volatility(const double& v,const mesh& grid)
	:
	 m_init_vol(v),
	 m_grid(grid)
	{
		//m_vector_time = m_grid.Getvector_time();
		//m_vector_stock = m_grid.Getvector_stock();
		//m_nb_time = m_vector_time.size();
		//m_nb_spot = m_vector_stock.size();	
	};
	
	//Here no computation since vol is constant --> just return init vol 
	double volatility::compute_vol(const double& S,const double& t,const double& x,const double& y)
	{
		return m_init_vol;
	};
	
	std::vector<std::vector<double>> volatility::vector_vol()
	{
		std::size_t nb_step = m_grid.Getvector_stock().size();
		std::size_t nb_step_time =  m_grid.Getvector_time().size();
		std::vector<double> volat_vector;
		
		for(std::size_t j = 0; j < nb_step_time ; ++j)
		{
			m_vol_const.push_back(volat_vector);
			for (std::size_t i = 0; i < nb_step ; ++i)
			{
				m_vol_const[j].push_back(compute_vol());
			}
		}
		return m_vol_const;
	};
//Classe qui hérite pour la surface de vol

	vol_surface::vol_surface(const double& v, mesh grid,const double& coeff_tps, const double& coeff_spot)
	:
	 volatility(v, grid), //Normalament on obtient ainsi les variables protected de volatility
	 m_coeff_spot(coeff_spot),
	 m_coeff_time(coeff_tps)
	{
	};
	
	//Here computation of a vol thanks to parameters
	double vol_surface::compute_vol(const double& S,const double& t,const double& x,const double& y)
	{
		double vol = x*t + y*S;
		
		return vol;
	};
	
	std::vector<std::vector<double>> vol_surface::vector_vol()
	{
		std::size_t nb_step = m_grid.Getvector_stock().size();
		std::size_t nb_step_time =  m_grid.Getvector_time().size();
		std::vector<double> volat_vector;
		
		//m_vol_const.resize(nb_step,0.);
		for(std::size_t j = 0; j < nb_step_time ; ++j)
		{
			m_vol_matrix.push_back(volat_vector);
			
			for (std::size_t i = 0; i < nb_step ; ++i)
			{
				m_vol_matrix[j].push_back(compute_vol(m_grid.Getvector_stock()[i],m_grid.Getvector_time()[j],m_coeff_time,m_coeff_spot));
			}
		}		
		return m_vol_matrix;
	};	
rate_surface::rate_surface(const double& r, mesh grid,const double& coeff_tps, const double& coeff_spot)
	:
	 volatility(r, grid), //Normalament on obtient ainsi les variables protected de volatility
	 m_coeff_spot(coeff_spot),
	 m_coeff_time(coeff_tps)
	{
	};
	
	//Here computation of a vol thanks to parameters
	double rate_surface::compute_vol(const double& S,const double& t,const double& x,const double& y)
	{
		double rate = x*t + y*S;
		
		return rate;
	};
	
	std::vector<std::vector<double>> rate_surface::vector_vol()
	{
		std::size_t nb_step = m_grid.Getvector_stock().size();
		std::size_t nb_step_time =  m_grid.Getvector_time().size();
		std::vector<double> rate_vector;
		
		//m_vol_const.resize(nb_step,0.);
		for(std::size_t j = 0; j < nb_step_time ; ++j)
		{
			m_rate_matrix.push_back(rate_vector);
			
			for (std::size_t i = 0; i < nb_step ; ++i)
			{
				m_rate_matrix[j].push_back(compute_vol(m_grid.Getvector_stock()[i],m_grid.Getvector_time()[j],m_coeff_time,m_coeff_spot));
			}
		}		
		return m_rate_matrix;
	};	

}