#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP
#include <iostream>
#include <vector>
#include <algorithm> 
#include "payoff.hpp"// Act on containers through iterators to apply modyfing/non_modifying operations
	
	
namespace project{
class VanillaOption 
	{
		 public:
		  PayOff* pay_off; //Pointer 

		  double K;
		  double r;
		  double T;
		  double sigma;

		  VanillaOption();
		  VanillaOption(double _K, double _r, double _T, 
						double _sigma, PayOff* _pay_off);
	 };
	
	double vanilla_payoff(double fwd, double strike, bool is_call);
    double bs_time_value(double fwd, double strike, double volatility, double maturity);
    double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call);

    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity);
    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call);

}
#endif
