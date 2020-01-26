
#include "Resolve.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>


namespace project{

 
 solver::solver(const mesh& grid, const double& theta, const std::vector<double>& boundaries, const std::vector<std::vector<double>>& vol_mat,const std::vector<std::vector<double>>& rate_mat)
 :m_mesh(grid)
 {
	
	 std::vector<std::vector<double>> vol(vol_mat);
	 std::vector<std::vector<double>> rate(rate_mat);
	 long T = grid.Getvector_time().size(); 
	 std::size_t m = grid.Getvector_stock().size();
	 std::vector<double> init_f = grid.get_init_vector(); // get the initial X_T
	 
	 init_f.erase(init_f.begin());
	 init_f.pop_back();
	 
	 std::vector<std::vector<double>> matrix_boundaries = transform_matrix(boundaries,m-2); //get the boundaries matrix 
	 std::vector<double> cond_diff(m-2);
	 std::vector<double> up_A;
	 std::vector<double> low_A;
	 std::vector<double> mid_A;
	 std::vector<double> up_B;
	 std::vector<double> low_B;
	 std::vector<double> mid_B;
	 
	 
	 //print(init_f);
	  std::vector<double> solution_vector(init_f); //init of the solution vector 
	  m_results.insert(m_results.begin(),solution_vector);
	 
	 //initialisation de la matrice B à maturité 
	 up_B = Upper_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 low_B = Lower_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 mid_B = Mid_diag_coeff(grid,false,theta,vol.back(), rate.back());
	 
	 vol.pop_back();
	 rate.pop_back();
	 
	up_A = Upper_diag_coeff(grid,true,theta,vol.back(), rate.back());
	low_A = Lower_diag_coeff(grid,true,theta,vol.back(), rate.back());
	mid_A = Mid_diag_coeff(grid,true,theta,vol.back(), rate.back());
	
	 
	 for(size_t s = 0; s<m-2; ++s)
	 {
		 
		 if(s==0){
			 
			 double BetaN = low_A[0];
			 double BetaN1 = low_B[0];
			 
			 cond_diff[s] = BetaN1*matrix_boundaries[s].back()-BetaN*matrix_boundaries[s][T-1];
		 }
		 else if( s ==m-3){
			 
			 double AlphaN = up_A.back();
			 double AlphaN1 = up_B.back();
			 
			 cond_diff[s] = AlphaN1*matrix_boundaries[s].back()-AlphaN*matrix_boundaries[s][T-1];
			 
		 }
		 
		 
		 else{cond_diff[s] =0;};
	 }; //initialisation du vecteur de C(n+1) - C(n) 
	 
	
	std::vector<double> B = BX_vector(up_B,mid_B,low_B,cond_diff,init_f);
	
	thomas_algorithm(up_A, mid_A, low_A, B, solution_vector);
	 
	
	 for(int i = grid.Getvector_time().size() - 2; i != 0; --i)
	 {

		
		m_results.insert(m_results.begin(),solution_vector);
		
		up_B = Upper_diag_coeff(grid, false,theta,vol.back(), rate.back());
	    low_B = Lower_diag_coeff(grid, false,theta,vol.back(), rate.back());
		mid_B = Mid_diag_coeff(grid, false,theta,vol.back(), rate.back());
		
		rate.pop_back();
		vol.pop_back();
		
		up_A = Upper_diag_coeff(grid,true,theta,vol.back(), rate.back());
		low_A = Lower_diag_coeff(grid,true,theta,vol.back(), rate.back());
		mid_A = Mid_diag_coeff(grid,true,theta,vol.back(), rate.back());
		
		 //create a matrix containing all the prices computed by the solver in order to be displayed afterwards
		
		
		for(size_t s = 0; s<m-2; ++s){
		 	
		if(s==0){
			 
			 double BetaN = low_A[0];

			 double BetaN1 = low_B[0];
			 
			 cond_diff[s] = BetaN1*matrix_boundaries[s][i]-BetaN*matrix_boundaries[s][i-1];
		 }
		 else if( s==m-3){
			
			 double AlphaN = up_A.back();
			 double AlphaN1 = up_B.back();
			 
			 cond_diff[s] = AlphaN1*matrix_boundaries[s][i]-AlphaN*matrix_boundaries[s][i-1];
			 
		 }
		 else{ cond_diff[s] =0.0;};
	     }; 
		 
		 //print(cond_diff);
		//std::cout << cond_diff.size() << std::endl;
		
		//print(solution_vector);
		
	    B = BX_vector(up_B,mid_B,low_B,cond_diff,solution_vector);
		
		thomas_algorithm(up_A, mid_A, low_A, B, solution_vector);
		 
	 }; 

 }; 
 
std::vector<std::vector<double>> solver::get_vector_price(){
	
	return m_results;
};

const double solver::get_price(const double& n){
	
	return m_results[0][n];
	
	//return the price of the option at time 0 and for a certain point of the mesh;
	
};
 
 std::vector<double> solver::BX_vector(const std::vector<double>& upper, const std::vector<double>& mid, const std::vector<double>& low,const std::vector<double>& bound_diff,const std::vector<double>& Fn1){
//this procedure creates the right_hand vector of the AX = D equation based on BX(n+1) + C(n+1) - C(n)
	std::size_t N = mid.size(); // as the resolution si between f1 and fN-1
	
	std::vector<double> BX;
	
	
	BX.push_back(mid[0]*Fn1[0] + upper[0]*Fn1[1] +bound_diff[0]);
	
	BX[0] = mid[0]*Fn1[0] + upper[0]*Fn1[1] +bound_diff[0];
	
	for(size_t i = 1; i < N-1; i++){

   	BX.push_back(mid[i]*Fn1[i] + upper[i]*Fn1[i+1] +low[i-1]*Fn1[i-1] + bound_diff[i]);
	
	
	};
	
	
	BX.push_back(mid.back()*Fn1.back() + low.back()*Fn1[N-2] + bound_diff.back());
	
	return BX;
};
	
//set of function to define the A matrix (3 better than just one huge matrix ?) 
std::vector<double> solver::Mid_diag_coeff(const mesh& grid, const bool& A,const double& theta, const std::vector<double>& sigma, const std::vector<double>& rate){
	
	std::vector<double> vol(sigma);
	std::vector<double> tx(rate);
	vol.pop_back();
	vol.erase(vol.begin());
	
	tx.pop_back();
	tx.erase(tx.begin());
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	double size_gamma = grid.Getvector_stock().size()-2; //no minus-1 as we are on diagonal 
	
	
	
		//create the vector that holds the diagonal 
	std::vector<double> gamma_coefficient(size_gamma);
	
	for (std::size_t i = 0; i < size_gamma; ++i){
		
		gamma_coefficient[i] =  std::pow(vol[i],2)/std::pow(dx,2) + tx[i];
		
		if (A==false){
			
			gamma_coefficient[i] = -dt*(1-theta)*gamma_coefficient[i] + 1.0;
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			gamma_coefficient[i] = dt*theta*gamma_coefficient[i] + 1.0;
		};
	};
	
	return gamma_coefficient;
	
}; 


std::vector<double> solver::Upper_diag_coeff(const mesh& grid,const bool& A,const double& theta,const std::vector<double>& sigma, const std::vector<double>& rate){
	
	std::vector<double> vol(sigma);
	std::vector<double> tx(rate);
	vol.pop_back();
	vol.pop_back();
	vol.erase(vol.begin());
	
	tx.pop_back();
	tx.pop_back();
	tx.erase(tx.begin());
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	double size_alpha = grid.Getvector_stock().size()-3; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> alpha_coefficient(size_alpha);
	
	for (std::size_t i = 0; i < size_alpha; ++i){
		
		alpha_coefficient[i] =  (-std::pow(vol[i],2))/(2*std::pow(dx,2)) + (std::pow(vol[i],2))/(4*dx) - (tx[i])/(2*dx);
		
		 
		if (A==false){
			
			alpha_coefficient[i] = -dt*(1-theta)*alpha_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			alpha_coefficient[i] = dt*theta*alpha_coefficient[i];
		};
	};
	
	return alpha_coefficient;
	
};


std::vector<double> solver::Lower_diag_coeff(const mesh& grid, const bool& A,const double& theta, const std::vector<double>& sigma, const std::vector<double>& rate){
	
	std::vector<double> vol(sigma);
	std::vector<double> tx(rate);
	vol.pop_back();
	vol.erase(vol.begin());
	vol.erase(vol.begin());
	
	tx.pop_back();
	tx.erase(tx.begin());
	tx.erase(tx.begin());
	double dt = grid.getdt(); //need the time step 
	double dx = grid.getdx(); //need the stock step 
	double size_beta = grid.Getvector_stock().size()-3; //minus 1 because we are on the uppdiag 
	
	//create the vector that holds the diagonal 
	std::vector<double> beta_coefficient(size_beta);
	
	for (std::size_t i = 0; i < size_beta; ++i){
		
		beta_coefficient[i] =  (-std::pow(vol[i],2))/(2*std::pow(dx,2)) + std::pow(vol[i],2)/(4*dx) + (tx[i])/(2*dx);
		
		
		if (A==false){
			
			beta_coefficient[i] = -dt*(1-theta)*beta_coefficient[i];
			//this condition checks if we are computing the A matrix or the B matrix as B is only the opposite of A 
		}
		
		else{
			
			beta_coefficient[i] = dt*theta*beta_coefficient[i];
		};
	};
	
	return beta_coefficient;
	
};

//this is the thomas algo for inverting the matrix  
void solver::thomas_algorithm(const std::vector<double>& upper_diag, const std::vector<double>& mid_diag, const std::vector<double>& lower_diag, const std::vector<double>& f_n1, std::vector<double>& x) {
  size_t nb_spot = f_n1.size();
  
  std::vector<double>  a(lower_diag);
  std::vector<double>  c(upper_diag);
  std::vector<double>  b(mid_diag);
  x = f_n1;
  
  c.push_back(0.0);
  a.insert(a.begin(),0.0);
  // Create the temprary vectors to store new coef                                                                                                                                                                                                                                                                                                                                                

//std::cout << "vecteur init" << f_n1.size() << std::endl;
//print(f_n1);
//Step 1  
  c[0] = c[0] / b[0];
  x[0] = x[0] / b[0];

//Steps 2
  //forward sweep                                                                                                                                                  
  for (std::size_t i=1; i<nb_spot; i++) 
  {
    const double m = 1.0 / (b[i] - a[i] * c[i-1]);
    c[i] = c[i] * m;
    x[i] = (x[i] - a[i] * x[i-1]) * m;
  }
  
 //Step 3

  //reverse sweep, used to update the solution vector f                                                                                                                                                 
  for (std::size_t i=nb_spot-1; i-- > 0; ) 
  {
    x[i] -= c[i] * x[i+1];
  }
};

solver::~solver(){}; 

}