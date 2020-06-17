/////////////////////////////////////
// Author: Rei Yatsuhashi
// Institution: American School In Japan
// Date Created: 2020/06/17
// Date Modified: 2020/06/17
// Purpose: Class for Constraints in Redistricting Adjacency
/////////////////////////////////////

// Header files
#include <RcppArmadillo.h>
#include <redist_aList.cpp>

using namespace Rcpp;

// Class to consolidate methods relating to the constraints and tempering for redistricting

class redist_aList_beta: public redist_aList {
  
  private:
    
    double pct_dist_parity;
  
	  NumericVector beta_sequence;
  
	  NumericVector beta_weights;
  
	  NumericMatrix ssdmat;
    
    NumericVector betas = {0.0, 0.0, 0.0, 0.0};
    
    IntegerVector anneals = {0, 0, 0, 0};
  
    List constraint_vals;

    /* Inputs to function:
     
     pct_dist_parity: strength of population parity requirement
     
     beta_sequence: sequence of betas to anneal over
     
     beta_weights: prior weights on the beta sequence
     
     ssdmat: matrix of squared distances between geographic units. For constraining
     on compactness
     
     betas: {beta_population, beta_compact, beta_segregation, beta_similar}
     
       beta_population: strength of constraint for achieving population parity
     
       beta_compact: strength of constraint for achieving district compactness
     
       beta_segregation: strength of constraint for packing group into district
     
       beta_similar: strength of constraint for examining plans similar to original district
     
     anneals: {anneal_beta_population, anneal_beta_compact, anneal_beta_segregation, anneal_beta_similar}
     
       anneal_beta_population: flag for whether to anneal the beta pop parameter
     
       anneal_beta_compact: flag for whether to anneal the beta compactness parameter
     
       anneal_beta_segregation: flag for whether to anneal the beta segregation parameter
     
       anneal_beta_similar: flag for whether to anneal the beta similarity parameter
     
     */
  
  public:
  
    // Constructor for constraint-related values
    void init_constraints(double p, NumericVector b_s, NumericVector b_w, NumericMatrix ssd);
    void init_betavals(NumericVector b);
    void init_annealvals(IntegerVector a);
    
    // Modifiers for constraint-related values
  
    // Function to calculate the strength of the beta constraint for population
    List calc_betapop(arma::vec current_dists, arma::vec new_dists,
		  NumericVector pops,
		  double beta_population,
		  NumericVector distswitch)
  

}

// Constructor
void redist_aList_beta::init_constraints(double p, NumericVector b_s, NumericVector b_w, NumericMatrix ssd) 
{
  
  pct_dist_parity = p;
  beta_sequence = b_s;
  beta_weights = b_w;
  ssdmat = ssd;
  
}

void redist_aList_beta::init_betavals(NumericVector b)
{
  
  betas = b;
  
}

void redist_aList_beta::init_annealvals(IntegerVector a)
{
  
  anneals = a;
  
}







