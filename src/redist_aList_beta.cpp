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
    
    List betas = List::create(_["population"] = 0.0, _["compact"] = 0.0, _["segregation"] = 0.0, _["similar"] = 0.0);
    
    List anneals = List::create(_["population"] = 0, _["compact"] = 0, _["segregation"] = 0, _["similar"] = 0);
  
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
    List calc_betapop(arma::vec new_dists)
      
    // 

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

// Function to calculate the strength of the beta constraint for population
List calc_betapop(arma::vec new_dists)
{

  /* Inputs to function
  
     current_dists: vector of the current cong district assignments
     
     new_dists: vector of the new cong district assignments
     
   */
	
  beta_population = betas["population"];
  
  NumericVector distswitch;
  
  for(int i = 0; i < cdvec.size(); i++){
    if(is_true(any(distswitch == cdvec(i))) == FALSE){
      distswitch.push_back(cdvec(i));
    }
  }
	  
  // Calculate parity
  double parity = (double) sum(popvec) / (max(current_dists) + 1);

  // Log_e(2)
  double loge2 = log(2.0);

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Population objects
    int pop_new = 0;
    int pop_old = 0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(cdvec == distswitch(i));

    // Get population of the old districts
    for(int j = 0; j < new_cds.size(); j++){
      pop_new += popvec(new_cds(j));
    }
    for(int j = 0; j < current_cds.size(); j++){
      pop_old += popvec(current_cds(j));
    }

    // Calculate the penalty
    psi_new += (double)std::abs((pop_new / parity) - 1);
    psi_old += (double)std::abs((pop_old / parity) - 1);

  }

  // Calculate the ratio
  double ratio = (double) exp(beta_population * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["pop_ratio"] = ratio;
  out["pop_new_psi"] = psi_new;
  out["pop_old_psi"] = psi_old;

  return out;

}






