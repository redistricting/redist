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
    
    NumericVector betas = NumericVector::create(Named("population") = 0.0, Named("compact") = 0.0, 
                                                Named("segregation") = 0.0, Named("similar") = 0.0);

    NumericVector anneals = NumericVector::create(Named("population") = 0, Named("compact") = 0, 
                                                Named("segregation") = 0, Named("similar") = 0);
  
    List constraint_vals;
  
    NumericVector current_dists;
  
    NumericVector weights;
  
    int adjswap = 1

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
     
     current_dists: current vector of congressional district assignments
     
     weights: prior weights on the beta sequence
     
     adjswap: flag - do we want adjacent swaps? default to 1
     
     */
  
  public:
  
    // Constructor for constraint-related values
    void init_constraints(double p, NumericVector b_s, NumericVector b_w, NumericMatrix ssd);
    void init_betavals(NumericVector b);
    void init_annealvals(IntegerVector a);
    
    // Modifiers for constraint-related values
  
    // Function to calculate the strength of the beta constraint for population
    List calc_betapop(arma::vec new_dists);
    
    // Function that applies the Geyer Thompson algorithm for simulated tempering
    List changeBeta(arma::vec betavec,
		                double beta,
		                double constraint,
		                NumericVector weights,
		                int adjswap = 1)
    
      
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
List redist_aList_beta::calc_betapop(arma::vec new_dists)
{

  /* Inputs to function
     
     new_dists: vector of the new cong district assignments
     
   */
	
  beta_population = betas["population"];
  
  NumericVector distswitch;
  
  for(int i = 0; i < current_dists.size(); i++){
    if(is_true(any(distswitch == current_dists(i))) == FALSE){
      distswitch.push_back(current_dists(i));
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
    arma::uvec current_cds = find(current_dists == distswitch(i));

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

// Function that applies the Geyer Thompson algorithm for simulated tempering
List redist_aList_beta::changeBeta(double beta, double constraint)
{
  
  /* Inputs to function 
     
     beta: current value of the beta constraint
     
     constraint: the evaluation of the constraint on the current plan
     
   */
  
  // Find beta in betas
  arma::uvec findBetaVec = find(betas == beta);
  int findBeta = findBetaVec(0);

  // Object to test whether beta is at RHS of vector
  int betaLoc = betas.size() - 1;

  // Get transition probabilities and propose a new beta
  double qij;
  double qji;
  double wi;
  double wj;
  double propBeta;

  // Procedure if conducting adjacent swaps
  if(adjswap == 1){
    if(findBeta == 0){ // At first element in betas
      qij = 1;
      qji = .5;
      wi = weights(0);
      wj = weights(1);
      propBeta = betas(1);
    } else if(findBeta == betaLoc){ // At last element in betas
      qij = 1;
      qji = .5;
      wi = weights(betaLoc);
      wj = weights(betaLoc - 1);
      propBeta = betas(betaLoc - 1);
    } else{ // Anywhere in the middle of betas
      qij = .5;
      qji = .5;
      wi = weights(findBeta);
      arma::vec betaswitch = runif(1);
      if(betaswitch(0) < .5){
	      propBeta = betas(findBeta - 1);
	      wj = weights(findBeta - 1);
      }
      if(betaswitch(0) >= .5){
	      propBeta = betas(findBeta + 1);
	      wj = weights(findBeta + 1);
      }
    }
  } else{
    // Procedure if not conducting adjacent swaps
    // qij = qji in non-adjacent framework, don't have to worry abt end units
    qij = 1;
    qji = 1;

    // Draw element from betavec
    arma::vec rand_randindex = runif(1, 0, 1000000000);
    int randindex = fmod(rand_randindex(0), betaLoc);

    // Weight wi 
    wi = weights(findBeta);

    // Draw the proposed beta value
    if(randindex < findBeta){
      propBeta = betas(randindex);
      wj = weights(randindex);
    } else{
      propBeta = betas(randindex + 1);
      wj = weights(randindex + 1);
    }

  }

  // Accept or reject the proposal
  double mhprobGT = (double) exp(constraint * (propBeta - beta)) * wj / wi * qji / qij;
  if(mhprobGT > 1){
    mhprobGT = 1;
  }
  arma::vec testkeepGT = runif(1);
  int decision = 0;
  if(testkeepGT(0) <= mhprobGT){
    decision++;
    beta = propBeta;
  }

  // Create output
  List out;
  out["beta"] = beta;
  out["mh_decision"] = decision;
  out["mh_prob"] = mhprobGT;

  return out;

}







