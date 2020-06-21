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
  
    NumericVector grouppopvec;
  
    NumericVector beta_sequence;
  
    NumericVector beta_weights;
    
    NumericMatrix ssdmat;
    
    NumericVector betas = NumericVector::create(Named("population") = 0.0, Named("compact") = 0.0, 
                                                Named("segregation") = 0.0, Named("similar") = 0.0);

    NumericVector anneals = NumericVector::create(Named("population") = 0, Named("compact") = 0, 
                                                  Named("segregation") = 0, Named("similar") = 0);
  
    List constraint_vals;
  
    NumericVector current_dists;
  
    NumericVector distswitch;
  
    int adjswap = 1;

    /* Inputs to function:
     
     pct_dist_parity: strength of population parity requirement
     
     grouppopvec: vector of subgroup populations for each geographic unit
     
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
     
     distswitch: vector containing the old district, and the proposed new district
     
     adjswap: flag - do we want adjacent swaps? default to 1
     
     */
  
  public:
  
    // Constructor for constraint-related values
    void init_constraints(double p, NumericVector b_s, NumericVector b_w, NumericMatrix ssd);
    void init_betavals(NumericVector b);
    void init_annealvals(IntegerVector a);
    
    void update_current_dists(NumericVector c);
    void update_distswitch();

    // Modifiers for constraint-related values

    // Function that applies the Geyer Thompson algorithm for simulated tempering
    List changeBeta(double beta, double constraint);
  
    // Function to calculate the strength of the beta constraint for population
    List calc_betapop(arma::vec new_dists);
    
    // Function to calculate the strength of the beta constraint for compactness
    // Fryer and Holden 2011 RPI index
    List calc_betacompact(arma::vec new_dists, NumericMatrix ssdmat, double denominator = 1.0);
	
    // Function to constrain by segregating a group
    List calc_betasegregation(arma::vec new_dists);
  
    // Function to constrain on plan similarity to original plan
    List calc_betasimilar(arma::vec new_dists);

}

// Constructor
void redist_aList_beta::init_constraints(double p, NumericVector b_s, NumericVector b_w, NumericMatrix ssd) 
{
  
  pct_dist_parity = p;
  beta_sequence = b_s;
  beta_weights= b_w;
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

void redist_aList_beta::update_current_dists(NumericVector c)
{
  
  current_dists = c;
  
}

void redist_aList_beta::update_distswitch()
{

  for(int i = 0; i < current_dists.size(); i++){
    if(is_true(any(distswitch == current_dists(i))) == FALSE){
      distswitch.push_back(current_dists(i));
    }
  }
  
}

// Function that applies the Geyer Thompson algorithm for simulated tempering
List redist_aList_beta::changeBeta(double beta, double constraint)
{
  
  /* Inputs to function 
     
     beta: current value of the beta constraint
     
     constraint: the evaluation of the constraint on the current plan
     
   */
  
  // Find beta in betas
  arma::uvec findBetaVec = find(beta_sequence == beta);
  int findBeta = findBetaVec(0);

  // Object to test whether beta is at RHS of vector
  int betaLoc = beta_sequence.size() - 1;

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
      wi = beta_weights(0);
      wj = beta_weights(1);
      propBeta = beta_sequence(1);
    } else if(findBeta == betaLoc){ // At last element in betas
      qij = 1;
      qji = .5;
      wi = beta_weights(betaLoc);
      wj = beta_weights(betaLoc - 1);
      propBeta = beta_sequence(betaLoc - 1);
    } else{ // Anywhere in the middle of betas
      qij = .5;
      qji = .5;
      wi = beta_weights(findBeta);
      arma::vec betaswitch = runif(1);
      if(betaswitch(0) < .5){
        propBeta = beta_sequence(findBeta - 1);
        wj = beta_weights(findBeta - 1);
      }
      if(betaswitch(0) >= .5){
        propBeta = beta_sequence(findBeta + 1);
        wj = beta_weights(findBeta + 1);
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
    wi = beta_weights(findBeta);

    // Draw the proposed beta value
    if(randindex < findBeta){
      propBeta = beta_sequence(randindex);
      wj = beta_weights(randindex);
    } else{
      propBeta = beta_sequence(randindex + 1);
      wj = beta_weights(randindex + 1);
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

// Function to calculate the strength of the beta constraint for population
List redist_aList_beta::calc_betapop(arma::vec new_dists)
{

  /* Inputs to function
     
     new_dists: vector of the new cong district assignments
     
   */
	
  beta_population = betas["population"];
  
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

// Function to calculate the strength of the beta constraint for compactness
// Fryer and Holden 2011 RPI index
List calc_betacompact(arma::vec new_dists,
		      NumericMatrix ssdmat,
		      double denominator = 1.0){

  /* Inputs to function:
  
     new_dists: vector of the new cong district assignments
     
     ssdmat: squared distance matrix
     
     denominator: normalizing constant for rpi
     
   */
  
  beta_compact = betas["compact"];
  
  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Log_e(2)
  double loge2 = log(2.0);

  // Loop over the congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    double ssd_new = 0.0;
    double ssd_old = 0.0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // SSD for new partition
    for(int j = 0; j < new_cds.size(); j++){
      for(int k = j + 1; k < new_cds.size(); k++){
	      ssd_new += (double) ssdmat(new_cds(j),new_cds(k)) *
	        popvec(new_cds(j)) * popvec(new_cds(k));
      }
    }

    // SSD for old partition
    for(int j = 0; j < current_cds.size(); j++){
      for(int k = j + 1; k < current_cds.size(); k++){
	      ssd_old += (double) ssdmat(current_cds(j),current_cds(k)) *
	      popvec(current_cds(j)) * popvec(current_cds(k));
      }
    }

    // Add to psi
    psi_new += ssd_new;
    psi_old += ssd_old;

  }

  // Normalize psi
  psi_new = (double) psi_new / denominator;
  psi_old = (double) psi_new / denominator;

  // Calculate ratio
  double ratio = (double) exp(beta_compact * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["compact_ratio"] = ratio;
  out["compact_new_psi"] = psi_new;
  out["compact_old_psi"] = psi_old;

  return out;

}

// Function to constrain by segregating a group
List redist_aList_beta::calc_betasegregation(arma::vec new_dists)
{

  /* Inputs to function:

     new_dists: vector of the new cong district assignments
     
  */
  
  beta_segregation = betas["segregation"];

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Log_e(2)
  double loge2 = log(2.0);

  // Initialize denominator
  int T = sum(popvec);
  double pAll = (double) sum(grouppopvec) / T;
  double denom = (double)2 * T * pAll * (1 - pAll);
  
  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    int oldpopall = 0;
    int newpopall = 0;
    int oldpopgroup = 0;
    int newpopgroup = 0;
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));
  
    // Segregation for proposed assignments
    for(int j = 0; j < new_cds.size(); j++){
      newpopall += popvec(new_cds(j));
      newpopgroup += grouppopvec(new_cds(j));
    }
  
    // Segregation for current assignments
    for(int j = 0; j < current_cds.size(); j++){
      oldpopall += popvec(current_cds(j));
      oldpopgroup += grouppopvec(current_cds(j));
    }
  
    // Calculate proportions
    // Rcout << "old population group " << oldpopgroup << std::endl;
    // Rcout << "old population all " << oldpopall << std::endl;
    double oldgroupprop = (double) oldpopgroup / oldpopall;
    // Rcout << "old proportion group " << oldgroupprop << std::endl;
    double newgroupprop = (double) newpopgroup / newpopall;

    // Get dissimilarity index
    psi_new += (double)(newpopall * std::abs(newgroupprop - pAll));
    psi_old += (double)(oldpopall * std::abs(oldgroupprop - pAll));

  }
  
  // Standardize psi
  psi_new = (double) psi_new / denom;
  psi_old = (double) psi_old / denom;

  // Get mh ratio
  double ratio = (double) exp(beta_segregation * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["segregation_ratio"] = ratio;
  out["segregation_new_psi"] = psi_new;
  out["segregation_old_psi"] = psi_old;

  return out;

}

// Function to constrain on plan similarity to original plan
List redist_aList_beta::calc_betasimilar(arma::vec new_dists)
{

  /* Inputs to function:
  
     new_dists: vector of the new cong district assignments

   */
  
  beta_similar = betas["similar"];
  
  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Log_e(2)
  double loge2 = log(2.0);

  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    int new_count = 0;
    int old_count = 0;
    NumericVector orig_cds = wrap(find(cdvec == distswitch(i)));
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // Similarity measure for proposed assignments
    for(int j = 0; j < new_cds.size(); j++){
      if(any(cdvec == new_cds(j)).is_true()){
    	new_count++;
      }
    }

    // Similarity measure for current assignments
    for(int j = 0; j < current_cds.size(); j++){
      if(any(cdvec == current_cds(j)).is_true()){
    	old_count++;
      }
    }

    // Calculate proportions
    double old_count_prop = (double) old_count / cdvec.size();
    double new_count_prop = (double) new_count / cdvec.size();
    
    // Add to psi
    psi_new += (double) std::abs(new_count_prop - 1);
    psi_old += (double) std::abs(old_count_prop - 1);

  }

  // Normalize by dividing by number of congressional districts
  psi_new = psi_new / distswitch.size();
  psi_old = psi_old / distswitch.size();

  // Get MH ratio
  double ratio = (double) exp(beta_similar * loge2 * (psi_new - psi_old));

  // Create return object
  List out;
  out["similar_ratio"] = ratio;
  out["similar_new_psi"] = psi_new;
  out["similar_old_psi"] = psi_old;

  return out;
}



