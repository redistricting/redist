/////////////////////////////////////
// Author: Rei Yatsuhashi
// Institution: American School In Japan
// Date Created: 2020/06/27
// Date Modified: 2020/06/27
// Purpose: Header file for redist_aList_beta class
/////////////////////////////////////

#ifndef REDIST_ALIST_B
#define REDIST_ALIST_B

#include <RcppArmadillo.h>
#include "redist.h"

// Class to consolidate methods relating to the constraints and tempering for redistricting

class redist_aList_beta: public redist_aList {
  
  protected:
    
    double pct_dist_parity;
  
    Rcpp::NumericVector grouppopvec;
  
    Rcpp::NumericVector beta_sequence;
    
    Rcpp::NumericMatrix ssdmat;
    
    Rcpp::NumericVector beta_weights = Rcpp::NumericVector::create(Rcpp::Named("population") = 0.0, Rcpp::Named("compact") = 0.0, 
                                                Rcpp::Named("segregation") = 0.0, Rcpp::Named("similar") = 0.0);
  
    Rcpp::NumericVector distswitch;
  
    int adjswap = 1;

    /* Inputs to function:
     
     pct_dist_parity: strength of population parity requirement
     
     grouppopvec: vector of subgroup populations for each geographic unit
     
     beta_sequence: sequence of betas to anneal over
     
     ssdmat: matrix of squared distances between geographic units. For constraining
     on compactness
     
     beta_weights: {beta_population, beta_compact, beta_segregation, beta_similar}
     
       beta_population: strength of constraint for achieving population parity
     
       beta_compact: strength of constraint for achieving district compactness
     
       beta_segregation: strength of constraint for packing group into district
     
       beta_similar: strength of constraint for examining plans similar to original district
     
     distswitch: vector containing the old district, and the proposed new district
     
     adjswap: flag - do we want adjacent swaps? default to 1
     
     */
  
  public:
  
    // Constructor for constraint-related values
    redist_aList_beta(double p, Rcpp::NumericVector b_s, Rcpp::NumericVector b_w, Rcpp::NumericMatrix ssd, Rcpp::NumericVector d);
    
    void update_current_dists(Rcpp::NumericVector c);
    void update_distswitch(); 
    Rcpp::NumericVector get_grouppopvec();
    Rcpp::NumericVector get_ssdmat();
    Rcpp::NumericVector get_beta_sequence();
    Rcpp::NumericVector get_beta_weights();
    double get_pct_dist_parity();

    void update_weights(double b, std::string s);

    // Function that applies the Geyer Thompson algorithm for simulated tempering
    Rcpp::List changeBeta(double beta, double constraint);
  
    // Function to calculate the strength of the beta constraint for population
    Rcpp::List calc_betapop(arma::vec new_dists);
    
    // Function to calculate the strength of the beta constraint for compactness
    // Fryer and Holden 2011 RPI index
    Rcpp::List calc_betacompact(arma::vec new_dists, double denominator = 1.0);
	
    // Function to constrain by segregating a group
    Rcpp::List calc_betasegregation(arma::vec new_dists);
  
    // Function to constrain on plan similarity to original plan
    Rcpp::List calc_betasimilar(arma::vec new_dists);

};

#endif
