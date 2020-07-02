///////////////////////////////////////////////
// Author: Ben Fifield
// Institution: Princeton University
// Date Created: 2014/12/26
// Date Last Modified: 2015/02/26
// Purpose: Contains functions to run calculate beta constraints
/////////////////////////////////////////////// 

// Header files
#include <RcppArmadillo.h>
#include "make_swaps_helper.h"

using namespace Rcpp;

// Function to calculate the strength of the beta constraint for population
List calc_psipop(arma::vec current_dists,
		 arma::vec new_dists,
		 NumericVector pops,
		 NumericVector distswitch)
{

  /* Inputs to function 
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     weight_population: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
  */

  // Calculate parity
  double parity = (double)sum(pops) / (max(current_dists) + 1);

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
      pop_new += pops(new_cds(j));
    }
    for(int j = 0; j < current_cds.size(); j++){
      pop_old += pops(current_cds(j));
    }

    // Calculate the penalty
    psi_new += std::pow(pop_new / parity - 1.0, 2.0);
    psi_old += std::pow(pop_old / parity - 1.0, 2.0);

  }

  // Create return object
  List out;
  out["pop_new_psi"] = std::sqrt(psi_new);
  out["pop_old_psi"] = std::sqrt(psi_old);

  return out;

}

// Function to calculate the strength of the beta constraint for compactness
// Currently implemented: Fryer and Holden 2011 RPI index, Polsby-Popper
List calc_psicompact(arma::vec current_dists,
		     arma::vec new_dists,
		     NumericVector distswitch,
		     std::string measure,
		     // For Polsby-Popper
		     List aList,
		     NumericVector areas_vec,
		     arma::mat borderlength_mat,
		     bool discrete,
		     // For Fryer Holden
		     NumericVector pops,
		     NumericMatrix ssdmat,
		     double denominator = 1.0){

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     beta_compact: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
     ssdmat: squared distance matrix
     denominator: normalizing constant for rpi
  */
  
  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Initialize lists and boundary vectors
  List aList_new;
  List aList_current;
  NumericVector boundarylist_new(new_dists.size());
  NumericVector boundarylist_current(current_dists.size());
  if(measure == "polsby-popper"){
    aList_new = genAlConn(aList, NumericVector(new_dists.begin(), new_dists.end()));
    aList_current = genAlConn(aList, NumericVector(current_dists.begin(), current_dists.end()));
    boundarylist_new = findBoundary(aList, aList_new);
    boundarylist_current = findBoundary(aList, aList_current);
  }

  // Loop over the congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    if(measure == "fryer-holden"){

      List fh_out = fh_compact(new_cds, current_cds, pops, ssdmat, denominator);
      
      // Add to psi
      psi_new += as<double>(fh_out["ssd_new"]);
      psi_old += as<double>(fh_out["ssd_old"]);
      
    }else if(measure == "polsby-popper"){

      List pp_out = pp_compact(new_cds, current_cds,
			       as<arma::vec>(areas_vec),
			       as<arma::vec>(boundarylist_new),
			       as<arma::vec>(boundarylist_current),
			       borderlength_mat,
			       pops,
			       aList,
			       discrete);

      // Add to psi
      psi_new += as<double>(pp_out["pp_new"]);
      psi_old += as<double>(pp_out["pp_old"]);
      
    }

  }

  // Create return object
  List out;
  out["compact_new_psi"] = psi_new;
  out["compact_old_psi"] = psi_old;

  return out;

}

// Function to constrain by segregating a group
List calc_psisegregation(arma::vec current_dists,
			 arma::vec new_dists,
			 NumericVector pops,
			 NumericVector distswitch,
			 NumericVector grouppop)
{

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     pops: vector of district populations
     beta_segregation: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
     grouppop: vector of subgroup district populations
     
  */

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Initialize denominator
  int T = sum(pops);
  double pAll = (double)sum(grouppop) / T;
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
      newpopall += pops(new_cds(j));
      newpopgroup += grouppop(new_cds(j));
    }
  
    // Segregation for current assignments
    for(int j = 0; j < current_cds.size(); j++){
      oldpopall += pops(current_cds(j));
      oldpopgroup += grouppop(current_cds(j));
    }
  
    // Calculate proportions
    // Rcout << "old population group " << oldpopgroup << std::endl;
    // Rcout << "old population all " << oldpopall << std::endl;
    double oldgroupprop = (double)oldpopgroup / oldpopall;
    // Rcout << "old proportion group " << oldgroupprop << std::endl;
    double newgroupprop = (double)newpopgroup / newpopall;

    // Get dissimilarity index
    psi_new += (double)(newpopall * std::abs(newgroupprop - pAll));
    psi_old += (double)(oldpopall * std::abs(oldgroupprop - pAll));

  }
  
  // Standardize psi
  psi_new = (double)psi_new / denom;
  psi_old = (double)psi_old / denom;

  // Create return object
  List out;
  out["segregation_new_psi"] = psi_new;
  out["segregation_old_psi"] = psi_old;

  return out;

}

// Function to constrain on plan similarity to original plan
List calc_psisimilar(arma::vec current_dists,
		     arma::vec new_dists,
		     arma::vec orig_dists,
		     NumericVector distswitch)
{

  /* Inputs to function:
     current_dists: vector of the current cong district assignments
     new_dists: vector of the new cong district assignments
     orig_dists: vector of the true congressional district assignments
     beta_similar: strength of the beta constraint
     distswitch: vector containing the old district, and the proposed new district
  */

  // Initialize psi values
  double psi_new = 0.0;
  double psi_old = 0.0;

  // Loop over congressional districts
  for(int i = 0; i < distswitch.size(); i++){

    // Initialize objects
    int new_count = 0;
    int old_count = 0;
    NumericVector orig_cds = wrap(find(orig_dists == distswitch(i)));
    arma::uvec new_cds = find(new_dists == distswitch(i));
    arma::uvec current_cds = find(current_dists == distswitch(i));

    // Similarity measure for proposed assignments
    for(int j = 0; j < new_cds.size(); j++){
      if(any(orig_cds == new_cds(j)).is_true()){
    	new_count++;
      }
    }

    // Similarity measure for current assignments
    for(int j = 0; j < current_cds.size(); j++){
      if(any(orig_cds == current_cds(j)).is_true()){
    	old_count++;
      }
    }

    // Calculate proportions
    double old_count_prop = (double)old_count / orig_cds.size();
    double new_count_prop = (double)new_count / orig_cds.size();
    
    // Add to psi
    psi_new += (double)std::abs(new_count_prop - 1);
    psi_old += (double)std::abs(old_count_prop - 1);

  }

  // Normalize by dividing by number of congressional districts
  psi_new = psi_new / distswitch.size();
  psi_old = psi_old / distswitch.size();

  // Create return object
  List out;
  out["similar_new_psi"] = psi_new;
  out["similar_old_psi"] = psi_old;

  return out;
}
