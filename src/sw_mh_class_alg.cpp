/////////////////////////////////////
// Author: Rei Yatsuhashi
// Institution: American School In Japan
// Date Created: 2020/06/26
// Date Modified: 2020/06/26
// Purpose: swMH() function modified with new redist_aList_beta class
/////////////////////////////////////

// Header files
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <time.h>
#include <R.h>
/* convert to header files: #include "redist_aList.cpp"
#include "redist_aList_beta.cpp" */
#include "sw_mh_helper.h"
#include "make_swaps_helper.h"
#include "constraint_calc_helper.h"
#include "redist_analysis.h"

using namespace Rcpp; 

/* Primary function to run redistricting algorithm. An implementation of 
   Algorithm 1 in Barbu and Zhu (2005) using classes. */
// [[Rcpp::export]]

List swMH(redist_aList_beta region,
	  int nsims,
	  int anneal_beta_population = 0,
	  int anneal_beta_compact = 0,
	  int anneal_beta_segregation = 0,
	  int anneal_beta_similar = 0,
	  int exact_mh = 0,
	  int adapt_eprob = 0,
	  int adapt_lambda = 0)
{

/* Inputs to function:
     
     region: the region, in the form of a class containing all necessary information
     
     nsims: number of simulations
     
     anneal_beta_population: flag for whether to anneal the beta pop parameter
     
     anneal_beta_compact: flag for whether to anneal the beta compactness parameter
     
     anneal_beta_segregation: flag for whether to anneal the beta segregation parameter
     
     anneal_beta_similar: flag for whether to anneal the beta similarity parameter
     
     exact_mh: flag for whether to calculate the exact metropolis-hastings w boundary correction
     
     adapt_eprob:
     
     adapt_lambda: 
     
  */
  
  cdvec = region.get_cdvec();
  cdorigvec = region.get_cdvec();
  popvec = region.get_popvec();
  
  // Preprocess vector of congressional district assignments
  if(min(cdvec) == 1){
    for(int i = 0; i < cdvec.size(); i++){
      cdvec(i)--;
    }
  }

  // Preprocess vector of original congressional district assignments
  if(min(cdorigvec) == 1){
    for(int i = 0; i < cdorigvec.size(); i++){
      cdorigvec(i)--;
    }
  }  
  
  // Get populations of districts
  NumericVector district_pops = region.init_pop(cdvec);
  
  // Get beta strengths
  betas = region.get_betas();
  
  // Get vector of unique district ids
  NumericVector uniquedists;
  for(int i = 0; i < cdvec.size(); i++){
    if(is_true(any(uniquedists == cdvec(i))) == FALSE){
      uniquedists.push_back(cdvec(i));
    }
  }
  
  // Get ssd denominator
  double ssd_denom;
  if(betas["compact"] != 0.0 || anneal_beta_compact == 1){
    ssd_denom = as<double>(region.calc_betacompact(uniquedists,
					    1.0)["compact_new_psi"]);
  }else{
    ssd_denom = 1;
  }
  
  pct_dist_parity = region.get_pct_dist_parity();
  
  // Define parity, min and max popoulations
  double parity = sum(popvec) / (max(cdvec) + 1);
  double dist_parity = parity * pct_dist_parity;
  double min_parity = parity - dist_parity;
  double max_parity = parity + dist_parity;

  // Set counter variable
  int k = 0;
  // For storing beta sequence
  int z = 0;
  // For printing progress
  int nsims_10pct = ceil((double)nsims / 10);
  
  // Store outputted congressional districts
  NumericMatrix cd_store(cdvec.size(), nsims);
  
  // Store metropolis-hastings decisions for swaps
  NumericVector decision_store(nsims);
  NumericVector mhprob_store(nsims);
  int decision_counter = 0;
  
  // Store value of psi for all constraints
  NumericVector psipop_store(nsims);
  NumericVector psicompact_store(nsims);
  NumericVector psisegregation_store(nsims);
  NumericVector psisimilar_store(nsims);
  
  // Store value of p, lambda, weights for all simulations
  NumericVector pparam_store(nsims);
  
  // Store number of connected components along the boundary, boundary weights
  NumericVector boundarypartitions_store(nsims);
  NumericVector boundaryratio_store(nsims);
  
  // Store sequence of betas - geyer thompson
  NumericVector betaseq_store(nsims);
  
  if(anneal_beta_population == 1){
    betaseq_store[z] = betas["population"];
  }
  if(anneal_beta_compact == 1){
    betaseq_store[z] = betas["compact"];
  }
  if(anneal_beta_segregation == 1){
    betaseq_store[z] = betas["segregation"];
  }
  if(anneal_beta_similar == 1){
    betaseq_store[z] = betas["similar"];
  }

  // Iterate up z
  z++;
  
  if(anneal_beta_population == 0 && beta != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), betas["population"]);
  }
  if(anneal_beta_compact == 0 && beta_compact != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), betas["compact"]);
  }
  if(anneal_beta_segregation == 0 && beta_segregation != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), betas["segregation"]);
  }
  if(anneal_beta_similar == 0 && beta_similar != 0.0){
    std::fill(betaseq_store.begin(), betaseq_store.end(), betas["similar"]);
  }

  // Store metropolis-hastings decisions - geyer thompson
  NumericVector decision_betaseq_store(nsims);
  NumericVector mhprob_betaseq_store(nsims);
  
  // Initialize objects
  List aList_con; NumericVector boundary; List swap_partitions;
  List boundary_partitions; List cutedge_lists; int p; List aList_con_prop;
  NumericVector boundary_prop; List boundary_partitions_prop; int decision;
  List get_constraint; List gt_out; NumericVector cdvec_prop; int i;
  
  // Open the simulations
  while(k < nsims){

    /////////////////////////////////////
    // First: determine boundary cases //
    /////////////////////////////////////
    // Modify aList list for connected components within cd
    aList_con = region.genAlConn(cdvec);

    // Get vector of boundary units
    boundary = region.findBoundary(aList_con);

    ///////////////////////////////////////////////////////////////////////////
    // Second: within each congressional district, turn on edges with pr = p //
    ///////////////////////////////////////////////////////////////////////////
    // Continue trying until you get p good swaps
    do{
      
      // First element is connected adjlist, second element is cut adjlist
      cutedge_lists = region.cut_edges(aList_con);
      
      ////////////////////////////////////////////////////////////////////
      // Third: generate a list of connected components within each cd //
      ///////////////////////////////////////////////////////////////////
      /* List of connected partitions after edgecuts - first element is list of 
	      partitions, second element is number of partitions */
      boundary_partitions = region.bsearch_boundary(cutedge_lists["connectedlist"],
					     boundary);

      ///////////////////////////////////////////////////////////////////////
      // Fourth - select several connected components w/ unif distribution //
      ///////////////////////////////////////////////////////////////////////
      // Draw parameter p (number of swaps for iteration of alg) from pois(lambda)
      p = region.draw_p();
      
      /* Loop over p, draw p connected components
      swap_partitions = make_swaps(boundary_partitions["bsearch"], 
				   aList, 
				   cdvec,
				   cdorigvec,
				   popvec,
				   district_pops,
				   grouppopvec,
				   ssdmat,
				   min_parity,
				   max_parity,
				   p,
				   eprob,
				   beta_population,
				   beta_compact,
				   beta_segregation,
				   beta_similar,
				   ssd_denom);

    }while(as<int>(swap_partitions["goodprop"]) == 0); */
      
    // Get new boundary, then get number of partitions
    if(exact_mh == 1){
      aList_con_prop = region.genAlConn(as<NumericVector>(swap_partitions["proposed_partition"]));
      boundary_prop = region.findBoundary(aList_con_prop);
      boundary_partitions_prop = region.bsearch_boundary(cutedge_lists["connectedlist"],
						    boundary_prop);
      
      // Correct npartitions to only include boundary partitions that don't break contiguity
      int nvalid_current = region.count_valid(boundary_partitions["bsearch"], cdvec);
      int nvalid_prop = region.count_valid(boundary_partitions_prop["bsearch"],
			      swap_partitions["proposed_partition"]);
      
      // Construct multiple swaps term
      double p_0;
      double F_pi;
      double F_pi_prime;
      lambda = region.get_lambda();
        
      if(lambda > 0){
        p_0 = R::ppois(0, lambda, 1, 0);         
        F_pi = R::ppois(nvalid_current, lambda, 1, 0) - p_0;
        F_pi_prime = R::ppois(nvalid_prop, lambda, 1, 0) - p_0;
      }else{
        F_pi = 1.0;
        F_pi_prime = 1.0;
      }
      
      // Modify metropolis-hastings ratio
      swap_partitions["mh_prob"] = as<double>(swap_partitions["mh_prob"]) *
        pow((double)nvalid_current / nvalid_prop, (double)p) * (F_pi / F_pi_prime);
      boundaryratio_store(k) = pow((double)nvalid_current / nvalid_prop, (double)p);
        
    }
      
    //////////////////////////////////////////
    // Fifth - Accept with some probability //
    //////////////////////////////////////////
    decision = region.mh_decision(as<double>(swap_partitions["mh_prob"]));

    /////////////////////////////////////////////////////////////
    // Also - for simulated tempering, propose a possible swap //
    /////////////////////////////////////////////////////////////
    if((anneal_beta_population == 1) || (betas["population"] != 0.0)){ 
      get_constraint = region.calc_betapop(as<NumericVector>(swap_partitions["proposed_partition"]));
 
      // Store psi value
      if(decision == 1){
	     psipop_store[k] = as<double>(get_constraint["pop_new_psi"]);
      }else{
	     psipop_store[k] = as<double>(get_constraint["pop_old_psi"]);
      }
    
      if(anneal_beta_population == 1){ // Tempering step
	     if(decision == 1){ // Using value of psi if accepted
	
	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["population"],
			    as<double>(get_constraint["pop_new_psi"]));
	
	     }else{ // Using value of psi if not accepted
	
	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["population"],
			      as<double>(get_constraint["pop_old_psi"]));
	
	     }
      
	     // Change beta
	     region.update_betas(as<double>(gt_out["beta"]), "population");
             betas["population"] = as<double>(gt_out["beta"]);
      
	     // Store the output of geyer thompson
	     if(k < nsims){
	       betaseq_store[z] = beta_population;
	     }
	     decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	     mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);
      
      }
    
    }
      
    if((anneal_beta_compact == 1) || (betas["compact"] != 0.0)){ // If constraining, get value of constraints
      get_constraint = region.calc_betacompact(as<NumericVector>(swap_partitions["proposed_partition"]));
      
      // Get psi value
      if(decision == 1){
	     psicompact_store[k] = as<double>(get_constraint["compact_new_psi"]);
      }else{
	     psicompact_store[k] = as<double>(get_constraint["compact_old_psi"]);
      }
    
      if(anneal_beta_compact == 1){ // Annealing step
	     if(decision == 1){ // Using value of psi if accepted
	
	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["compact"], as<double>(get_constraint["compact_new_psi"]));
	
	     }else{ // Using value of psi if not accepted
	  
	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["population"], as<double>(get_constraint["compact_old_psi"]));
	  
	     }
	
	     // Change beta
       betas["compact"] = as<double>(gt_out["beta"]);
       region.update_betas(betas["compact"], "compact");
	
	     // Store the output of geyer thompson
	     if(k < nsims){
	       betaseq_store[z] = beta_compact;
	     }
	     decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	     mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);

      }

    }
      
    if((anneal_beta_segregation == 1) || (betas["segregation"] != 0.0)){ // If constraining, get value of constraint
      get_constraint = region.calc_betasegregation(as<NumericVector>(swap_partitions["proposed_partition"]));
      
      // Get psi value
      if(decision == 1){
	     psisegregation_store[k] = as<double>(get_constraint["segregation_new_psi"]);
      }else{
	     psisegregation_store[k] = as<double>(get_constraint["segregation_old_psi"]);
      }

      if(anneal_beta_segregation == 1){ // Annealing step
	     if(decision == 1){

	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["segregation"],
			    as<double>(get_constraint["segregation_new_psi"]));
	  
	     }else{ // Use value of psi if not accepted

	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["segregation"],
			     as<double>(get_constraint["segregation_old_psi"]));
	  
	     }
          
	     // Change beta
	     betas["segregation"] = as<double>(gt_out["beta"]);
             region.update_betas(betas["segregation"], "segregation");

	     // Store output of geyer thompson
	     if(k < nsims){
	       betaseq_store[k] = beta_segregation;
	     }
	     decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	     mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);

      }

    }
      
    if((anneal_beta_similar == 1) || (betas["similar"] != 0.0)) { // If constraining, get value of constraint
      get_constraint = region.calc_betasimilar(as<NumericVector>(swap_partitions["proposed_partition"]));

      // Get psi value
      if(decision == 1){
	     psisimilar_store[k] = as<double>(get_constraint["similar_new_psi"]);
      }else{
	     psisimilar_store[k] = as<double>(get_constraint["similar_old_psi"]);
      }

      if(anneal_beta_similar == 1){ // Annealing step
	     if(decision == 1){

	       // Propose swapping beta
	       gt_out = region.changeBeta(betas["similar"]
			    as<double>(get_constraint["similar_new_psi"]));
	  
	     }else{ // Use value of psi if not accepted

	       // Propose swapping beta
	       gt_out = changeBeta(betas["similar"],
			    as<double>(get_constraint["similar_new_psi"]));

	     }

	     // Change beta
	     betas["similar"] = as<double>(gt_out["beta"]);
       region.update_betas(betas["similar"], "similar");

	     // Store output of geyer thompson
	     if(k < nsims){
	       betaseq_store[k] = beta_similar;
	     }
	     decision_betaseq_store[k] = as<int>(gt_out["mh_decision"]);
	     mhprob_betaseq_store[k] = as<double>(gt_out["mh_prob"]);
	
      }

    }
      
    //////////////////////////////////////
    // Six = clean up and store results //
    //////////////////////////////////////
    cdvec_prop = clone(as<NumericVector>(swap_partitions["proposed_partition"]));
    if(decision == 1){
      // Update cds to proposed cds
      cdvec = clone(as<NumericVector>(swap_partitions["proposed_partition"]));
	    region.set_cdvec(cdvec);
      // Update district_pops to proposed district pops
      district_pops = clone(as<NumericVector>(swap_partitions["updated_cd_pops"]));
      region.set_cd_pop_vec(district_pops);
      // Store number of boundary partitions
      boundarypartitions_store[k] = boundary_partitions["npartitions"];
    } else{
      boundarypartitions_store[k] = boundarypartitions_store[k-1];
    }
      
    // Store previous iteration
    for(i = 0; i < cdvec.size(); i++){
      cd_store[k * cdvec.size() + i] = cdvec(i);
    }

    // Store p
    pparam_store[k] = p;

    // Store the decision
    decision_store[k] = decision;
    decision_counter += decision;
    
    mhprob_store[k] = as<double>(swap_partitions["mh_prob"]);

    // Advance k, z
    k++;
    z++;

    // Print Progress
    if(k % nsims_10pct == 0){
      Rcout << (double)k / nsims_10pct * 10 << " percent done." << std::endl;
      if(adapt_lambda == 1){
	     Rcout << "Lambda: " << lambda << std::endl;
      }
      if(adapt_eprob == 1){
	     Rcout << "Edgecut Probability: " << eprob << std::endl;
      }
      Rcout << "Metropolis acceptance ratio: "<< (double) decision_counter / (k-1) << std::endl << std::endl;
    }
  
    // Change eprob, lambda if adaptive
    if(adapt_eprob == 1 || adapt_lambda == 1){
      if(k % 50 == 0){
	     if((double)decision_counter / (k-1) > .4){
	       if(adapt_lambda == 1 && lambda < floor((double)aList.size() / 10)){
	         lambda++;
           region.set_lambda(lambda);
	       }
	       if(adapt_eprob == 1 && eprob < .5){
	         eprob = eprob + .01;
           region.set_eprob(eprob);
	       }
	     }
	     if((double)decision_counter / (k-1) < .2){
	       if(adapt_lambda == 1 && lambda > 0){
	         lambda--;
           region.set_lambda(lambda);
	       }
	       if(adapt_eprob == 1 && eprob > 0){
	         eprob = eprob - .01;
           region.set_eprob(eprob);
	       }
	     }
      }
    }
        
  }
    
  // Get distance from parity of each partition
  NumericVector dist_parity_vec = distParity(cd_store, popvec);

  NumericVector dist_orig_vec = diff_origcds(cd_store, cdorigvec);
  
  // Create list, store outputx
  List out;
  out["partitions"] = cd_store;
  out["distance_parity"] = dist_parity_vec;
  out["distance_original"] = dist_orig_vec;
  out["mhdecisions"] = decision_store;
  out["mhprob"] = mhprob_store;
  out["pparam"] = pparam_store;
  out["beta_sequence"] = betaseq_store;
  out["constraint_pop"] = psipop_store;
  out["constraint_compact"] = psicompact_store;
  out["constraint_segregation"] = psisegregation_store;
  out["constraint_similar"] = psisimilar_store;
  out["boundary_partitions"] = boundarypartitions_store;
  out["boundaryratio"] = boundaryratio_store;
  if((anneal_beta_population == 1) || (anneal_beta_compact == 1) ||
     (anneal_beta_segregation == 1) || (anneal_beta_similar == 1)){
    out["mhdecisions_beta"] = decision_betaseq_store;
    out["mhprob_beta"] = mhprob_betaseq_store;
  }
  if(adapt_eprob == 1){
    out["final_eprob"] = eprob;
  }
  if(adapt_lambda == 1){
    out["final_lambda"] = lambda;
  }
  
  return out;
  
}
    

}



