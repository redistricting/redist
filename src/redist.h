/////////////////////////////////////
// Author: Rei Yatsuhashi
// Institution: American School In Japan
// Date Created: 2020/06/11
// Date Modified: 2020/06/25
// Purpose: Header file for redist_aList class
/////////////////////////////////////

#ifndef REDIST_ALIST
#define REDIST_ALIST
#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#pragma diagnostic clang pop

class redist_aList;
RCPP_EXPOSED_CLASS(redist_aList)

class redist_aList {
	
  protected: 

    Rcpp::List aList;
    Rcpp::NumericVector cdorigvec;
    Rcpp::NumericVector cdvec;
    Rcpp::NumericVector popvec;
    Rcpp::NumericVector cd_pop_vec;
    double eprob; 
    double mh_prob = 1.0;
    int lambda = 0;
  
  /* Members of class:
	
     aList: full adjacency list of geographic units
     
     cdorigvec: initial vector of congressional district assignments
     
     cdvec: current vector of congressional district assignments
     
     popvec: vector of populations
     
     cd_pop_vec: congressional district populations 
		 
     eprob: edgecut probability
     
     mh_prob: current metropolis-hastings probability
		 
     lambda: parameter for conducting multiple swaps
		 
   */
	
  public:
	
    // Constructor
    void init_values(Rcpp::List a, Rcpp::NumericVector c, Rcpp::NumericVector c2, Rcpp::NumericVector p, Rcpp::NumericVector c3, double e, double m, int l);
    void init_values(Rcpp::List a, Rcpp::NumericVector c, Rcpp::NumericVector c2, Rcpp::NumericVector p, Rcpp::NumericVector c3, double e);

    void set_cdvec(Rcpp::NumericVector c);
    void set_cd_pop_vec(Rcpp::NumericVector c);
    void set_eprob(double e);
    void set_lambda(int l); 
  
    Rcpp::List get_aList();
    Rcpp::NumericVector get_cdvec();
    Rcpp::NumericVector get_cdorigvec();
    Rcpp::NumericVector get_popvec();
    int get_lambda();
    double get_eprob();
        
    // Function to generate initial vector of populations
    Rcpp::NumericVector init_pop(arma::vec cds);
			
    // Function to modify adjacency list to reflect adjacency only within a particular congressional district
    Rcpp::List genAlConn(Rcpp::NumericVector cds);
		
    // Function to identify which precincts lie on the boundary of a congressional district
    Rcpp::NumericVector findBoundary(Rcpp::List conList);
	
    // Function to make unidirectional adjacency list bidirectional
    Rcpp::List add_ties(Rcpp::List adj_list);
	
    // Function to cut edges of adjacency list probabilistically
    Rcpp::List cut_edges(Rcpp::List aList_con);
	
    // Function to generate adjacency graph and count clusters
    int countpartitions(Rcpp::List adj_list);
			
    /* Function to run breadth-first search, returning only sets of 
       connected components that reside on the boundary of the districts
     */
    Rcpp::List bsearch_boundary(Rcpp::List adj_list, arma::vec boundary);
	
    // Function to count number of valid partitions to swap
    int count_valid(Rcpp::List boundarypart, Rcpp::NumericVector prop_cdvec);
	
    // Function to draw p for the number of connected components 
    int draw_p();
		
    // Select a connected component as a proposal
    Rcpp::List propose_partition(Rcpp::List boundary_cc, Rcpp::NumericVector cds_prop);
	
    /* Function to check adjacency of a randomly selected connected component 
       against connected components already selected. Also gets candidate
       congressional district swaps if valid
     */
     Rcpp::List adjcheck_propcd(Rcpp::NumericVector prop_partitions, 
                          Rcpp::NumericVector accepted_partitions,
                          Rcpp::NumericVector cds);
	
    // Function to do the elimination check
    int elim_check(Rcpp::NumericVector prop_partition, Rcpp::NumericVector cds);
	
    // Function to do the split check
    int split_check(Rcpp::NumericVector prop_partitions, Rcpp::List adjcheck_out, Rcpp::NumericVector cds_prop);
	
    // Function to update district populations
    void update_cd_pop_vec(Rcpp::NumericVector prop_partition,
                        Rcpp::NumericVector unitpop_vec,
                        int prop_cd,
                        int curr_cd);
  
    // Function to update the edgecut probability
    void update_eprob(double e); 
  
    // Function to update the metropolis-hastings probability for a swap
    void update_mhprob(Rcpp::NumericVector prop_partition,
                       arma::vec cds,		     
                       int prop_cd);
	  
    // Function to update lambda
    void update_lambda(double l); 
  
    // Function to accept or reject swaps
    int mh_decision(double prob);
		
};

#endif 

// Expose classes to R:
RCPP_MODULE(redist_aList_cpp){
  using namespace Rcpp;
  
  class_<redist_aList>( "redist_aList" )
    .default_constructor("Default Constructor")
    .method("init_pop", &redist_aList::init_pop)
    .method("genAlConn", &redist_aList::genAlConn)
    .method("findBoundary",&redist_aList::findBoundary)
    .method("add_ties", &redist_aList::add_ties)
    .method("cut_edges", &redist_aList::cut_edges)
    .method("countpartitions",&redist_aList::countpartitions)
    .method("bsearch_boundary", &redist_aList::bsearch_boundary)
    .method("count_valid", &redist_aList::count_valid)
    .method("draw_p", &redist_aList::draw_p)
    .method("propose_partition", &redist_aList::propose_partition)
    .method("adjcheck_propcd", &redist_aList::adjcheck_propcd)
    .method("elim_check", &redist_aList::elim_check)
    .method("split_check", &redist_aList::split_check)
    .method("update_cd_pop_vec", &redist_aList::update_cd_pop_vec)
    .method("update_eprob", &redist_aList::update_eprob)
    .method("update_mhprob", &redist_aList::update_mhprob)
    .method("update_lambda", &redist_aList::update_lambda)
    .method("mh_decision", &redist_aList::mh_decision)
  ;
}
//.method("",redist_aList::)
