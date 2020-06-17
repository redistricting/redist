/////////////////////////////////////
// Author: Rei Yatsuhashi
// Institution: American School In Japan
// Date Created: 2020/06/11
// Date Modified: 2020/06/15
// Purpose: Preliminary Class for Redistricting Adjacency
/////////////////////////////////////

// Header files
#include <RcppArmadillo.h>

using namespace Rcpp;

// Class to consolidate methods relating to redistricting and adjacency

class redist_aList {
	
  private: 

    List aList;
    NumericVector cdvec;
    NumericVector popvec;
    NumericVector cd_pop_vec;
    double eprob; 
    double mh_prob = 1.0;
    int lambda = 0;
  
  /* Members of class:
	
     aList: full adjacency list of geographic units
     
     cdvec: initial vector of congressional district assignments
     
     popvec: vector of populations
     
     cd_pop_vec: congressional district populations 
		 
     eprob: edgecut probability
     
     mh_prob: current metropolis-hastings probability
		 
     lambda: parameter for conducting multiple swaps
		 
   */
	
  public:
	
    // Constructor
    void init_values(List a, NumericVector c, NumericVector p, NumericVector c2, double e, double m, int l);
    void init_values(List a, NumericVector c, NumericVector p, NumericVector c2, double e);
	
    // Function to generate initial vector of populations
    NumericVector init_pop(arma::vec cds);
			
    // Function to modify adjacency list to reflect adjacency only within a particular congressional district
    List genAlConn(NumericVector cds);
		
    // Function to identify which precincts lie on the boundary of a congressional district
    NumericVector findBoundary(List conList);
	
    // Function to make unidirectional adjacency list bidirectional
    List add_ties(List adj_list);
	
    // Function to cut edges of adjacency list probabilistically
    List cut_edges(List aList_con);
	
    // Function to generate adjacency graph and count clusters
    int countpartitions(List adj_list);
			
    /* Function to run breadth-first search, returning only sets of 
       connected components that reside on the boundary of the districts
     */
    List bsearch_boundary(List adj_list, arma::vec boundary);
	
    // Function to count number of valid partitions to swap
    int count_valid(List boundarypart, NumericVector prop_cdvec);
	
    // Function to draw p for the number of connected components 
    int draw_p();
		
    // Select a connected component as a proposal
    List propose_partition(List boundary_cc);
	
    /* Function to check adjacency of a randomly selected connected component 
       against connected components already selected. Also gets candidate
       congressional district swaps if valid
     */
     List adjcheck_propcd(NumericVector prop_partitions, 
                          NumericVector accepted_partitions,
                          NumericVector cds);
	
    // Function to do the elimination check
    int elim_check(NumericVector prop_partition, NumericVector cds);
	
    // Function to do the split check
    int split_check();
	
    // Function to update district populations
    void update_cd_pop_vec(NumericVector prop_partition,
                        int prop_cd,
                        int curr_cd);
  
    // Function to update the edgecut probability
    void update_eprob(double e); 
  
    // Function to update the metropolis-hastings probability for a swap
    void update_mhprob(NumericVector prop_partition,
                       arma::vec cds,		     
                       int prop_cd);
	  
    // Function to update lambda
    void update_lambda(double l); 
  
    // Function to accept or reject swaps
    int mh_decision();
		
}

// Constructor
void redist_aList::init_values(List a, NumericVector c, NumericVector p, NumericVector c2, double e, double m, int l) 
{

  aList = a;
  cdvec = c;
  popvec = p;
  cd_pop_vec = c2;
  eprob = e;
  mh_prob = m;
  lambda = l;
  
}

void redist_aList::init_values(List a, NumericVector c, NumericVector p, NumericVector c2, double e) 
{

  aList = a;
  cdvec = c;
  popvec = p;
  cd_pop_vec = c2;
  eprob = e;

}

// Function to generate initial vector of populations
NumericVector redist_aList::init_pop(arma::vec cds)
{

  /* Inputs to function:
	
     cds: Vector of congressional district populations
		 
   */ 

  // Get number of cds
  int ncds = cds.max() + 1;

  // Create container vector
  NumericVector distpop(ncds);

  // Initialize
  int i; int pop; arma::uvec cd_i_ind; int j;

  // Loop through cd assignments
  for(i = 0; i < ncds; i++){

    // Initialize population count
    pop = 0;

    // Get indices of cds 
    cd_i_ind = find(cds == i);
    
    // Loop through cd_i_ind, get population values
    for(j = 0; j < cd_i_ind.n_elem; j++){
      pop += popvec(cd_i_ind(j));
    }

    // Put in distpop
    distpop(i) = pop;

  }

  return distpop;

}

// Function to modify adjacency list to reflect adjacency only within a particular congressional district
List redist_aList::genAlConn(NumericVector cds) 
{

  /* Inputs to function:
	
     cds: vector of congressional district assignments
		 
   */
  
  // Initialize container list
  List alConnected(cds.size());

  // Initialize
  int i; NumericVector avec; int cd_i; int j;

  // Loop through precincts
  for(i = 0; i < cds.size(); i++){

    // For precinct i, get adjacent precincts
    avec = aList(i);
    
    // Get precinct i's congressional district
    cd_i = cds(i);

    // Initialize empty vector
    NumericVector avec_cd;

    // Loop through avec to identify which are in same cd
    for(j = 0; j < avec.size(); j++){
      
      // Check if j'th entry in avec is same cd, add to avec_cd if so
      if(cds(avec(j)) == cd_i) {
        avec_cd.push_back(avec(j));
      }

    }

    // Add to alConnected list
    alConnected(i) = avec_cd;

  }

  return alConnected;

}

// Function to identify which precincts lie on the boundary of a congressional district
NumericVector redist_aList::findBoundary(List conList)
{
	
  /* Inputs to function:
	
     conList: adjacency list of geographic units within cong district
		 
   */
	
  fullList = aList;

  // Initialize container vector of 0's (not boundary) and 1's (boundary)
  NumericVector isBoundary(fullList.size());

  // Initialize inside loop
  NumericVector full; NumericVector conn; int i;

  // Loop through aList
  for(i = 0; i < fullList.size(); i++){

    // Get vectors of full and cd-connected components for precinct i
    full = fullList(i);
    conn = conList(i);

    // Compare lengths - if conn < full, then boundary unit
    if(full.size() > conn.size()){
      isBoundary(i) = 1;
    }
    
  }

  return isBoundary;

}

// Function to make unidirectional adjacency list bidirectional
List redist_aList::add_ties(List adj_list){

  // Initialize
  int i; NumericVector list1; int j; NumericVector list2;
  
  // Loop through vectors in aList
  for(i = 0; i < adj_list.size(); i++){

    // Get i'th entry in list
    list1 = adj_list(i);
    
    // Loop through elements in list1
    for(j = 0; j < list1.size(); j++){

      // Extract adjacency vector for j'th element of i's adjacency list
      list2 = adj_list(list1(j));
      
      // Check if list 2 includes i
      if(is_true(any(list2 == i)) == FALSE){

        // If not included, add to adjacency vector
        list2.push_back(i);

        // Modify aList to include new adjacency vector
        adj_list(list1(j)) = list2;

      }
    
    }
 
  }

  return adj_list;

}


// Function to cut edges of adjacency list probabilistically
List redist_aList::cut_edges(List aList_con)
{

  /* Inputs to function:
	
     aList_con: adjacency list within cong district
		 
   */

  // Create threshold
  double threshold_prob = 1 - eprob;

  // Define lists to store cut-edge and uncut-edge vectors
  List aList_uncut(aList_con.size());
  List aList_cut(aList_con.size());

  // Initialize inside loop
  int i; NumericVector cc_vec_i_all; NumericVector cc_vec_i;
  arma::vec draws;

  // Define list to store output of both lists

  // Loop through elements of aList_con
  for(i = 0; i < aList_con.size(); i++){

    // Extract i'th vector in list
    cc_vec_i_all = aList_con(i);

    // Subset cc_vec_i to elements > i
    cc_vec_i = cc_vec_i_all[cc_vec_i_all > i];

    // For each element in vector, take random draw from [0,1] uniform
    draws = runif(cc_vec_i.size());

    // Create container vectors of cut and uncut edges
    NumericVector cut;
    NumericVector uncut;

    // Loop through elements of cc_vec_i and compare to entry in draws
    for(int j = 0; j < cc_vec_i.size(); j++){
      
      // Compare to threshold_prob - if draws < thresh, cut edge, else uncut
      if(draws(j) < threshold_prob){
        cut.push_back(cc_vec_i(j));
      } else{
        uncut.push_back(cc_vec_i(j));
      }

    }

    // Store vectors in container lists
    aList_uncut(i) = uncut;
    aList_cut(i) = cut;

  }
  
  // Add ties to aList_uncut, aList_cut
  List aList_uncut_bd = add_ties(aList_uncut);
  List aList_cut_bd = add_ties(aList_cut);

  // Return contents
  List out;
  out["connectedlist"] = aList_uncut_bd;
  out["cutedgelist"] = aList_cut_bd;
  
  return out;

}

// Function to generate adjacency graph and count clusters
int redist_aList::countpartitions(List adj_list) 
{   

  //Takes an adjacency list,
  //The vector of subset nodes
  //The number of subset nodes
						
  //initialize connCompVec   
  //Initialize visited indices
  IntegerVector visitedInd(adj_list.size());
  int indexVisit = 0;
  
  //Initialize connected components
  IntegerVector currConnComp(adj_list.size());

  //Initialize the number of connected components
  int numConnComp = 0;
  
  //Loop over nodes
  for(int i = 0; i < adj_list.size(); i++){
    
    //If i has not been visited...
    if(visitedInd[i] == 0){
      
      //List i as visited
      visitedInd[i] = 1;

      //Increase the number of connected components
      numConnComp++;

      //Add i to the connected component list
      currConnComp[indexVisit] = i;
      
      //increase index visit
      indexVisit++;
      
      //Count the number of nodes in the current connected component
      int nodeCount = indexVisit - 1;
      
      //Initialize a stopping variable:
      int toStop = 0;

      //While we don't stop
      while(toStop == 0){
	
        //get the neighbors of the next current comp
        IntegerVector listNeighs = adj_list[currConnComp[nodeCount]];
	
        //If listNeighs does not have length zero...
        int listLength = listNeighs.size();
        if(listLength > 0){
	  
          //Add nodes of listLength to currConnComp
          //and mark nodes as visited
          for(int j = 0; j < listLength; j++){
            if(visitedInd[listNeighs[j]] == 0){
              currConnComp[indexVisit] = listNeighs[j];
              visitedInd[listNeighs[j]] = 1;

              //Increment indexVisit
              indexVisit++;
            }
          }
        }
	
        //Increment nodeCount
        nodeCount++;

        //If currConnComp[nodeCount] is zero, then we must have new connected component
        //Also stop if we have too many guys.
        if(nodeCount == adj_list.size()){
          toStop = 1;
        }
        else if(currConnComp[nodeCount] == 0){
          toStop = 1;
        }
      }
    }
  }
  
  return numConnComp;
  
}

// Function to run breadth-first search, returning only sets of connected components that reside on the boundary of the districts
List redist_aList::bsearch_boundary(List adj_list, arma::vec boundary)
{

  /* Inputs to function:
	
     adj_list: adjacency list
		 
     boundary: vector of boundary element indicators (as arma)
		 
   */

  // Get indices of boundary units
  arma::uvec boundary_indices = find(boundary == 1);

  // Container - outputted of breadth search, a list
  List bsearch;

  // Container - partition vector, gets added to bsearch when queue is empty
  NumericVector partition;

  // Set mark vector - ledger of which indices have been reached
  NumericVector mark(adj_list.size());

  // Set queue vector
  NumericVector q;

  // Initialize breadth search with first element in boundary_indices
  mark(boundary_indices(0)) = boundary_indices(0);
  partition.push_back(boundary_indices(0));
  q = adj_list(boundary_indices(0));

  // Initialize objects inside loop
  int u; bool in_part; NumericVector adj_u; int i; int v; 

  // Begin do{} loop - run until number of elements in boundary_indices is 0
  do{

    // Begin while{} loop - run until q is empty
    while(q.size() > 0){
      
      // Dequeue first element in queue
      u = q(0);

      // Mark that element in ledger
      mark(u) = u;

      // Check if element is in the partition - add to partition if false
      in_part = is_true(any(partition == u));
      if(in_part == false){
        partition.push_back(u);
      }
      
      // Get adjacency vector for unit u
      adj_u = adj_list(u);

      // Loop through elements of adj_u, add to queue and mark if not reached
      if(adj_u.size() > 0){
	
        // Start loop
        for(i = 0; i < adj_u.size(); i++){

          // Reach element v
          v = adj_u(i);

          /* Check if already reached - if false, mark, add to partition, and
             add to queue */
          if(is_true(any(mark == v)) == FALSE){
            mark(v) = v;
            partition.push_back(v);
            q.push_back(v);
          }

        }

      }

      // Erase dequeued element from queue when done searching
      q.erase(q.begin());

    }

    // Handling an empty queue
    if(q.size() == 0){

     /* First, find boundary units that are in the reached partition and
        remove them from boundary_units vector */
       for(i = boundary_indices.n_elem - 1; i >= 0; i--){
         if(is_true(any(partition == boundary_indices(i))) == TRUE){
           boundary_indices.shed_row(i);
         }
       }
   	
       // Store the partition, clear partition vector
       bsearch.push_back(partition);
       partition.erase(partition.begin(), partition.end());

       // Re-initialize breadth search from new starting value if nonempty
       if(boundary_indices.n_elem > 0){
         q = adj_list(boundary_indices(0));
         mark(boundary_indices(0)) = boundary_indices(0);
         partition.push_back(boundary_indices(0));
       }
     }

  } while(boundary_indices.n_elem > 0);

  // Get breadth search size
  int bsearch_size = bsearch.size();

  // Get weight_boundary vector
  double weight_boundary = (double) countpartitions(adj_list) / bsearch_size;

  List out;
  out["bsearch"] = bsearch;
  out["npartitions"] = bsearch_size;
  out["weight_boundary"] = weight_boundary;

  return out;

}

// Function to count number of valid partitions to swap
int redist_aList::count_valid(List boundarypart, NumericVector prop_cdvec)
{

  int cd_boundary; arma::vec part; int j; int i;
  arma::uvec find_cds; int counter = 0;

  for(i = 0; i < boundarypart.size(); i++){
    
    // Get the partition
    part = as<arma::vec>(boundarypart(i));
    
    // Get the congressional district of the boundary
    cd_boundary = prop_cdvec(part(0));
    
    // Find indices within that congressional district
    find_cds = find(as<arma::vec>(prop_cdvec) == cd_boundary);
    
    // Remove elements in the partition from that cd
    NumericVector cd_less_boundary;
    for(j = 0; j < find_cds.n_elem; j++){
      if(any(part == find_cds(j)) == false){
         cd_less_boundary.push_back(find_cds(j));
      }
    }

    // If cd_less_boundary empty, then continue
    // Eliminates district so invalid partition
    if(cd_less_boundary.size() == 0){
      continue;
    }
    
    // Create new adjacency list
    List newadj(cd_less_boundary.size());
    for(j = 0; j < newadj.size(); j++){
      
      // Extract vector from adjacency list
      NumericVector getadjvec = aList(cd_less_boundary(j));
      
      // Subset down to elements in cd_less_boundary
      NumericVector getadjvec_sub;
      for(int k = 0; k < getadjvec.size(); k++){
        if(any(as<arma::vec>(cd_less_boundary) == getadjvec(k))){
           getadjvec_sub.push_back(getadjvec(k));
         }
      }
      
      // Change indices
      NumericVector getadjvec_new;
      for(int k = 0; k < getadjvec_sub.size(); k++){
        arma::uvec ind = find(as<arma::vec>(cd_less_boundary) == getadjvec_sub(k));
        getadjvec_new.push_back(ind(0));
      }
      
      // Add to newadj
      newadj(j) = getadjvec_new;
      
    }
    
    // Calculate number of partitions
    int nparts = countpartitions(newadj);
    if(nparts == 1){
      counter++;
    }
    
  }
  
  return counter;
  
}

// Function to draw p for the number of connected components 
int redist_aList::draw_p()
{
	
  int p;
  if(lambda > 0){
    p = R::rpois(lambda);
    p++;
  } else{
    p = 1;
  }

  return p;

}

// Function to propose a connected component for partition
List redist_aList::propose_partition(List boundary_cc)
{

  arma::vec rand_sample_index = runif(1, 0, 1000000000);
  int sample_index = fmod(rand_sample_index(0), boundary_cc.size());
	
  NumericVector prop_partitions = boundary_cc(sample_index);
  boundary_cc.erase(sample_index);
  curr_cd = cds_prop(prop_partitions(0));

  List out;
  out["prop_partition"] = prop_partitions;
  out["curr_cd"] = curr_cd;
	
  return out;
	
}

List redist_aList::adjcheck_propcd(NumericVector prop_partitions,
                                   NumericVector accepted_partitions,
                                   NumericVector cds)
{  
  /* Inputs to function:
	
     prop_partitions: The proposed partition to be swapped
		 
     accepted_partitions: Vector of district ID's that have been acccepted
		 
     cds: Vector of cd assignments - has to be cds_prop to 
     avoid complications with not recognizing splitting, elimination
		 
   */
  
  // Initialize adjacency check value
  int adj_check = 0;

  // Initialize vector of proposed congressional cds
  NumericVector prop_cds;
  
  // Get current cd of the proposed partition
  int current_cd = cds(prop_partitions(0));

  /* Loop over units in prop_partitions, test to see if any are adjacent to 
     any units in accepted_partitions */
  for(int i = 0; i < prop_partitions.size(); i++){

    // For i'th unit in prop_partitions, get the adjacency list
    NumericVector adj_units = aList(prop_partitions(i));

    // Loop to see if any of the indices in adj_units are in accepted_partitions
    for(int j = 0; j < adj_units.size(); j++){

      // See if element j of adj_units is in accepted_partitions
      bool test_adj = is_true(any(accepted_partitions == adj_units(j)));

      /* If true, iterate adj_check to 1 and break the loop to throw out 
	 			 the partition */
      if(test_adj == TRUE){
        adj_check++;
        break;
      }

      // If not true, then look at the congressional districts of adjacent units
      // Is unit j's congressional district equal to current_cd?
      bool same_cd = cds(adj_units(j)) == current_cd;
      // Would this be a new addition to prop_cds?
      bool new_cd = is_true(any(prop_cds == cds(adj_units(j))));

      // If both conditions are false, add to prop_cds
      if((same_cd == FALSE) && (new_cd == FALSE)){
        prop_cds.push_back(cds(adj_units(j)));
      }

    }

    // If partition is found to be adjacent, do not test any more
    if(adj_check == 1){
      break;
    }

  }

  // Create output from function
  List out;
  out["adjacency_check"] = adj_check;
  out["proposed_cds"] = prop_cds;

  return out;

}

// Function to do the elimination check
int redist_aList::elim_check(NumericVector prop_partition, NumericVector cds)
{

  /* Inputs to function:
	
     prop_partition: Proposed partition
		 
     cds: Vector of congressional districts - using accepted partitions
		 
   */ 
  
  // Indicator for elimimation
  int elimcheck = 0;

  // Get current congressional district
  int current_cd = cds(prop_partition(0));

  // Get length of cds that is of that cd
  NumericVector subcd = cds[cds == current_cd];
  int current_cd_size = subcd.size();

  // If prop_partition size is equal to number of units of that cd, then elim.
  if(current_cd_size == prop_partition.size()){
    elimcheck++;
  }

  return elimcheck;

}

// Function to do the split check
int redist_aList::split_check(List adjcheck_out, NumericVector cds_prop)
{
	
  // Indicator for elimination
  int splitcheck = 0;
	
  NumericVector possible_cd_swaps = as<NumericVector>(adjcheck_out["proposed_cds"]);
  NumericVector cds_splittest = clone(cds_prop);

  for(int j = 0; j < prop_partitions.size(); j++){
    cds_splittest(prop_partitions(j)) = possible_cd_swaps(0);
  }
	
  // Get adjacency list
  List aList_testsplit = genAlConn(aList, cds_splittest);
	
  // Get number of connected components
  int num_cds = countpartitions(aList_testsplit);
  
  if(num_cds != ndists){
    splitcheck = 1;
  }
	
}

// Function to update district populations
void redist_aList::update_cd_pop_vec(NumericVector prop_partition,
			          int prop_cd,
			          int curr_cd)
{

  /* Inputs to function:
	
     prop_partition: Proposed partition to be swapped
		 
     prop_cd: proposed cong district for prop_partition
		 
     curr_cd: old cong district for prop_partition
		 
   */

  // Clone distpop_vec
  NumericVector distpop_vec_clone = clone(cd_pop_vec);
  
  // Current population, proposed district population
  int currpop = distpop_vec_clone(curr_cd);
  int proppop = distpop_vec_clone(prop_cd);

  // Loop through prop_partition
  for(int i = 0; i < prop_partition.size(); i++){
    currpop -= unitpop_vec(prop_partition(i));
    proppop += unitpop_vec(prop_partition(i));
  }

  // Put back in distpop_vec
  distpop_vec_clone(curr_cd) = currpop;
  distpop_vec_clone(prop_cd) = proppop;

  cd_pop_vec = distpop_vec_clone;

}

// Function to update the edgecut probability
void update_eprob(double e)
{
  
  eprob = e;
  
}

// Function to update the metropolis-hastings probability for a swap
void redist_aList::update_mhprob(NumericVector prop_partition,
		                 arma::vec cds,		     
		                 int prop_cd)
{

  /* Inputs to function:
	
     prop_partition: Proposed partition to swap
		 
     cds: original vec of cong districts
		 
     prop_cd: proposed cong district
		 
   */

  // Initialize c1, c2
  int c1 = 0;
  int c2 = 0;
	
  // Loop through prop_partition
  for(int i = 0; i < prop_partition.size(); i++){

    // Get adjacency vector
    NumericVector adj_vec = aList(prop_partition(i));

    // Loop throgh elements of adj_vec
    for(int j = 0; j < adj_vec.size(); j++){

      /* Calculate C(V_0, V_l' \ V_0) - add 1 if adjacent to switched
	 partition, and your old cd assignment is same as proposed switch */
      if(cds(adj_vec(j)) == prop_cd){
        c1++;
      }

      /* Calculate C(V_0, V_l \ V_0) - add 1 if you are adjacent to the 
	 proposed switch, if your cd assignment is the same as the old cong
	 district, and you are not in the switch partition */
      if((cds(adj_vec(j)) == cds(prop_partition(0))) && (is_true(any(prop_partition == adj_vec(j))) == FALSE)){
        c2++;
      }

    }

  }

  // Recalculate mh probability
  mh_prob = (double) mh_prob * ((double)pow(1 - eprob, c1) / pow(1 - eprob, c2));

}

// Function to update lambda
void update_lambda(double l)
{
  
  lambda = l;
  
}

// Function to accept or reject swaps
int redist_aList::mh_decision()
{
  
  // Initialize decision
  int decision = 0;

  // Get acceptance probability
  double acc_prob;
  if(mh_prob < 1){
    acc_prob = mh_prob;
  } else{
    acc_prob = 1;
  }

  // Draw from uniform
  arma::vec draw_prob = runif(1);

  // Make decision
  if(draw_prob(0) <= acc_prob){
    decision++;
  }

  return decision;

}
