#' Identify Cores of a District (Heuristic)
#'
#' Creates a grouping ID to unite geographies and perform analysis on a smaller
#' set of precincts. It identifies all precincts more than \code{boundary} edges
#' of a district district boundary. Each contiguous group of precincts more than 
#' \code{boundary} steps away from another district gets it own group. Some 
#' districts may have multiple, disconnected components that make up the core, 
#' but each of these is assigned a separate grouping id so that a call to 
#' \code{sf::st_union()} would produce only connected pieces.
#'
#' This is a loose interpretation of the
#' \href{https://www.ncsl.org/research/redistricting/redistricting-criteria.aspx}{NCSL's summary}
#' of redistricting criteria to preserve the cores of prior districts. Using the
#' adjacency graph for a given plan, it will locate the precincts on the
#' boundary of the district, within \code{boundary} steps of the edge. Each of
#' these  is given their own group. Each remaining entry that is not near the
#' boundary of the district is given an id that can be used to group the
#' remainder of the district by connected component. This portion is deemed the
#' core of the district.
#'
#' @param adj zero indexed adjacency list.
#' @param plan An integer vector or matrix column of district assignments.
#' @param boundary Number of steps to check for. Defaults to 1.
#' @param focus Optional. Integer. A single district to focus on.
#' @param simplify Optional. Logical. Whether to return extra information or just grouping ID.
#' @param adjacency Deprecated, use adj. A zero indexed adjacency list.
#' @param district_membership Deprecated, use plan. An integer vector or matrix column of district assignments.
#'
#' @return integer vector (if simplify is false). Otherwise it returns a tibble with the grouping
#' variable as \code{group_id} and additional information on connected components.
#'
#' @importFrom dplyr row_number cur_group_id
#'
#' @seealso [redist.plot.cores()] for a plotting function
#' @concept prepare
#' @export
#'
#' @examples
#' data(fl250)
#' adj <- redist.adjacency(fl250)
#' plan <- redist.smc(adj = adj, total_pop = fl250$pop,
#'                   nsims = 1, ndists = 4)$plans
#' core <- redist.identify.cores(adj = adj, plan = plan)
#' redist.plot.cores(shp = fl250, plan = plan, core = core)
#'
redist.identify.cores <- function(adj, plan, boundary = 1, focus = NULL, simplify = TRUE,
                                  adjacency, district_membership){
  if(!missing(adjacency)){
    adj <- adjacency
    .Deprecated(new = 'adj', old = 'adjacency')
  }
  if(!missing(district_membership)){
    plan <- district_membership
    .Deprecated(new = 'plan', old = 'district_membership')
  }


  if(missing(adj)){
    stop('Please provide an object to adj.')
  }

  if(! 'list' %in% class(adj)){
    stop('adj must be an adj list.')
  }

  if(missing(plan)){
    stop('Please provide an object to plan.')
  }
  if('matrix' %in% class(plan)){
    plan <- plan[,1]
  }


  # init a nice empty list
  cd_within_k <- lapply(1:length(plan), FUN = function(x){integer(0)})

  core <- cores(adj = adj, dm = plan, k = boundary, cd_within_k = cd_within_k)

  if(!is.null(focus)){
    idx <- unlist(lapply(core$cd_within_k, FUN = function(x){focus %in% x})) | plan == focus

    core$k <- ifelse(idx, core$k, 0)

    conncomp <- update_conncomp(dm  = core$dm, kvec = core$k, adj = adj)
    core$conncomp <- conncomp
  }

  tb <- tibble(dm = plan, boundary = core$k, cc = core$conncomp) %>%
    group_by(dm, boundary) %>%
    mutate(gid = row_number()) %>%
    ungroup() %>%
    mutate(gid = ifelse(boundary == 0, cc, gid)) %>%
    mutate(gid = paste0(dm, '-', boundary, '-', gid)) %>% group_by(gid) %>%
    mutate(group = cur_group_id()) %>% ungroup()

  gid <- tb$group


  if(simplify){
    return(gid)
  } else{
    core$group_id <- gid
    return(core)
  }


}

#' Uncoarsen a District Matrix
#'
#' After a cores analysis or other form of coarsening, sometimes you need
#' to be at the original geography level to be comparable. This takes in a
#' coarsened matrix and uncoarsens it to the original level
#'
#' @param plans A coarsened matrix of plans.
#' @param district_membership Deprecated, use plans. A coarsened district membership matrix
#' @param group_index The index used to coarsen the shape.
#'
#' @return matrix
#'
#' @concept post
#' @export
redist.uncoarsen <- function(plans, district_membership, group_index){
  if(!missing(district_membership)){
    plan <- district_membership
    .Deprecated(new = 'plan', old = 'district_membership')
  }


  uncoarse <- matrix(nrow = length(group_index),
                     ncol = ncol(plans))

  remain <- sort(unique(group_index))


  for(i in 1:length(group_index)){
    uncoarse[i,] <- plans[which(group_index[i] == remain),]
  }

  return(uncoarse)
}


globalVariables(c('dm', 'cc'))
