author <- "Ron Ammar <ron.ammar@bms.com>"
company <- "Bristol-Myers Squibb Co."
date <- as.Date("2015-11-10")
version <- "0.1.0"
summary <- "Collective Influence algorithm of Morone & Makse (2015) in R"

library(assertthat)
library(igraph)
library(parallel)

plotg <- function(g, layout=layout_with_kk(g), community=FALSE) {
  # Convenience function for plotting graphs.
  if (community) {
    plot(multilevel.community(g), g, vertex.size=20, vertex.label.cex=0.8,
         edge.width=2, layout=layout)
  } else {
    plot(g, vertex.size=20, vertex.color="darkolivegreen3", 
         vertex.label.cex=0.8, edge.width=2, layout=layout)
  }
}

exampleGraph <- function() {
  # Create example graph from Kovacs and Barabasi "News & Views" image.
  # The influencers are (in descending order of influence): n10, n17, n5, n20
  g <- graph(
    c(1,4, 2,4, 3,4, 3,5, 4,5, 4,7, 5,6, 5,29, 6,7, 6,10, 6,20, 7,8, 7,9, 7,10, 10,11, 10,30, 11,12, 11,17, 12,13, 12,14, 12,15, 12,16, 14,15, 16,17, 17,18, 17,29, 18,19, 18,20, 20,21, 21,22, 21,23, 21,24, 21,27, 21,30, 22,23, 24,25, 24,26, 24,27, 24,29, 26,27, 27,28, 27,29),
    directed=FALSE
  )
  
  # Give the nodes names. This is useful because as we remove nodes, the numbering
  # of nodes will change but the names will remain the same.
  g <- set_vertex_attr(g, "name", value=paste("n", as.vector(V(g)), sep=""))
  
  return(g)
}

collectiveInfluence <- function(g, d=3) {
  # Collective Influence as described by Morone & Makse (2015).
  # In simple terms, it is the product of the reduced degree (degree - 1) of a
  # node and the total reduced degree of all nodes at a distance d from the
  # node.
  #
  # Args:
  #   g: the igraph graph
  #   d: the distance, expressed in # of steps, from a given node (default=3)
  #      Distance must be > 0.
  #      According to Morone & Makse, optimal results can be reached at d=3,4,
  #      but this depends on the size/"radius" of the network
  #      NOTE: the distance d is not inclusive. This means that nodes at a
  #      distance of 3 from our node-of-interest do not include nodes at
  #      distances 1 and 2. Only 3.
  #
  # Returns:
  #   A vector of collective influence for each vertex of the graph
  #   corresponding to the order of vertices output by V(g).
  
  assert_that(d > 0)

  ci <- vector(mode="integer", length=length(V(g)))  # collective influence output
  
  reducedDegreeAll <- degree(g) - 1

  # Only identify nodes at distance d
  nodesAtD <- neighborhood(g, order=d, mindist=d)
  
  for (i in 1:length(nodesAtD)) {
    reducedDegree <- reducedDegreeAll[i]  # i is the index of the node
    rdList <- reducedDegreeAll[as_ids(nodesAtD[[i]])]
    # Setting 0 as default in case the graph doesn't have 2nd-order neighbours
    reducedDegreeNeighbours <- ifelse(length(rdList) > 0, sum(rdList), 0)
    ci[i] <- reducedDegree * reducedDegreeNeighbours
  }

  return(ci)
}

collectiveInfluenceParallel <- function(g, d=3, cores=detectCores()) {
  # Collective Influence as described by Morone & Makse (2015).
  # In simple terms, it is the product of the reduced degree (degree - 1) of a
  # node and the total reduced degree of all nodes at a distance d from the
  # node.
  # This version of the algorithm is parallelized, which can improve processing
  # time for CI on large networks when multiple processing cores are available.
  #
  # Args:
  #   g: the igraph graph
  #   d: the distance, expressed in # of steps, from a given node (default=3)
  #      Distance must be > 0.
  #      According to Morone & Makse, optimal results can be reached at d=3,4,
  #      but this depends on the size/"radius" of the network
  #      NOTE: the distance d is not inclusive. This means that nodes at a
  #      distance of 3 from our node-of-interest do not include nodes at
  #      distances 1 and 2. Only 3.
  #   cores: the number of processing cores to use (default is the system max)
  #
  # Returns:
  #   A vector of collective influence for each vertex of the graph
  #   corresponding to the order of vertices output by V(g).
  
  assert_that(d > 0)
  
  # Store in case multiple V(g) calls takes more time
  allNodes <- V(g)
  
  reducedDegreeAll <- degree(g) - 1
  
  ciParallel <- function(i) {
    # Only identify nodes at distance d
    nodesAtD <- neighborhood(g, nodes=allNodes[i], order=d, mindist=d)
    reducedDegree <- reducedDegreeAll[i]
    rdList <- reducedDegreeAll[as_ids(nodesAtD[[1]])]
    # Setting 0 as default in case the graph doesn't have 2nd-order neighbours
    reducedDegreeNeighbours <- ifelse(length(rdList) > 0, sum(rdList), 0)
    return(reducedDegree * reducedDegreeNeighbours)
  }
  
  ci <- mclapply(1:length(allNodes), ciParallel, mc.cores=cores)
  
  return(unlist(ci))
}

approxLargestEigenvalue <- function(g, d, ci=NULL) {
  # Based on correspondence with Flaviano Morone, the giant component of the
  # network is destroyed when the largest eigenvalue of the non-backtracking
  # matrix is equal to 1 (as proved in their paper). This function approximates
  # the largest eigenvalue.
  #
  # Args:
  #   g: the igraph graph
  #   d: the distance, expressed in # of steps, from a given node. Distance must
  #      be > 0. This is passed to collectiveInfluence().
  #   ci: the precomputed collective influence for the graph; optional
  # 
  # Returns:
  #   An approximation of the largest eigenvalue.
  
  if (is.null(ci)) {
    ci <- collectiveInfluence(g, d)
  } 
  k <- degree(g)
  lambda <- (mean(ci) / mean(k)) ^ (1 / (d + 1))
  return(lambda)
}

getInfluencers <- function(g, d=3, verbose=FALSE, plot=FALSE, layout=FALSE,
                           cores=1) {
  # Runs the collective influence algorithm by successively removing the 
  # maximal influencer from the graph until the giant component has been
  # destroyed.
  #
  # Args:
  #   g: the igraph graph. Nodes of the graph must be named. Names must be
  #      unique.
  #   d: the distance, expressed in # of steps, from a given node (default=3)
  #      Distance must be > 0.
  #      According to Morone & Makse, optimal results can be reached at d=3,4,
  #      but this depends on the size/"radius" of the network
  #      NOTE: the distance d is not inclusive. This means that nodes at a
  #      distance of 3 from our node-of-interest do not include nodes at
  #      distances 1 and 2. Only 3.
  #   verbose: print status messages for each iteration of CI (default=FALSE)
  #   plot: plot each graph after removal of the influencer
  #   layout: track the layout of the graph to preserve node positions even as
  #       influencers are being removed from the network
  #   cores: the number of processing cores to use (default = 1)
  #
  # Returns:
  #   A list containing 2 objects:
  #   1) The updated graph without the influencers
  #   2) A complete list of influencers that, once removed, destroy the giant
  #   component of the network. Returned as a character vector of node names.
  
  assert_that(!is.null(V(g)$name))
  assert_that(!any(duplicated(V(g)$name)))
  assert_that(d > 0)
  assert_that(cores >= 1)
  
  # Used to keep nodes in the same spot for each subsequent plot
  if (layout | plot) {
    fixedLayout <- layout_with_kk(g)
    if (plot) {
      plotg(g, layout=fixedLayout)
    }
  }
  
  if(cores == 1) {
    ci <- collectiveInfluence(g, d)
  } else {
    ci <- collectiveInfluenceParallel(g, d, cores)
  }
  
  ev <- approxLargestEigenvalue(g, d, ci)
  influencers <- vector("character")
  i <- 1
  
  # Keep removing the maximal CI nodes until the eigenvalue = 1. For practical
  # implementation purposes, we continue removing influencers until ev <= 1.
  while (ev > 1) {
    # update this for each iteration since nodes have been removed
    nodes <- V(g)$name
    
    # If there is a tie for collective influence, just choose the first node
    maxCIIndex <- which(ci == max(ci))[1]
    toBeRemoved <- nodes[maxCIIndex]
    influencers[i] <- toBeRemoved
    
    # New network
    g <- delete_vertices(g, toBeRemoved)
    
    # remove coordinates for the influencer node that was deleted
    if (layout | plot) {
      fixedLayout <- fixedLayout[-maxCIIndex, ]
      # Plot each step
      if (plot) {
        plotg(g, layout=fixedLayout)
      }
    }
    
    # Recompute the CI and largest eigenvalue approximation
    if(cores == 1) {
      ci <- collectiveInfluence(g, d)
    } else {
      ci <- collectiveInfluenceParallel(g, d, cores)
    }
    
    ev <- approxLargestEigenvalue(g, d, ci)

    # Optional verbose message
    if (verbose) {
      print(paste("Influencer=", toBeRemoved, "Eigenvalue=", ev, sep=" "))
    }
        
    i <- i + 1
  }
  
  out <- list(graph=g, influencers=influencers)
  if (layout) {
    out$fixedLayout <- fixedLayout
  }
  
  return(out)
}