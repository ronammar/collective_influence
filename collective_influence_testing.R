# Implementing Collective Influence algorithm/calculator in R

library(assertthat)
library(igraph)

# Convenience function for plotting graphs
plotg <- function(g) {
  plot(g, vertex.size=20, vertex.color="darkolivegreen3", vertex.label.cex=0.8)
}

# Graph from Kovacs and Barabasi "News & Views" image
g <- graph(
  c(1,4, 2,4, 3,4, 3,5, 4,5, 4,7, 5,6, 5,29, 6,7, 6,10, 6,20, 7,8, 7,9, 7,10, 10,11, 10,30, 11,12, 11,17, 12,13, 12,14, 12,15, 12,16, 14,15, 16,17, 17,18, 17,29, 18,19, 18,20, 20,21, 21,22, 21,23, 21,24, 21,27, 21,30, 22,23, 24,25, 24,26, 24,27, 24,29, 26,27, 27,28, 27,29),
  directed=FALSE
)

# Give the nodes names. This is useful because as we remove nodes, the numbering
# of nodes will change but the names will remain the same.
g <- set_vertex_attr(g, "name", value=paste("n", as.vector(V(g)), sep=""))

plotg(g)

collective_influence <- function(g, d=3) {
  # Collective Influcence as described by Morone & Makse (2015).
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


d <- data.frame(node=V(g)$name, PageRank=page_rank(g)$vector, degree=degree(g),
                stringsAsFactors=FALSE)

# Using PageRank centrality or degree, we identify the same top 6 nodes from the
# figure of Kovacs and Barabasi.

# Now we identify the nodes with the highest collective influence
# See Morone & Makse for formal definition of CI.

d$CI <- collective_influence(g, 2)

View(d)


# Now it's time to identify the influencers and destroy the giant component of
# the network, breaking the network apart into smaller non-interconnected 
# subgraphs.

approxLargestEigenvalue <- function(g, d) {
  # Based on correspondence with Flaviano Morone, the giant component of the
  # network is destroyed when the largest eigenvalue of the non-backtracking
  # matrix is equal to 1 (as proved in their paper). This function approximates
  # the largest eigenvalue.
  #
  # Args:
  #   g: the igraph graph
  #   d: the distance, expressed in # of steps, from a given node. Distance must
  #      be > 0. This is passed to collective_influence().
  # 
  # Returns:
  #   An approximation of the largest eigenvalue.
  
  ci <- collective_influence(g, d)
  k <- degree(g)
  lambda <- (mean(ci) / mean(k)) ^ (1 / (d + 1))
  return(lambda)
}

get_influencers <- function(g, d=3) {
  # Runs the collective influence algorithm by successively removing the 
  # maximal influencer from the graph until the giant component has been
  # destroyed.
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
  #   A list containing 2 objects:
  #   1) 
  #   2) A complete list of influencers that, once removed, destroy the giant
  #   component of the network. Returned as a character vector of node names.
 
  ci <- collective_influence(g, d)
  ev <- approxLargestEigenvalue(g, d)
  influencers <- vector("character")
  i <- 1
  
  # Keep removing the maximal CI nodes until the eigenvalue = 1. For practical
  # implementation purposes, we continue removing influencers until ev <= 1.
  while (ev > 1) {
    # update this for each iteration since nodes have been removed
    nodes <- V(g)$name
    
    # If there is a tie for collective influence, just choose the first node
    toBeRemoved <- nodes[which(ci == max(ci))][1]
    influencers[i] <- toBeRemoved
    
    # New network
    g <- delete_vertices(g, toBeRemoved)
    
    # Recompute the CI and largest eigenvalue approximation
    ci <- collective_influence(g, d)
    ev <- approxLargestEigenvalue(g, d)
    
    i <- i + 1
  }
  
  return(list(graph=g, influencers=influencers))
}

out <- get_influencers(g, 2)

influencers <- out[["influencers"]]
g <- out[["graph"]]

print(paste("Final largest eigenvalue =", ev, collapse=" "))
print(paste("The influencers are:", paste(influencers, collapse=" "), collapse=" "))
plotg(g)