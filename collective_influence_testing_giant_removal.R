# This test script is old; it relied on sourcing of the collective_infuence()
# function from another script.

# Graph from Kovacs and Barabasi "News & Views" image
g <- graph(
  c(1,4, 2,4, 3,4, 3,5, 4,5, 4,7, 5,6, 5,29, 6,7, 6,10, 6,20, 7,8, 7,9, 7,10, 10,11, 10,30, 11,12, 11,17, 12,13, 12,14, 12,15, 12,16, 14,15, 16,17, 17,18, 17,29, 18,19, 18,20, 20,21, 21,22, 21,23, 21,24, 21,27, 21,30, 22,23, 24,25, 24,26, 24,27, 24,29, 26,27, 27,28, 27,29),
  directed=FALSE
)
# Give the nodes names. This is useful because as we remove nodes, the numbering
# of nodes will change but the names will remain the same.
g <- set_vertex_attr(g, "name", value=paste("n", as.vector(V(g)), sep=""))

plotg <- function(g) {
  plot(g, vertex.size=20, vertex.color="darkolivegreen3", vertex.label.cex=0.8)
}

# Based on correspondence with Flaviano Morone, the largest eigenvalue is equal
# to 1 when the giant component is destroyed. This eigenvalue is calculated
# below.
approxLargestEigenvalue <- function(g, d) {
  ci <- collective_influence(g, d)
  k <- degree(g)
  lambda <- (mean(ci) / mean(k)) ^ (1 / (d + 1))
  return(lambda)
}

# plot initial network
plotg(g)

ci <- collective_influence(g, 2)
ev <- approxLargestEigenvalue(g, 2)
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
  ci <- collective_influence(g, 2)
  ev <- approxLargestEigenvalue(g, 2)
  
  i <- i + 1
}

print(paste("Final largest eigenvalue =", ev, collapse=" "))
print(paste("The influencers are:", paste(influencers, collapse=" "), collapse=" "))
plotg(g)
