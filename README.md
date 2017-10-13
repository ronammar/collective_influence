# Collective Influence

Implementation of the collective influence algorithm as presented by [*Morone and Makse* (2015)](http://www.ncbi.nlm.nih.gov/pubmed/26131931). The collective influence computation step has been parallelized for use with large networks.

Depends on [*igraph*](http://igraph.org/).

### Quick sourcing of the script

If you want to test this out quickly:
```r
source("collective_influence_algorithm.R")
plotg(exampleGraph(), community=T)
# find influencers, remove them, and plot the results after each removal
g <- getInfluencers(exampleGraph(), d=2, plot=T)
```

Note that `getInfluencers()` is parallelized , which can improve processing time for collective influence on large networks when multiple processing cores are available. Simply specify the number of cores to use.

### More information

Example from [*Kovacs and Barabasi* News & Views article (2015)](http://www.ncbi.nlm.nih.gov/pubmed/26245576):
(A) the initial network with highest collective influence node in red, highest degree node in yellow. (B) resultant network, with giant component intact, after removal of the 6 nodes with the highest degree. (C) resultant network after removal of top 4 influencers.
![](img/Kovacs and Barabasi graph.jpg)

The CI algorithm at work:

![](img/collective_influence_animation.gif)
