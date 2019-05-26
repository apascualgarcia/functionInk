
# README for functionInk project

## Description

The aim of this project is to detect communities in complex networks for
which we have not only a weight and/or direction in the edges, but also
a qualitative label describing the kind of links, i.e. what is known as 
multidimensional network.

The strategy is to find a partition of the network such that members
of the same cluster share the same neighbours. We say that two nodes A and B
share a neighbour C if both are connected with C with an edge, and the
direction and type of the edges AC and BC is of the same type.

For more information about the method, see _paper_.

## Pipeline

The strategy consists of several steps. 

* Compute from your network a similarity measure between the nodes which will be 
higher if the nodes share the same neighbours with the same type of links. This is performed with the algorithm ```NodeSimilarity.pl```,
to see the valid formats, options and output files please use the flag ```-h```.

* With the similarity matrix obtained and the original network you can cluster
the nodes with the algorithm ```NodeLinkage.pl```. Again use the flag ```-h``` to see
the format, options and output files. This method should be typically run twice:
 * In the first run, you will cluster the nodes until no significant relations are found. This
should be done to obtain the different measures of the _partition density_ which will be used
to determine the optimal stopping point of the clustering.
 * The analysis of the partition density can be done with the script ```nodeLinkage_analysis.R```, from
which you will obtain the step at which the total, external and internal partition densities have a maxima and a summary graphic. In our
experience the maxima of the internal partition density brings a partition more similar to the one that would be
found with traditional methods maximizing the modularity, while the external partition density would bring you
communities that we interpreted in _paper_ as _guilds_.
 * After identifying the desired stopping point, you can run again ```NodeLinkage.pl``` with a flag
that indicates the step or similarity threshold were you want it to stop, and you will obtain at that point
the partition of the network and the description of the communities.



