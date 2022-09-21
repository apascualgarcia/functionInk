---
layout: page
title: "-functionink- Functional groups linkage"
subtitle: A method for community detection in multiplex networks.
cover-img: ["/assets/img/NodesSimilarity.png" : "Nodes similarity","/assets/img/BlockModelling_Vs_Amaral.png" : "Comparison of methods","/assets/img/Guilds.png" : "Guilds","/assets/img/Guild_Vs_Module.png" : "Guild Vs. Module","/assets/img/PartitionDensity.png" : "Partition density","/assets/img/Entangled_Bank.png" : "The bank disentangled"]
---

## Description

The aim of this project is to detect communities in complex networks in
which we have not only a weight and/or direction in the edges, but also
a qualitative label describing the type of links. Having different types of links
increases the complexity in the detection of communities, because our definition of
community may require that members of the same community have links of the same type only.

Considering types of links also facilitates to naturally encode _multiplex networks_, in which nodes are distributed
in different layers. Links within layers will be encoded with one type of link for each layer, and
links between layers by simply adding a new type for each set of links connecting a pair of layers. This leads
to a very flexible way to find communities in otherwise very complicated networks.

The strategy is to find a partition of the network such that members
of the same commmunity (approximately) share the same neighbours with connections of the same type and similar weights. Therefore, we say that two nodes A and B
share a neighbour C if both are connected with C with an edge, and the
direction and type of the edges AC and BC is of the same type.


## Overview of the pipeline

The strategy to detect communities in a network considers several steps. For more details please follow the [Vignette](./_pages/Vignette)

* Compute from your network a similarity measure between the nodes. This similarity will be  higher if the nodes share the same neighbours with the same type of links. For this task, we use the script  `NodeSimilarity.pl`.  See the [help](_pages/help) page to explore the options.
* With the similarity matrix obtained we cluster the nodes with an agglomerative algorithm twice:
    * In the first run, we cluster the nodes until no significant relations are found. During this computation, the algorithm will compute different measures of _partition density_ which will be used to estimate the optimal stopping point of the clustering. Depending on the measure used, we can identify an optimal point for the two types of communities that the method identifies, _guilds_ or  _modules_ (see publication below for details). We use the script `NodeLinkage.pl`. See the [help](_pages/help) page to explore the options.
    * Analyse the partition densities. We identify the maxima of these quantities, which determine the clustering stopping point in which we want to obtain the communities. As we mentioned, the stopping point will correspond to one of the maxima in the partition densities. In our experience, the maxima of the internal partition density brings a partition more similar to the one that would be found with traditional methods maximizing the modularity, while the external partition density would bring you communities that we interpreted in _Pascual-García & Bell_ as _guilds_. To extract the maxima please use the R function  `extractPartDensity.R`, whose use is illustrated in the script `nodeLinkage_analysis.R`
    * After identifying the desired stopping point, we run again the agglomerative clustering with an option to indicate the step or similarity threshold in which we want to stop the algorithm. We will obtain at that point the partition of the network (i.e. an id identifying the community each node belongs) and the description of the communities.


## Data

* In the folder `fake_example` a network is provided to test the scripts. Results and figures are also provided. See the [Install](_pages/Install) page for details.
* In the folder `data` you will find a set of bipartite synthetic matrices used in the publication of the method to test its performance. The networks are labelled with the fields "Nest" and "Conn" indicating the nestedness and connectance of the matrix connecting both pools of nodes, and the label "CompConn" indicates the connectance within the pools. See the [Vignette](_pages/Vignette) for details.



## Citation

If you use this method please cite: 

**Pascual-García, A. & Bell, T.**  (2020). functionInk: An efficient method to detect functional groups in multidimensional networks reveals the hidden structure of ecological communities. _Methods in Ecology and Evolution_, 11(7), 804-817.


