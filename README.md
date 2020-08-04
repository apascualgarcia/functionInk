
# README for functionInk project

For details of the method please see: 

**Pascual-García, A. & Bell, T.** "functionInk: An efficient method to detect functional groups in multidimensional networks reveals the hidden structure of ecological communities". _Methods in Ecology and Evolution_ (2020)


## Description

The aim of this project is to detect communities in complex networks for
which we have not only a weight and/or direction in the edges, but also
a qualitative label describing the kind of links, i.e. what is known as 
multidimensional network.

The strategy is to find a partition of the network such that members
of the same cluster share the same neighbours. We say that two nodes A and B
share a neighbour C if both are connected with C with an edge, and the
direction and type of the edges AC and BC is of the same type.

## Pipeline

The detection of communities in a network considers several steps. One of the reasons why the method is not fully automatic, is that the analysis of the three partition densities is an important step in which the user will get a feeling on the topology of her network, and she will decide which is(are) the relevant type(s) of communities that wants to extract. The strategy considers these steps:

* Compute from your network a similarity measure between the nodes. This similarity will be  higher if the nodes share the same neighbours with the same type of links. This is performed with the algorithm ```NodeSimilarity.pl```, to see the valid formats, options and output files please use the flag ```-h```.

* With the similarity matrix obtained and again considering the original network, you can cluster the nodes with the algorithm ```NodeLinkage.pl```. Please use the flag ```-h``` to see the format, options and output files. This method should be typically used following three steps:
    * In the first run, you will cluster the nodes until no significant relations are found. This should be done to obtain the different measures of the _partition density_ which will be used to determine the optimal stopping point of the clustering.
    * Then analyse the partition densities. This analysis can be done with the script ```nodeLinkage_analysis.R```, from which you will obtain the step at which the total, external and internal partition densities have a maxima and a summary graphic. In our experience, the maxima of the internal partition density brings a partition more similar to the one that would be found with traditional methods maximizing the modularity, while the external partition density would bring you communities that we interpreted in _Pascual-García & Bell_ as _guilds_.
    * After identifying the desired stopping point, you can run again ```NodeLinkage.pl``` with a flag that indicates the step or similarity threshold in which you want to stop the algorithm, and you will obtain at that point the partition of the network and the description of the communities.

## Data

* In the folder `fake_example` a network is provided to test the scripts. Results and figures are also provided.
* In the folder `data` you will find a set of bipartite synthetic matrices used in the publication of the method to test its performance. The networks are labelled with the fields "Nest" and "Conn" indicating the nestedness and connectance of the matrix connecting both pools of nodes, and the label "CompConn" indicates the connectance within the pools.

## Install

The scripts do not require any installation, but the computer should have a Perl interpreter. Most Unix distributions come with a Perl interpreter, if you are running a different OS you can find more information [in this page](https://perldoc.perl.org/5.32.0/perlfaq2.html#What-machines-support-Perl%3f-Where-do-I-get-it%3f). Similarly, the script for the analysis requires R.



