---
layout: page
title: Vignette
subtitle: Last update November 2022
toc: true
toc_label: "Contents"
toc_icon: "cog"
---


In this vignette, we analyze one of the networks provided in the repository. We assume that you
cloned the repository in your computer and that the scripts have permissions to be executed  (see [Install](../Install) page). 
All results obtained can be found in the directory `vignette_example`.


# Input data format

In the directory `data` we have a set of bipartite synthetic matrices used in the publication of the method to test its performance. The networks are labelled with the fields "Nest" and "Conn" indicating the nestedness and connectance of the matrix connecting both pools of nodes, and the label "CompConn" indicates the connectance within the pools. Each pool may represent a pool of species such as plants and their pollinator species (labelled animals) with links representing competition within the pools and mutualistic links between the pools.

We start setting the root directory of the repo as working directory.

```
$> cd path_to_the_repository
```

We are going to work with the network ```Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt```. We have a look at the format:

```
$> less data/Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt
#SpeciesA	SpeciesB	interaction
Pl_1	An_1	1
Pl_2	An_1	1
Pl_3	An_1	1
Pl_4	An_1	1
Pl_5	An_1	1
Pl_6	An_1	1
Pl_7	An_1	1
Pl_8	An_1	1
Pl_11	An_1	1
(...) long file, press q to exit
```

This is a tab-separated file in which the header must start with "#". We can see that in the column  `interaction` we have a value equal to 1 for the interactions between plants and animals, and a negative value for the interactions between pools (animal-animal or plant-plant). Therefore, this column represents the weight of the interactions. There are interactions absent, and should **not** be included with (e.g.) a weight value equal to zero. Since competitive and mutualistic interactions are qualitatively different, we want the algorithm to differentiate them as different types of links. One possibility is to add one column to define the different types (in the next section we show an alternative if you have only positive and negative links). In this file, we simply used a number which is equal to 1 for mutualistic links and to 2 for competitive links in the column "type":

```
$> less data/Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.format2.txt
#SpeciesA	SpeciesB	interaction	type
Pl_1	An_1	1	1
Pl_2	An_1	1	1
Pl_3	An_1	1	1
Pl_4	An_1	1	1
Pl_5	An_1	1	1
Pl_6	An_1	1	1
Pl_7	An_1	1	1
Pl_8	An_1	1	1
Pl_11	An_1	1	1
(...) long file, press q to exit
```

Note that these formats are hard-coded, i.e. the order in which the specific columns should be presented is fixed: nodeA, nodeB, weight, type.

There are more complicated situations depending on whether the network is directed, etc. that can be formatted in different ways, please see the help page of the script [NodeSimilarity.pl](../help) for more details on the format.

# Lazy pipeline

In the folder `scripts/analysisR` we provide a wrapper function (`run_pipeline.R`) that automatically performs the whole pipeline, allowing the user to specify most parameters. In addition, the script `nodeLinkage_pipeline.R` provides an example on how to use it. This is the same example discussed in the "Detailed pipeline" below. We recommend using this wrapper once the pipeline is understood. In particular, an explanation of the output files is provided in the section "Detailed pipeline".


# Detailed pipeline

Here, we describe the different steps that the pipeline includes and how to execute the scripts and interpret the data step by step.

## Similarity between nodes (`NodeSimilarity.pl`)
We now start the search of communities by computing the similarity between nodes, considering the two formats discussed above. We start with the second format, in which we explicitly have a column for the type of interaction.

Remember that the following commands have the root directory of the repo as working directory:

```
$> cd path_to_the_repository
```

Then, we run:

```
./NodeSimilarity.pl -w 1 -d 0 -t 1 -f data/Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.format2.txt
```

In the options we indicated that the file has weights (`-w 1`), it is not directed (`-d 0`), and has different types of links (`-t 1`). See the help page of the script [NodeSimilarity.pl](../help) for more details.
The algorithm prints some information about the number of nodes, links, etc. and it finally returns the file `Nodes-Similarities_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.format2.txt`. Before inspecting this file, let's see how we can run the file for the first format, in which we do not have the column "type" but we have positive and negative links, and we want to interpret the interactions with different signs as different types. This situation is so frequent in ecological networks that we implemented a specific option `-w 2`, in which the algorithm will interpret positive and negative values as different types, with no need of a "type" column:

```
./NodeSimilarity.pl -w 2 -d 0 -t 1 -f data/Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt
```

Note that `-t` must still be given and fixed to one. 

Both files provide the same results (possibly ordered differently). The output looks like this:

```
# >>1NodeA, 2NodeB, 3TanimotoCoeff, 4JaccardCoeff, 5SharedNeighs, 6NeighsA, 7NeighsB
Pl_14   Pl_3    0.0811608092753452      0.0769230769230769      2       7       19
Pl_14   An_46   0       0       0       7       5
Pl_14   An_17   0       0       0       7       9
Pl_14   An_37   0       0       0       7       7
Pl_14   Pl_21   0.0781017924746104      0.0769230769230769      1       7       6
```

where we have, the identities of nodeA and nodeB in the first two columns, the Tanimoto and Jaccard coefficients in columns 3 and 4, the number of shared neighbours in column 5, and the number of neighbours of nodeA and nodeB in columns 6 and 7.

## Clustering nodes (`NodeLinkage.pl`)

The next step is clustering nodes using any of the similarity metrics we computed. In a first run, we will cluster nodes until we have all of them clustered in a single cluster. This will allow us to compute the _partition density_ metrics, which we will use to determine the optimal partition, i.e. the optimal number of communities (clusters). The algorithm requires the original network and the similarity matrix we just computed as inputs, we run:

```
./NodeLinkage.pl -fn data/Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt -fs Nodes-Similarities_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt
```

with no further options.  See the help page of the script [NodeLinkage.pl](../help) for more details of the options that the algorithm offers. 

We obtain two files. The first file `HistExtend-NL_Average_NoStop_$label` is a detailed description of the clustering. The name of the file tell us that the default method of clustering was used (Average Linkage) with no stopping criteria. $label stands for the name of the network. We inspect the file:

```
head  HistExtend-NL_Average_NoStop_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt
#CODE4VALUES_1Step, 2Similarity, 3nodeA, 4nodeB, 5NumIntNodesA, 6NumIntNodesB, 7NumExtNodesA, 8NumExtNodesB, 9NumIntEdgesA, 10NumIntEdgesB, 11NumExtEdgesA, 12NumExtEdgesB
VALUES	1	0.387207196911377	An_29	An_16	0	0	10	9	0	0	10	9
NODES_A	An_29
NODES_B	An_16
EDGES_A	An_9XXXXAn_29	An_16XXXXAn_29	An_17XXXXAn_29	An_18XXXXAn_29	An_21XXXXAn_29	An_28XXXXAn_29	Pl_26XXXXAn_29	An_29XXXXAn_35	An_29XXXXAn_39	An_29XXXXAn_40
EDGES_B	Pl_1XXXXAn_16	Pl_10XXXXAn_16	An_16XXXXAn_17	An_16XXXXAn_21	An_16XXXXAn_28	An_16XXXXAn_29	An_16XXXXAn_35	An_16XXXXAn_40	An_16XXXXAn_45
NODES_AB	An_29	An_16
EDGES_AB	An_9XXXXAn_29	An_16XXXXAn_45	Pl_1XXXXAn_16	An_16XXXXAn_28	Pl_26XXXXAn_29	An_16XXXXAn_29	An_17XXXXAn_29	Pl_10XXXXAn_16	An_21XXXXAn_29	An_16XXXXAn_40	An_16XXXXAn_35	An_29XXXXAn_40	An_29XXXXAn_35	An_29XXXXAn_39	An_16XXXXAn_21	An_28XXXXAn_29	An_18XXXXAn_29	An_16XXXXAn_17
```

There are different types of lines we can easily parse:

* **VALUES** This line appears at the beginning of each clustering step, and it give us a summary of what happened at that step. The fields included are:
    *  `1Step`: Step of the clustering algorithm.
    *  `2Similarity`: Similarity value in which the clusters are joined.
    *  `3nodeA, 4nodeB`: Identities of the clusters joined. After joining, the identity of nodeA will be given to the new cluster.
    *  `5NumIntNodesA, 6NumIntNodesB`: Number of elements within cluster A (B).
    *  `7NumExtNodesA, 8NumExtNodesB`: Number of neighbours (other clusters) connected with cluster A (B).
    *  `9NumIntEdgesA, 10NumIntEdgesB`: Number of edges that cluster A (B) have within its members.
    *  `11NumExtEdgesA, 12NumExtEdgesB`: Number of edges that cluster A (B) have with other clusters.
* **NODES_A, NODES_B, NODES_AB**: Identity of the elements within cluster A, (B or AB). Where AB is the new cluster.
* **EDGES_A, EDGES_B, EDGES_AB**: Links between the elements in cluster A (B or AB) and other clusters. A link (edge) has the format Source_NodeXXXXTarget_Node.

This format will be repeated every time two clusters are joined. The second file `HistCompact-NL_Average_NoStop_$label` is a more compact description of the clustering, which basically includes the same quantities present in the line VALUES of the extended description, plus the _partition density_ values in the columns _Density_, _DensityInt_ and _DensityExt_, representing the total partition density, internal partition density, and external partition density. These are the quantities we will inspect to estimate the optimal partition.

## Identifying the optimal partition (`extractPartDensity.R`).

To analyze the partition density, it is provided in the directory `scripts/analysisR` a function (`extractPartDensity.R`) that extracts the value of the maxima for the three types of partition density quantities, the step in which these maxima are located, and optionally returns a plot with the values at each time step. To illustrate its usage, we also provide the script `nodeLinkage_analysis.R`. Using this method in our example, we obtain the following results:

```
 "-- The maximum value of the total partition density is 0.2251 found at step = 47"
 "-- The maximum value of the internal partition density is 0.1152 found at step = 92"
 "-- The maximum value of the external partition density is 0.1885 found at step = 41"
```

indicating that the optimal partition is found at step 47, where the total partition density peaks. We also observe that the external partition density peaks nearby, and that its maximum is clearly higher than the internal partition density. This suggests that, in this community, guilds are more relevant than modules. Indeed, the plot that we retrieve shows that the internal partition density does not have an important contribution except in the last steps. This makes sense, because the main modules in this network are plants and pollinators, which will be retrieved when there are only two clusters. 

![ Partition Densities](../../_images/Plot_PartitionDensityVsStep.png)


## Obtaining communities

We finally run the clustering with additional arguments to stop it at the desired point. We will select the step = 47, where the total partition density peaks:


```
./NodeLinkage.pl -fn data/Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt -fs Nodes-Similarities_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt -s step -v 47
```

where the flag `-s` indicates that we are using the "step" as stopping criteria, and `-v` the value of the step to stop (47). We obtain four more files. The history files (`HistExt` and `HistCompact`) are identical to those obtained in the first run, except that only contain until step 47. Then we have a new file describing the content of each cluster:

```
$> head Clusters-NL_Average_StopStep-47_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt
# Clusters at Step/Threshold: 47 0.20356097028712
CLUS_1_INFO NumberNodes=        6       NumberEdges=    45
CLUS_1_NODES An_40      An_39   An_21   An_17   An_29   An_16
CLUS_1_EDGES Pl_26XXXXAn_29     An_29XXXXAn_39  An_18XXXXAn_29  An_17XXXXAn_39  An_36XXXXAn_39  Pl_1XXXXAn_39   An_38XXXXAn_40  An_13XXXXAn_21  An_13XXXXAn_17  Pl_1XXXXAn_21   An_40XXXXAn_44  An_16XXXXAn_35  An_33XXXXAn_40  An_28XXXXAn_40  An_8XXXXAn_40   An_40XXXXAn_46  An_17XXXXAn_40  An_17XXXXAn_29  An_16XXXXAn_45  An_21XXXXAn_31  An_16XXXXAn_40  Pl_1XXXXAn_40   An_14XXXXAn_21  An_21XXXXAn_39  Pl_1XXXXAn_17   An_39XXXXAn_44  Pl_1XXXXAn_16   An_28XXXXAn_29  Pl_6XXXXAn_17   An_9XXXXAn_29   An_29XXXXAn_40  An_21XXXXAn_29  An_19XXXXAn_39  An_16XXXXAn_29  An_17XXXXAn_41  Pl_10XXXXAn_16  An_16XXXXAn_17  An_15XXXXAn_40  An_16XXXXAn_28  An_17XXXXAn_37  An_32XXXXAn_40  An_4XXXXAn_40   An_21XXXXAn_28  An_16XXXXAn_21  An_29XXXXAn_35
```

Each cluster is described by three lines:

* **CLUS_$ID_INFO** Tell us the number of nodes and both internal and external edges in the cluster with id $ID.
* **CLUS_$ID_NODES** Identity of the nodes within the cluster.
* **CLUS_$ID_EDGES** Links of the nodes in the cluster connecting them or connecting nodes in external clusters.

The last file the algorithm will generate simply describes the cluster each node belongs:

```
$> head Partition-NL_Average_StopStep-47_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt
# Partitions at Step/Threshold: 47 0.20356097028712
An_40   1
An_39   1
An_21   1
An_17   1
An_29   1
An_16   1
Pl_38   2
An_20   3
An_19   3

```

This file will allow us to separate the nodes in their communities when we represent a network. To see how see the page [Visualization](../Visualization)
