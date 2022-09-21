---
layout: page
title: Help
subtitle: Last update September 2022
---

# Options scripts

_Return to [HOME](/)_

The following are the options of the two main scripts. These options are printed using the -h flag, for example:

```
$> ./NodeSimilarity.pl -h
```

### Options for NodeSimilarity.pl

Generic call:

```
$> ./NodeSimilarity.pl -w $Option1 -d $Option2 -t $Option4 -f $File
```
Options:

``` 
 FLAGS: -h  Prints a help message 
        -w  equal to 0 if the network is not weighted, to 1 otherwise (required).
        -d  equal to 0 if the network is undirected, to 1 otherwise (required).
        -t  equal to 0 if the links are of the same type, to 1 if they are of different types (required). 
        -f  flag to include the input file (required).

 INPUT: A TAB-separated file describing a network with the format:
        [1] If the network is undirected the general format is: 
                   NodeA   NodeB    Weight  Type
            where: 
            ... "NodeA" is the source node and "NodeB" the target node.
            ... "Weight" is a real value indicating the strength of the link  (can be positive or negative, if
                negative the absolute value will be taken in the computation of the Tanimoto coefficients).
            ... "Type" is a string indicating the type of the link (e.g. mutualistic=0, competitive=1).
            ... The script accepts an indefinite number of header lines starting with # 

        [2] If the network is directed the general format is the same, but it is assumed that 
            the direction is encoded in the order in which the names of the nodes appear, i.e. 
            NodeA   NodeB  vs. NodeB   NodeA. There is no specific assumption on which node is 
            is source or target, since the algorithm just needs to know that the order matters 
            to consider them as different types of links. This means that it is possible to   
            encode the presence of directed links using the field Type, explained in section [4]

        [3] If the network has both directed and undirected links, then the flag -d 1 should be
            used (i.e. as if it would be directed) and those nodes linked with
            an undirected link should appear twice in both directions and with the same weight:
                            NodeA    NodeB    Weight   Type
                            NodeB    NodeA    Weight   Type
            An alternative possibility to encode this situation is explained in section [4]

        [4] An alternative to encode directed links is to consider each direction 
            (or the lack of direction if there are also undirected links) as an attribute for the field Type.
            If each link then has an additional qualitative attribute, it should additionally be 
            considered. For instance, consider you have directed and undirected links with an 
            attribute that can be White or Black: 
                       NodeA NodeB Weight White   (undirected)   
                       NodeB NodeC Weight White   (directed)   
                       NodeA NodeC Weight Black   (directed)   
            we could transform the field type into a format like this:   
                       NodeA NodeB Weight UndirWhite     
                       NodeB NodeC Weight DirWhite      
                       NodeA NodeC Weight DirBlack      
            and then we use the options needed for an undirected network (section [1]) 

        [5] If the network has no weights, either you use the general format (with -w 1) and all the weights are equal to one,
            or you use -w 0, and then the file can be simply formatted as:
                            NodeA   NodeB   Type
        [6] If the network has no types,  either you use the general format (with -t 1)  and all your types are the same,
            or if you use -t 0 then the file can be simply formatted as:
                            NodeA   NodeB   Weight
        [7] If the network has no types and no weights, either you use the general format (with -t 1 and -w 1) ,
           and all your types are the same and weights equal to one or you use -t 0 and -w 0,
           in which case the file can be simply formatted as:
                            NodeA   NodeB  

 OUTPUT: A file describing a similarity matrix of the format:
         NodeA   NodeB   TanimotoCoeff  JaccardCoeff

 EXAMPLE USAGE: ./NodeSimilarity -w 1 -d 1 -t 1 -f path2network

 COMMENTS: In addition, if you want to change the order of the input columns you can code it 
        in the function "readParameters".
```


### Options for Linkage.pl

Generic call:

```
$>  ./NodeLinkage -fs path2SimilarityMatrix -fn path2OriginalNetwork -s $option1 -v $option2 -a $option3 -c $option4 
```
Options:

``` 
  INPUT: All the inputs, required and optional, are introduced with a flag: 
  - Required arguments 
  
       -fs path_to_file 
           A tab-separated long-formatted-matrix with the all-against all topological similarity between nodes, as computed by  
           NodeSimilarity.pl, with the format:     
                  NodeA NodeB  Similarity1(A,B) Similarity2(A,B) .... 
                  NodeA NodeC  Similarity1(A,C) Similarity2(A,C) .... 
                  .... 
  
       -fn  path_to_file 
          A tab-separated network from which the above similarity matrix was derived, i.e. the input of NodeSimilarity.pl, with the format: 
                  NodeA NodeB   Weight 
          If the name of the file is "Network"-label, "label" will be used for 
          the name of the output, otherwise you will find default names. 
  
        Note: Both input files accept a header starting with the character # 
  
   - Optional arguments: 
  
      -h 
          Prints this help and exits 
  
      -c integer 
          An integer with the column in which the similarity measure between nodes will be found, it 
          is given as an option because the ouptut of NodeSimilarity.pl provides Tanimoto and Jaccard 
          coefficients in different columns. Defaults to column 3 (Tanimoto), column 4 is Jaccard. 
  
      -a method 
          Where method determines the clustering method. Valid arguments are "Average" for Average linkage, 
          "Single" for single linkage and  "Complete" for complete linkage. Default is Average. 
  
      -s stop_criteria 
          Where stop_criteria is a string determining a criteria to stop the clustering. It may be a threshold in  
          the similarity (argument "thres") a stopping point ("step"), or we may want to cluster until 
          there is a single cluster, in which case the argument should be "none"). Default value is "none" 
  
      -v value 
          If a stop criteria is given, this flag must be used to include either a value that would 
          be either a threshold (for the similarity criteria) or a clustering step (for the step criteria). 
  
   OUTPUT:  
      A list of files: 
      - If it is not given any no stop criteria: 
           HistExtend.NoStop.InputFile 
           HistCompact.NoStop.InputFile 
      - With a stop criteria: 
           HistExtend.StopCriteria.InputFile 
           HistCompact.StopCriteria.InputFile 
           Clusters.StopCriteria.InputFile  
           Partition.StopCriteria.InputFile 
      -where: 
        -- HistExtend file: explicitly describes the clustering process, showing the two clusters 
           that are clustered (its nodes, edges, etc.) 
        -- HistCompact file: it only provides relevant quantities of the two clusters joined (number of edges, 
           number of elements contained), and of the clustering, (step, partition densities, etc.)	
        -- Clusters file: Describes the different clusters at the stopping point 
        -- Partition file: A vector assigning to each node the cluster id it belongs to. 
  EXAMPLE USAGE:  
   ./NodeLinkage.pl -fs path2SimilarityMatrix -fn path2OriginalNetwork -s step -v 145 -a single -c 4 
  
```

_Return to [HOME]()_
