#!/usr/bin/env perl 
# ******************************
# * NodeSimilarity.pl             *
# ******************************
#
# This script takes networks obtained from any method and
# it computes two similarity measures between its nodes. 
# A similar method was first proposed in Ahn et al. Nature (2010)
# for clustering links, and we have an implementation in Perl called
# EdgeLinkage.pl, over which this script has been built.
# We consider here a subtle modification of their similarity measure
# to take into account different types of edges. Basically,
# we aim to find communities taking into account that nodes belonging
# to the same community not only share the same partners, but they
# share the same type of relationships between these partners.
# We modify the similarity for both the Jaccard measure (which do
# not consider a weighted network) and Tanimoto coefficients (weighted).
#
# FLAGS: -h  Prints a help message
#        -w  equal to 0 if the network is not weighted, to 1 otherwise.
#        -d  equal to 0 if the network is undirected, to 1 otherwise
#        -f  flag to include the input file
#
# INPUT: A TAB-separated file describing an undirected network of the format:
#                   NodeA   NodeB    Weight  Type
#     -- "Weight" is the strength of the edge (can be positive or negative, if
#        negative the absolute value will be taken in the computation of the Tanimoto coefficients).
#     -- "Type" is an integer determining the type of edge (e.g. mutualistic=0, competitive=1).  
#     --  It accepts an indefinite number of header lines starting with '#'.
#     -- Networks with particular formats: 
#        --- If the network has no weights, then the file should be formatted simply as:
#                            NodeA   NodeB   Type
#        --- If the network is directed NodeA should be the source node and NodeB the target node.
#        --- If the network has both directed and undirected links the flag -d 1 should be
#            included (i.e. as if it would be directed) and those nodes linked with
#            an undirected link should appear twice in both directions and with the same weight:
#                            NodeA    NodeB    Weight   Type
#                            NodeB    NodeA    Weight   Type
#
# OUTPUT: A file describing a similarity matrix of the format:
#         NodeA   NodeB   TanimotoCoeff  JaccardCoeff
#
# USAGE: ./NodeSimilarity.pl -w 1 -d 1 -f $path2network
#
#        In addition, if you want to change the order of the input columns you can control it 
#        in the first section "parameters".
########################################
#
#########################################
# Silwood Park (Imperial College London)
# July 4th, 2016. Alberto Pascual-García 
# apascualgarcia.github.io
# Updated in September, 2022.
# Zürich, Theoretical Biology (ETH)
#########################################
#

use POSIX;
use Scalar::Util qw(looks_like_number);

print "  \n";
print "**************************************************  \n";
print "* Building nodes similarities NodeSimilarity.pl  *  \n";
print "**************************************************  \n";
print "  \n";
#
# ---  Read parameters

&readParameters(@ARGV); 

# --- Print the information collected

#&printParameters($fieldNodeA,$fieldNodeB,$fieldType,$Weighted,$fileIn);
&printParameters();

# --- Read the network and build hashes
print " \n";
print ">> Reading the network: \n";
print "~~~ The first lines for the fields read from file and after conversions are: \n";
$ctrl=-1;
foreach $line(@INTMP){ #  For each line "nodeA nodeB weight type"
    #print join(' ',' ... Reading: ',$line),"\n";  # DEBUG
    if((substr($line,0,1) eq '#')||(substr($line,1,3)eq 'row')){ # skip header and col names
	print join(' ','..Skip header: ',chomp($line)),"\n";
	next;
    }
    $ctrl+=1;
    chomp($line);
    @fields=split("\t",$line);
    $nodeA=$fields[$fieldNodeA];
    $nodeB=$fields[$fieldNodeB];
    #if(!looks_like_number($type)){ # Updated to consider them character
    #    print ">> Link types should be integers, I\'ve found this: $type \n";
    #    &abort();
    #}
    if($Weighted==0){	
	$weight=1;
	if($typed == 0){
	    $type = "1";
	}else{
	    $type=$fields[$fieldType];
	}
    }elsif($typed == 0){
	$type = "1";
	$weight=$fields[$fieldWeight];	    
    }else{
	if($Weighted==1){ # weight and type are provided
	    $weight=$fields[$fieldWeight];
	    $type=$fields[$fieldType];
	}else{ # only weight is provided, but the sign determines the type
	    $weight=$fields[$fieldWeight];
	    if($weight >= 0){		
		$type= "1";
	    }else{		
		$type= "0";
	    }
	}
    }
    if($ctrl < 10){
	print "~~~ Reading fields: nodeA = $nodeA,nodeB= $nodeB, type= $type, weight = $weight, \n";
	
    }
    #if($ctrl==0){
	#if(!$type || !$weight){print "~ Problems with file format, check type and weight fields. Exit...\n";exit;}		     
    #}
    $nodes2key{$nodeA}=1; # Store the nodes
    $nodes2key{$nodeB}=1;    
    $edgeTmp=$nodeA.'KKK'.$nodeB; # Create a single identifier for the edge (not used in this version)
    push(@{$edge2node{$edgeTmp}},$nodeA); # Relate the edge to their nodes, it may be done simply recovering  (not used in this version)
    push(@{$edge2node{$edgeTmp}},$nodeB); # them from the new edge identifier, but this will be faster  (not used in this version)
    if($edge2weight{$nodeA}{$nodeB}){ # If it already exists it means that it is an hybrid link
	if($directed==0){ # Not allowed for undirected networks
	    print ">> Repeated edge in input file: $nodeA $nodeB \n";
	    &abort();
	}
	$edge2dir{$nodeA}{$nodeB}=0; # Only fix the directions to be undirected, both the types and neighbours were already stored 
	$edge2dir{$nodeB}{$nodeA}=0;	
    }else{ # If it does not exist	
	$edge2weight{$nodeA}{$nodeB}=$weight; # Store the weight
	$edge2weight{$nodeB}{$nodeA}=$weight;
	$edge2type{$nodeA}{$nodeB}=$type; # Store the type
	$edge2type{$nodeB}{$nodeA}=$type;
	if($directed == 0){ # If it  is undirected
	    $edge2dir{$nodeA}{$nodeB}=0; # Both have the same directions
	    $edge2dir{$nodeB}{$nodeA}=0;
	}else{ # Otherwise
	    $edge2dir{$nodeA}{$nodeB}=1; # fix the directions
	    $edge2dir{$nodeB}{$nodeA}=0;
	}
	push(@{$neighbour{$nodeA}},$nodeB); # Store for every node its neighbours, below we will also introduce the node
	push(@{$neighbour{$nodeB}},$nodeA); # itself as its own neighbour   
    } # End foreach reading file
}
# --- Build lists and control some numbers

@Edges = keys%edge2node;
$Nedges = $#Edges;
@Nodes = keys%nodes2key;
$Nnodes = $#Nodes;

#print "  \n";
print join(' ','~~~ The number of edges is: ',$Nedges+1),"\n";
print join(' ','~~~ The number of nodes is: ',$Nnodes+1),"\n";
print "  \n";

# --- For Tanimoto coefficients, we need to compute some measures in advance.
# For every node_i, we define a vector a_i=(A_i1,...,A_ii,...,A_iN)
# where A_ii=mean(A_ij) (i != j)  and A_ij=W_ij. We need to compute from this vector
# the term A_ii and the term abs(a_i*a_i), where * stands for scalar product.
print " \n";
print " >> Pre-processing coefficients \n";
foreach $nodeTmp(@Nodes){ # For each node
    $AvWij=0;
    $Av2Wij=0;
    @neighs=@{$neighbour{$nodeTmp}}; # Just for peace of mind
    foreach $neighTmp(@neighs){ # Work with its neighbours
	$Wij=$edge2weight{$nodeTmp}{$neighTmp}; # Recover the weights
	$AvWij+=abs($Wij); # Store to compute the average of weights, no sign here (A_ii)
	$Av2Wij+=$Wij**2; # And the scalar product (a_i*a_i)
    } # end foreach @neighbour
    $Nneighs=$#neighs+1; # DEBUG Check this +1. If there is one neigh, the position is zero, so it is needed because still the node itself is not inside
    $edge2weight{$nodeTmp}{$nodeTmp}=$AvWij/$Nneighs; # Term  A_ii=mean(A_ij) (i != j)
    $Av2Wij+=($AvWij/$Nneighs)**2;
    $node2aa{$nodeTmp}=$Av2Wij; # term abs(a_i*a_i)
    push(@{$neighbour{$nodeTmp}},$nodeTmp); # Finally, introduce the node as its own neighbour
    #print join(' ',' *** Control: mean(A_ij),  abs(a_i*a_i)',$AvWij/$Nneighs,$Av2Wij),"\n"; # DEBUG
} # end foreach @Nodes
print "  \n";

# --- Start computing similarities between nodes and print the first output

# Header for the output -> This should be moved to a function

&theTime(); # Recover the time
$fileOut0='Nodes-Similarities_'.$fileIn;
open(OUT0, ">$fileOut0") || die "Couldn't open file $fileOut0"; # 

print OUT0 '# >> Output from NodeSimilarity.pl <<',"\n";
print OUT0 '# >> Input file: "',$fileIn,"\n";
print OUT0 '# >> Jaccard and Tanimoto similarity between nodes in the network ',"\n"; 
print OUT0 '# >> Running at date: ',$theTime,"\n";
print OUT0 '# >> Silwood Park, Imperial College London, A.P-G.', "\n";
print OUT0 '# >>1NodeA, 2NodeB, 3TanimotoCoeff, 4JaccardCoeff, 5SharedNeighs, 6NeighsA, 7NeighsB',"\n";

# Start computing
print " \n";
print " >> Start main computation  \n";
print '~~~ Computing similarities between nodes...',"\n";

$ctrl=1; # vars to control the timeline of the computation, no worries with this
$frac=floor($Nnodes/4);

for($i=0; $i<$Nnodes; $i++){ # Edge 1
    if($i==$ctrl*$frac){
	print '~~~ I have already processed ',25*$ctrl,'% of the nodes!',"\n";;
	$ctrl+=1;
    }
    $nodeA=$Nodes[$i]; # Recover its nodes   
    for($j=$i+1; $j<=$Nnodes; $j++){ # Edge 2
	$nodeB=$Nodes[$j];
	#... Compute their similarity
	$Wab=0;
	$Nab=0;
	#... First, check if A and B are linked themselves (I think I am computing this twice)
	# if(defined($edge2weight{$nodeA}{$nodeB})){
	#     $Wprod=$edge2weight{$nodeA}{$nodeB}**2;
	#     $Wab+=$Wprod;
	#     $Nab+=1;
	# }
	#... Next, check if they share any third
	@neighsA = @{$neighbour{$nodeA}}; # Again, peace of mind
	@neighsB = @{$neighbour{$nodeB}};
	$selfctrl=0;
	foreach $neighTmpA(@neighsA){ # Note that now each node is a neighbour of itself
	    foreach $neighTmpB(@neighsB){
		if($neighTmpA ne $neighTmpB){ # If they do not share the neighbour skip
		    next;
		}else{ # if it is the same, control if it is a link between both nodes (not with a neighbour) and if it is		       
		    if(($neighTmpA eq $nodeB)||($neighTmpB eq $nodeA)){ # count it just once
			if($selfctrl==0){
			    $selfctrl=1;
			    if($neighTmpA eq $nodeB){ # then, neighTmpB=nodeB so we change it
				$neighTmpB=$nodeA; # to point to nodeA (with selfctrl we will do this just once)
			    }else{ # the other way around if neighTmpB eq nodeA, neighTmpA is pointing to nodeA
				$neighTmpA=$nodeB; # so we make it point to nodeB
			    }
			}else{
			    next;
			}
		    }
		}
		# In any other situation, we just need to control that the type is the same
		if($edge2weight{$nodeA}{$neighTmpA}){ # This is implicit in the loop above and
		    if($edge2weight{$nodeB}{$neighTmpB}){ # is not needed. Just to ctrl everything is ok
			if($edge2type{$nodeA}{$neighTmpA} eq $edge2type{$nodeB}{$neighTmpB}){ # if it is a link between them or, being against a third, have the same type 
			    if($edge2dir{$nodeA}{$neighTmpA} == $edge2dir{$nodeB}{$neighTmpB}){ # and the same direction
				$Wac=$edge2weight{$nodeA}{$neighTmpA};
				$Wbc=$edge2weight{$nodeB}{$neighTmpB};
				$Wprod=abs($Wac)*abs($Wbc);
				$Wab+=$Wprod;
				$Nab+=1;				
			    }
			}
			#if(!$edge2type{$nodeB}{$neighTmpB}){
			#    print join(', ',' *** Control: nodeB',$nodeB,'tmpB',$neighTmpB,'Wbc',$edge2type{$nodeB}{$neighTmpB}),"\n"; # DEBUG
			#    #exit;
			#}
			#print join(', ',' *** Control: i',$i,'j',$j,'Wac',$Wac,'Wbc',$Wbc),"\n"; # DEBUG
		    }
		}	    
	    } # End foreach neighbour $j
	} # End foreach neighbour $i
	#print join(', ',' *** Control: nodeA',$nodeA,'Waa',$node2aa{$nodeA},' nodeB:',$nodeB,' Wbb',$node2aa{$nodeB},'Wab',$Wab),"\n"; # DEBUG
	$NneighsA=$#{$neighbour{$nodeA}};#+1; # As usual, check this +1. Now is not needed, because we have already one more neigh (the node itself)
	$NneighsB=$#{$neighbour{$nodeB}};#+1;	
	$Tanimoto=$Wab/($node2aa{$nodeA}+$node2aa{$nodeB}-$Wab);
	$Jaccard=$Nab/($NneighsA+$NneighsB);
	# if(($nodeC eq '"caulobacteraceae"')&&($Nab > 0)){
	#     print join(' ',' *** Control: neighbours of A',$Nab,$nodeA,$#{$neighbour{$nodeA}},@neighsA),"\n"; # DEBUG
	#     print join(' ',' *** Control: neighbours of B',$nodeB,$#{$neighbour{$nodeB}},@neighsB),"\n"; # DEBUG
	#     exit;
	# }
	#print join(' ',' *** Control: $NneighsA,$NneighsB,$node2aa{$nodeA}, $Wab,$Nab',$NneighsA,$NneighsB,$node2aa{$nodeA}, $Wab,$Nab),"\n"; # DEBUG
	print OUT0 join("\t",$nodeA,$nodeB,$Tanimoto,$Jaccard,$Nab,$NneighsA,$NneighsB),"\n";
    } # End for - similarity computation $j
} # End for - similarity computation $i

print '~~~ Done!',"\n";
print '  ',"\n";
print '****************************',"\n";
print '** Program finished',"\n";
print '** Check your results  ',"\n";
print '** Bye! ',"\n";
print '****************************',"\n";
print '  ',"\n";
print '  ',"\n";



###################################################
##                                       FUNCTIONS
###################################################
######################
#      readParameters
######################
# Print the different input parameters and choices to the standard output

sub readParameters{
    
    $messageOk{"-w"}="~~~ The network is weighted=1 or 2/unweighted=0? Value = ";
    $messageOk{"-d"}="~~~ The network is directed=1/undirected=0? Value = ";
    $messageOk{"-t"}="~~~ The network has different types=1/or a single type=0? Value = ";
    $messageOk{"-f"}="~~~ The network file is = ";
    $messageErr{"-w"}="~~~ The value for network weight argument (-w) is not numeric = ";
    $messageErr{"-d"}="~~~ The value for network directionality argument (-d) is not numeric = ";
    $messageErr{"-t"}="~~~ The value for the link types argument (-t) is not numeric = ";;
    $messageErr{"-f"}="~~~ This is not a valid name for a file = ";  
    $typeArg{"-d"}=1;
    $typeArg{"-w"}=1;
    $typeArg{"-t"}=1;
    $typeArg{"-f"}="string";

    
    if($ARGV[0] eq "-h"){ # If help is needed exit
	&helpme();	
    }
    $Nargs=$#ARGV;
    if($Nargs != 7){ # If there are missing arguments exit
	print " \n";
	print ">> Missing input arguments: \n";
	print ">> run  ./NodeSimilarity.pl -h for help \n";
	print " \n";
	exit;	
    }else{ # Otherwise
	print ">> Reading input arguments: \n";
	for($i=0; $i<$Nargs; $i=$i+2){ # For each flag 
	    if(looks_like_number($ARGV[$i+1]) eq looks_like_number($typeArg{$ARGV[$i]})){ # Check if the value passed for the argument makes sense
		$loadVar{$ARGV[$i]}=$ARGV[$i+1]; # load the value
		#print " \n";
		print $messageOk{$ARGV[$i]},$loadVar{$ARGV[$i]},"\n";
		#print " \n";
	    }else{
		print " \n";
		print $messageErr{$ARGV[$i]},$ARGV[$i+1],"\n";
		print ">> run  ./NodeSimilarity.pl -h for help \n";
		print ">> Leaving the script... \n";
		print " \n";
		exit;
	    }
	}
    }
    $Weighted=$loadVar{"-w"};
    $directed=$loadVar{"-d"};
    $typed=$loadVar{"-t"};
    $fileIn=$loadVar{"-f"};
    # --- Assign columns to different file types depending on the input file, and return values
    # ..... Independently of the options the first two columns are always the same
    $fieldNodeA=0; # Indicate the column where the first source node is found (minus one)
    $fieldNodeB=1; # where the target node is found (minus one)
    if($Weighted == 0){
	if($typed == 0){ # No weight no type
	    return $fileIn, $Weighted,$typed, $fieldNodeA, $fieldNodeB;
	}else{ # No weight only, but there is type
	    $fieldType=2; # where the type is found
	    return $fileIn, $Weighted, $typed, $fieldNodeA, $fieldNodeB,$fieldType;  
	}
    }elsif($typed == 0){ # No weight only
	$fieldWeight=2; # field for the interaction value (again minus one).
	return $fileIn, $Weighted, $typed, $fieldNodeA, $fieldNodeB,$fieldWeight;  
    }else{
	$fieldWeight=2; # and their interaction value (again minus one). 
	$fieldType=3; # where the type is found
	return $fileIn, $Weighted, $fieldNodeA, $fieldNodeB,$fieldWeight,$fieldType;
    }    
}

######################
#      printParameters
######################
# Print the different input parameters and choices to the standard output

sub printParameters{
    #my ($fieldNodeA,$fieldNodeB,$fieldType,$Weighted,$fileIn)=@_;
    print " \n";
    print ">> Processing input arguments: \n";
    print '~~~ Reading Node A from column ',$fieldNodeA+1,"\n";
    print '~~~ Reading Node B from column ',$fieldNodeB+1,"\n";
    if($Weighted==0){
	print '~~~ Working with an unweighted network -- Jaccard Similarity',"\n";
	if($typed==0){
	    print '~~~ Working with a network with a single type of link',"\n"; 
	}else{
	    print '~~~ Working with a network with different types of links',"\n"; 
	    print '~~~ Reading types from column ',$fieldType+1,"\n";
	}
    }elsif($typed == 0){
	print '~~~ Working with a weighted network -- Tanimoto coefficient',"\n";
	print '~~~ Reading weights from column ',$fieldWeight+1,"\n";
	print '~~~ Working with a network with a single type of link',"\n"; 
    }else{
	print '~~~ Working with a weighted network -- Tanimoto coefficient',"\n";
	print '~~~ Working with a network with different types of links',"\n";
	if($Weighted==1){
	    print '~~~ Reading weights from column ',$fieldWeight+1,"\n";
	    print '~~~ Reading types from column ',$fieldType+1,"\n";
	}else{
	    print '~~~ Reading weights from column ',$fieldWeight+1,"\n";
	    print '~~~ Inferring types from the sign of the weights',"\n";
	}
    }
    print '  ',"\n";
    
    $pathIn=$fileIn;
    print join(' ','~~~ Input path: ',$pathIn),"\n";
    open(INTMP,$pathIn)  or die "Can't find file $pathIn\n";
    @fields=split("/",$pathIn);
    $fileIn=$fields[$#fields]; # Take the last field of the path as the name of the file for outputs
    print join(' ','~~~ Input file: ',$fileIn),"\n";
    print '  ',"\n";
    
    @INTMP = <INTMP>;
    close(INTMP);
    return @INTMP,$pathIn,$fileIn;    
}

######################
#     theTime
######################
# Return the time where the script is being executed

sub theTime{
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset) = localtime();
    $year = 1900 + $yearOffset;
    $theTime = " (mm,dd,yy) $month, $dayOfMonth, $year, and time: $hour:$minute:$second,";
    return $theTime;
}

######################
#     helpme
######################
# Return a help message

sub helpme{
    print "\n";
    print " >> Help for NodeSimilarity.pl\n";
    print "\n";
    print " FLAGS: -h  Prints a help message \n";
    print "        -w  equal to 0 if the network is not weighted, to 1 otherwise (required). If the value is equal to 2 \n";
    print "            the sign of the weight will be used to create two types (type = A if weight >=0 and type = B otherwise).  \n";
    print "            The flag -t should be fixed to 1, since types exist, but no column should be provided.  \n";
    print "        -d  equal to 0 if the network is undirected, to 1 otherwise (required).\n";
    print "        -t  equal to 0 if the links are of the same type, to 1 if they are of different types (required). \n";
    print "        -f  flag to include the input file (required).\n";
    print "\n";
    print " INPUT: A TAB-separated file describing a network with the format:\n";
    print "        [1] If the network is undirected the general format is: \n";
    print "                   NodeA   NodeB    Weight  Type\n";
    print "            where: \n";               
    print "            ... \"NodeA\" is the source node and \"NodeB\" the target node.\n";
    print "            ... \"Weight\" is a real value indicating the strength of the link  (can be positive or negative, if\n";
    print "                negative the absolute value will be taken in the computation of the Tanimoto coefficients).\n";
    print "            ... \"Type\" is a string indicating the type of the link (e.g. mutualistic=0, competitive=1).\n";
    print "            ... The script accepts an indefinite number of header lines starting with # \n";
    print "\n";
    print "        [2] If the network is directed the general format is the same, but it is assumed that \n";
    print "            the direction is encoded in the order in which the names of the nodes appear, i.e. \n";
    print "            NodeA   NodeB  vs. NodeB   NodeA. There is no specific assumption on which node is \n";
    print "            is source or target, since the algorithm just needs to know that the order matters \n";
    print "            to consider them as different types of links. This means that it is possible to   \n";
    print "            encode the presence of directed links using the field Type, explained in section [4]\n";
    print "\n";
    print "        [3] If the network has both directed and undirected links, then the flag -d 1 should be\n";
    print "            used (i.e. as if it would be directed) and those nodes linked with\n";
    print "            an undirected link should appear twice in both directions and with the same weight:\n";
    print "                            NodeA    NodeB    Weight   Type\n";
    print "                            NodeB    NodeA    Weight   Type\n";
    print "            An alternative possibility to encode this situation is explained in section [4]\n";
    print "\n";
    print "        [4] An alternative to encode directed links is to consider each direction \n";
    print "            (or the lack of direction if there are also undirected links) as an attribute for the field Type.\n";
    print "            If each link then has an additional qualitative attribute, it should additionally be \n";
    print "            considered. For instance, consider you have directed and undirected links with an \n";
    print "            attribute that can be White or Black: \n";
    print "                       NodeA NodeB Weight White   (undirected)   \n";
    print "                       NodeB NodeC Weight White   (directed)   \n";
    print "                       NodeA NodeC Weight Black   (directed)   \n";
    print "            we could transform the field type into a format like this:   \n";
    print "                       NodeA NodeB Weight UndirWhite     \n";
    print "                       NodeB NodeC Weight DirWhite      \n";
    print "                       NodeA NodeC Weight DirBlack      \n";
    print "            and then we use the options needed for an undirected network (section [1]) \n";
    print "\n";
    print "        [5] If the network has no weights, either you use the general format (with -w 1) and all the weights are equal to one,\n";
    print "            or you use -w 0, and then the file can be simply formatted as:\n";
    print "                            NodeA   NodeB   Type\n";
    print "        [6] If the network has no types,  either you use the general format (with -t 1)  and all your types are the same,\n";
    print "            or if you use -t 0 then the file can be simply formatted as:\n";
    print "                            NodeA   NodeB   Weight\n";
    print "        [7] If the network has no types and no weights, either you use the general format (with -t 1 and -w 1) ,\n";
    print "           and all your types are the same and weights equal to one or you use -t 0 and -w 0,\n";
    print "           in which case the file can be simply formatted as:\n";
    print "                            NodeA   NodeB  \n";	
    print "\n";
    print " OUTPUT: A file describing a similarity matrix of the format:\n";
    print "         NodeA   NodeB   TanimotoCoeff  JaccardCoeff\n";
    print "\n";
    print " USAGE: ./NodeSimilarity.pl -w 1 -d 1 -t 1 -f path2network\n";
    print "\n";
    print " COMMENTS: In addition, if you want to change the order of the input columns you can code it \n";
    print "        in the function \"readParameters\".\n";
    print "\n";
    print "\n";
    exit;
}


######################
#    Abort
######################
# Abort execution

sub abort{
    print '~~~ I abort the execution...',"\n";
    print '~~~ Exit!',"\n";
    print '  ',"\n";
    exit
}
