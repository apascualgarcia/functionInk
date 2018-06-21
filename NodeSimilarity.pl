#!/usr/bin/perl -w
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
# If your network is not weighted, the weight will be fixed to 1
# and Tanimoto coefficients will be computed accordingly. You should, indicate
# this in the "parameters" section of the code. Finally, if your network
# is not weighted but you consider two types of nodes, you should
# provide as "Weight" two values (+1 and +1) and tell to the
# script that the network is weighted.
#
# INPUT: A TAB-separated file describing an undirected network of the format:
#        NodeA   NodeB    Weight 
#        where "Weight" may have sign, what will distinguish
#        between the different types of links. In future versions we
#        may add a field "Type" to consider further possibilities. It accepts
#        header lines starting with '#'.
#        
# OUTPUT: A file describing a similarity matrix of the format:
#         NodeA   NodeB   TanimotoCoeff  JaccardCoeff
#
# USAGE: ./NodeSimilarity $path2network
#        In addition, there are some parameters you should control, 
#        see the first section of the code.
########################################
#TO BE DONE: Generalizar el script para considerar distintas "clases" a la hora de hacer la similaridad en NodeSimilarity.pl. 
#Habría que añadir un campo de lectura (una columna más) que se utilice si quiere el usuario y, si no quiere, que use los signos de la interacción en el clustering.
#########################################
# Silwood Park (Imperial College London)
# July 4th, 2016. Alberto Pascual-García 
# alberto.pascual.garcia@gmail.com 
#########################################
#
use POSIX;

print "  \n";
print "**************************************************  \n";
print "* Building nodes similarities NodeSimilarity.pl  *  \n";
print "**************************************************  \n";
print "  \n";
#
# --- Fix parameters here

$Weighted=0; # One if there is a weighthed network, zero otherwise.
$fieldNodeA=0; # Indicate the column where the first source node is found (minus one)
$fieldNodeB=1; # where the target node is found (minus one)
$fieldWeight=2; # and their interaction value (again minus one). If "Weighted" is zero it will check still if there
 # is a sign, and it will assign a +1 or -1 value

# --- Print the information collected

&printParameters($fieldEdgeA,$fieldEdgeB,$Weighted,$ARGV[0]);

# --- Read the network and build hashes

foreach $line(@INTMP){ #  For each line "nodeA nodeB weight"
    #print join(' ',' ... Reading: ',$line),"\n";  # DEBUG
    if((substr($line,0,1) eq '#')||(substr($line,1,3)eq 'row')){ # skip header and col names
	print join(' ','..Skip header: ',chomp($line)),"\n";
	next;
    }
    chomp($line);
    @fields=split("\t",$line);
    $nodeA=$fields[$fieldNodeA];
    $nodeB=$fields[$fieldNodeB];
    $weight=$fields[$fieldWeight];
    if($Weighted==0){	
	if($weight >0){
	    $weight=1;
	}else{
	    $weight=-1;
	}	  
    }
    $nodes2key{$nodeA}=1; # Store the nodes
    $nodes2key{$nodeB}=1;    
    $edgeTmp=$nodeA.'KKK'.$nodeB; # Create a single identifier for the edge
    push(@{$edge2node{$edgeTmp}},$nodeA); # Relate the edge to their nodes, it may be done simply recovering 
    push(@{$edge2node{$edgeTmp}},$nodeB); # them from the new edge identifier, but this will be faster
    $edge2weight{$nodeA}{$nodeB}=$weight; # This structure basically codify every line
    $edge2weight{$nodeB}{$nodeA}=$weight;
    push(@{$neighbour{$nodeA}},$nodeB); # Store for every node its neighbours
    push(@{$neighbour{$nodeB}},$nodeA);    
} # End foreach reading file

# --- Build lists and control some numbers

@Edges = keys%edge2node;
$Nedges = $#Edges;
@Nodes = keys%nodes2key;
$Nnodes = $#Nodes;

print "  \n";
print join(' ','~~~ The number of edges is: ',$Nedges+1),"\n";
print join(' ','~~~ The number of nodes is: ',$Nnodes+1),"\n";
print "  \n";

# --- For Tanimoto coefficients, we need to compute some measures in advance.
# For every node_i, we define a vector a_i=(A_i1,...,A_ii,...,A_iN)
# where A_ii=mean(A_ij) (i != j)  and A_ij=W_ij. We need to compute from this vector
# the term A_ii and the term abs(a_i*a_i), where * stands for scalar product.

foreach $nodeTmp(@Nodes){ # For each node
    $AvWij=0;
    $Av2Wij=0;
    @neighs=@{$neighbour{$nodeTmp}}; # This is an excess of cleaningness...
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
	@neighsA = @{$neighbour{$nodeA}}; # Again, too clean? I do not like Perl brackets :-)
	@neighsB = @{$neighbour{$nodeB}};
	$selfctrl=0;
	foreach $neighTmpA(@neighsA){ # Note that now each node is a neighbour of itself
	    	foreach $neighTmpB(@neighsB){
		    if($neighTmpA ne $neighTmpB){
			next;
		    }else{ # we still need to control that the link between both nodes if it exist is		       
			if(($neighTmpA eq $nodeA)||($neighTmpA eq $nodeB)){ # counted just once
			    if($selfctrl==0){
				$selfctrl=1;
			    }else{
				next;
			    }
			}
		    }
		    if(defined($edge2weight{$nodeA}{$neighTmpA})){ # With the previous condition this
			 if(defined($edge2weight{$nodeB}{$neighTmpB})){ # is not needed. Just to ctrl everything is ok
			     $Wac=$edge2weight{$nodeA}{$neighTmpA};
			     $Wbc=$edge2weight{$nodeB}{$neighTmpB};
			     $Wprod=$Wac*$Wbc;
			     #print join(' ',' *** Control: i, j, Wac, Wbc',$i,$j,$Wac,$Wbc),"\n"; # DEBUG
			     if($Wprod>0){ # Both weights should have the same sign
				 $Wab+=$Wprod;
				 $Nab+=1;
			     }
			 }
		    }	    
		} # End foreach neighbour $j
	} # End foreach neighbour $i
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
	print OUT0 join(" ",$nodeA,$nodeB,$Tanimoto,$Jaccard,$Nab,$NneighsA,$NneighsB),"\n";
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
#      printParameters
######################
# Print the different input parameters and choices to the standard output

sub printParameters{
    my ($fieldEdgeA,$fieldEdgeB,$Weighted,$pathIn)=@_;
    print '~~~ Reading Node A from column ',$fieldNodeA+1,"\n";
    print '~~~ Reading Node B from column ',$fieldNodeB+1,"\n";
    if($Weighted==0){
	print '~~~ Working with an unweighted network -- Jaccard Similarity',"\n";
    }else{
	print '~~~ Working with a weighted network -- Tanimoto coefficient',"\n";
	print '~~~ Reading weights from column ',$fieldWeight+1,"\n";
    }
    print '  ',"\n";
    
    $pathIn=$ARGV[0];
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
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $year = 1900 + $yearOffset;
    $theTime = " (mm,dd,yy) $month, $dayOfMonth, $year, and time: $hour:$minute:$second,";
    return $theTime;
}
