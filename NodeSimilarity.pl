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
# and Tanimoto coefficients will be computed accordingly. You should indicate
# this in the "parameters" section of the code. Finally, if your network
# is not weighted but you consider two types of nodes, you should
# provide as "Weight" two values (+1 and +1) and tell to the
# script that the network is weighted.
#
# INPUT: A TAB-separated file describing an undirected network of the format:
#        NodeA   NodeB    Weight  Type
#        "Weight" is the strength of the edge (can be positive or negative, if
#        negative the absolute value will be taken in the computation of the Tanimoto coefficients).
#        "Type" is an integer determining the type of edge (e.g. mutualistic=0, competitive=1).  
#        It accepts an indefinite number of header lines starting with '#'.
#        If the network has no weights, then the file should be formatted simply as:
#        NodeA   NodeB   Type
#
# OPTIONS: -h  Prints a help message
#          -W  equal to 0 if the network is not weighted, to 1 otherwise.
# OUTPUT: A file describing a similarity matrix of the format:
#         NodeA   NodeB   TanimotoCoeff  JaccardCoeff
#
# USAGE: ./NodeSimilarity -W types $path2network
#
#        In addition, if you want to change the order of the input columns you can control it 
#        in the first section "parameters".
########################################
#
#########################################
# Silwood Park (Imperial College London)
# July 4th, 2016. Alberto Pascual-García 
# alberto.pascual.garcia@gmail.com
#
# Updated in June, 2018.
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
# ---  parameters here


&readParameters(@ARGV);


# --- Print the information collected

&printParameters($fieldNodeA,$fieldNodeB,$fieldType,$Weighted,$fileIn);

# --- Read the network and build hashes
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
    $type=$fields[$fieldType];
    if($ctrl==0){
	if(!$type || !$fields[$fieldWeight]){print "~ Problems with file format, check type and weight fields. Exit...\n";exit;}		     
    }
    if($Weighted==0){	
	$weight=1;	  
    }else{
	$weight=$fields[$fieldWeight];	    
    }
    $nodes2key{$nodeA}=1; # Store the nodes
    $nodes2key{$nodeB}=1;    
    $edgeTmp=$nodeA.'KKK'.$nodeB; # Create a single identifier for the edge
    push(@{$edge2node{$edgeTmp}},$nodeA); # Relate the edge to their nodes, it may be done simply recovering 
    push(@{$edge2node{$edgeTmp}},$nodeB); # them from the new edge identifier, but this will be faster
    $edge2weight{$nodeA}{$nodeB}=$weight; # This structure basically codify every line
    $edge2weight{$nodeB}{$nodeA}=$weight;
    $edge2type{$nodeA}{$nodeB}=$type; # This structure basically codify every line
    $edge2type{$nodeB}{$nodeA}=$type;
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
		    if(($neighTmpA eq $nodeA)||($neighTmpA eq $nodeB)){ # count it just once
			if($selfctrl==0){
			    $selfctrl=1;
			}else{
			    next;
			}
		    }
		}
		# In any other situation it is a neighbour, we just need to control that the type is the same
		if($edge2weight{$nodeA}{$neighTmpA}){ # This is implicit in the loop above and
		    if($edge2weight{$nodeB}{$neighTmpB}){ # is not needed. Just to ctrl everything is ok
			$pass=0; # control variable to count the shared neighbours
			if(($nodeA eq $neighTmpA)||($nodeB eq $neighTmpB)){ # If it is the link between them we want to count it
			    $pass=1;
			}elsif($edge2type{$nodeA}{$neighTmpA} == $edge2type{$nodeB}{$neighTmpB}){ # also f they belong to the same type
			    $pass=1;
			}
			if($pass == 1){ # Only those that pass the above filters
			    $Wac=$edge2weight{$nodeA}{$neighTmpA};
			    $Wbc=$edge2weight{$nodeB}{$neighTmpB};
			    $Wprod=abs($Wac)*abs($Wbc);
			    $Wab+=$Wprod;
			    $Nab+=1;				
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
#      readParameters
######################
# Print the different input parameters and choices to the standard output

sub readParameters{
    if($ARGV[0] eq "-h"){
	&helpme();
    }elsif($ARGV[0] eq "-W"){
	if(!looks_like_number($ARGV[1])){
	    print " \n";
	    print ">> Input argument error for  NodeSimilarity.pl: \n";
	    print $ARGV[1]," is not a valid Weight flag \n";
	    print ">> run  ./NodeSimilarity.pl -h for help \n";
	    print " \n";
	    exit;
	}else{
	    $Weighted=$ARGV[1];
	}
	if(!$ARGV[2]){
	    print " \n";
	    print ">> Input argument error for NodeSimilarity.pl: \n";
	    print ">> input file name missing \n";
	    exit;
	}else{	
	    $fileIn=$ARGV[2];
	}
    }else{
	if(looks_like_number($ARGV[0])){
	    print " \n";
	    print ">> Input argument error for NodeSimilarity.pl: \n";
	    print $ARGV[0]," is not a valid file name \n";
	    print ">> run  ./NodeSimilarity.pl -h for help \n";
	    print " \n";
	    exit;
	}else{
	    $fileIn=$ARGV[0];
	}
	if($ARGV[1] ne "-W"){
	    print " \n";
	    print ">> Input argument error for  NodeSimilarity.pl: \n";
	    print " -W flag expected instead of $ARGV[1] \n";
	    print ">> run  ./NodeSimilarity.pl -h for help \n";
	    print " \n";
	    exit;
	}else{
	    if(!looks_like_number($ARGV[2])){
		print " \n";
		print ">> Input argument error for  NodeSimilarity.pl: \n";
		print $ARGV[2]," is not a valid Weight flag \n";
		print ">> run  ./NodeSimilarity.pl -h for help \n";
		print " \n";
		exit;
	    }else{
		$Weighted=$ARGV[2];
	    }
	}
    }
    if($Weighted == 0){
	$fieldNodeA=0; # Indicate the column where the first source node is found (minus one)
	$fieldNodeB=1; # where the target node is found (minus one)
	$fieldType=2;
    }else{
	$fieldNodeA=0; # Indicate the column where the first source node is found (minus one)
	$fieldNodeB=1; # where the target node is found (minus one)
	$fieldWeight=2; # and their interaction value (again minus one). 
	$fieldType=3;
    }    
    if($Weighted==1){
	return $fileIn, $Weighted, $fieldNodeA, $fieldNodeB,$fieldWeight,$fieldType;
    }else{
	return $fileIn, $Weighted, $fieldNodeA, $fieldNodeB,$fieldType;  
    }
}

######################
#      printParameters
######################
# Print the different input parameters and choices to the standard output

sub printParameters{
    my ($fieldNodeA,$fieldNodeB,$fieldType,$Weighted,$fileIn)=@_;
    print '~~~ Reading Node A from column ',$fieldNodeA+1,"\n";
    print '~~~ Reading Node B from column ',$fieldNodeB+1,"\n";
    if($Weighted==0){
	print '~~~ Working with an unweighted network -- Jaccard Similarity',"\n";
    }else{
	print '~~~ Working with a weighted network -- Tanimoto coefficient',"\n";
	print '~~~ Reading weights from column ',$fieldWeight+1,"\n";
	print '~~~ Reading types from column ',$fieldType+1,"\n";
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
# Return the time where the script is being executed

sub helpme{
    print "\n";
    print " > Help for NodeSimilarity.pl\n";
    print " > USAGE: ./NodeSimilarity -W weighted path2network \n";
    print "       - weighted: equal to one if weighted network, zero otherwise \n";
    print "       - path2network: Path to a tab separated file with the following fields for weighted networks \n";
    print "          NodeA   NodeB    Weight  Type \n";
    print "         and, for unweighted networks formatted as: \n";
    print "           NodeA   NodeB   Type     \n";
    print "         where: \n";    
    print "           -- Weight: is the strength of the edge (can be positive or negative, if \n";
    print "                  negative the absolute value will be taken in the computation of the Tanimoto coefficients).\n";
    print "           -- Type: is an integer determining the type of edge (e.g. mutualistic=0, competitive=1). \n"; 
    print "       - It accepts an indefinite number of header lines starting with '#'.\n";
    print "\n";
    print "\n";
    exit;
}
