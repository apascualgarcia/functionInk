#!/usr/bin/perl -w
# ******************************
# * NodeLinkage.pl             *
# ******************************
#
# This script takes a similarity matrix between nodes in a 
# network,similar to the network proposed to cluster links by
# Ahn et al. Nature (2010), which we implement in EdgeLinkage.pl, and it
# performs an agglomerative linkage clustering computing additionally
# the partition density of the clustering, to find out an
# objective stopping point. The computation of this density
# is the main reason why we develop new code and we do not
# use previous implementations of classical agglomerative clusterings.
# The script is able to compute three types of clustering Single, Complete
# and Average Linkage. If not threshold or clustering step to stop is provided,
# it will return as output a matrix with the step and threshold at which
# the two elements where joined into the same cluster "Matrix2dendogram". This file is useful to
# build a dendogram. Another output is a "History"
# file printing the clusters that were joined at every step and the Partition
# Density. If not threshold or clustering step to stop is provided,
# it will return as output  two history files describing the clustering,
# printing the clusters that were joined at every step and the Partition
# Density. If a threshold or clustering step to stop is provided, it will print
# the different clusters "Clusters" and a matrix where, for every pair, it is indicated
# with a number the cluster they belong to "Matrix2clusters". This file is the one will be used
# to represent the network colouring the edges with the partition found.
# The idea is to run the algorithm with no stopping point and then analyse the Partition 
# Density looking for a maximum, which would be the stopping point. You will re-run
# then the algorithm.
#
# INPUT: A similarity matrix (all against all) between nodes obtained from 
#        the network topology (NodeSimilarity.pl), with the format:
#        "NodeA" "NodeB"   Similarity 
#        The network from which the above matrix was derived, with the format:
#        "NodeA" "NodeB"   Weight
#        If the name of the file is "Network"-$label, $label will be used for
#        the name of the output, otherwise you will find simple names.
# OUTPUT: "HistExtend"."NoStop".$InputFile
#         "HistCompact"."NoStop".$InputFile
#  - With stopping point:
#         "HistExtend"."StopCriteria".$InputFile
#         "HistCompact"."StopCriteria".$InputFile
#         "Matrix2Clusters_EdgeLinkage"."StopCriteria".$InputFile 
#         "Clusters_EdgeLinkage"."StopCriteria".$InputFile
#
# USAGE: ./NodeLinkage $path2SimilarityMatrix $path2OriginalNetwork
#        In addition, there are some parameters you should control, 
#        see the first section of the code.
#
#########################################
# Silwood Park (Imperial College London)
# July 4th, 2016. Alberto Pascual-García 
# alberto.pascual.garcia@gmail.com 
#########################################
#
# --- Modules used
use POSIX;
use List::Util qw( min max );
# use strict;
# use warnings;

# --- Fix parameters here:

#... clustering type
$LinkCode=2; # 0=Single Linkage; 1=Complete Linkage; 2=Average Linkage
$Linkage[0]="Single"; $Linkage[1]="Complete"; $Linkage[2]="Average";
$StopStep=129; # Zero to cluster until the last step, an integer to stop at some step
$StopThr=0; # Zero to cluster until the last step, a threshold to stop at some step

#... Fields to read from the first file (similarity between nodes due to the network structure)
$fieldNodeSimA=0; # Indicate the column where the first source node is found minus one
$fieldNodeSimB=1; # where the target node is found (minus one)
$fieldSim=2; # and their similarity (again minus one).  

#... Fields to read from the second file (original network)
$fieldNodeNetA=0; # Indicate the column where the first source node is found minus one
$fieldNodeNetB=1; # where the target node is found (minus one)
$fieldWeight=2; # and their similarity (again minus one).  

# --- Print the information collected to the standard output

print "  \n";
print "***********************************************  \n";
print "* Finding communities with nodes clustering   *  \n";
print "***********************************************  \n";
print "  \n";

&printParameters($fieldNodeSimA,$fieldNodeSimB,$fieldSim,$fieldNodeNetA,$fieldNodeNetB,$fieldWeight,$StopStep,$StopThr,$LinkCode,@Linkage);

# --- Open input file

print '~~ Opening the first input file: ',"\n";
&filesInSetUp($ARGV[0]); # It returns the array @INTMP, $fileIn and $labelIn
@INTMP1=@INTMP;
$fileIn1=$fileIn;
$labelIn1=$labelIn;

print '~~ Opening the second input file: ',"\n";
&filesInSetUp($ARGV[1]); # It returns the array @INTMP
@INTMP2=@INTMP;
$fileIn2=$fileIn;
$labelIn2=$labelIn;

# --- Open and build headers for the output files
&filesOutSetUp();

# --- Read the original network. This will be used to compute the partition density

foreach$line(@INTMP2){ #  For each line "nodeA nodeB weight"
    #print join(' ',' ... Reading: ',$line),"\n";  # DEBUG
    if((substr($line,0,1) eq '#')||(substr($line,1,3)eq 'row')){ # skip header and col names
	print join(' ','..Skip header: ',chomp($line)),"\n";
	next;
    }
    chomp($line);
    @fields=split("\t",$line);
    $nodeA=$fields[$fieldNodeNetA];
    $nodeB=$fields[$fieldNodeNetB];
    if(exists($edges2key{"${nodeB}XXXX$nodeA"})){ # It might happen that you have not only the upper triangular matrix
	next;  # but the whole matrix (e.g. if the network comes from SparCC). If the link was processed already we
    }          # will be able to control it here
    $weight=$fields[$fieldWeight];
    $nodes2netWeight{$nodeA}{$nodeB}=$weight;
    $nodes2netWeight{$nodeB}{$nodeA}=$weight;
    $nodes2key{$nodeA}=1; # Store the  nodes
    $nodes2key{$nodeB}=1; 
    if(!exists(${$cluster2nodes{$nodeA}}[0])){ # first time we find this node 
	push(@{$cluster2nodes{$nodeA}},$nodeA); # A single node is a cluster of size one
    }
    if(!exists(${$cluster2nodes{$nodeB}}[0])){ # Same for node B,
	push(@{$cluster2nodes{$nodeB}},$nodeB); # A single node is a cluster of size one
    }
    push(@{$cluster2edges{$nodeA}},"${nodeA}XXXX$nodeB"); 
    push(@{$cluster2edges{$nodeB}},"${nodeA}XXXX$nodeB"); # Keeping the order will allow us to eliminate redundancies later (see COMM1)   
    $edges2key{"${nodeA}XXXX$nodeB"}=1;
}
# --- Build list and control some numbers

@Nodes = keys%nodes2key;
$TotalNodes = $#Nodes+1;
@Edges = keys%edges2key;
$TotalEdges = $#Edges+1;

print "  \n";
print join(' ','~~~ The number of nodes to cluster is: ',$TotalNodes+1),"\n";
print join(' ','~~~ The number of edges in the network is: ',$TotalEdges+1),"\n";


# --- Read the similarity between nodes in the network (Jaccard or Tanimoto) and build hashes, 
#... this will be the measure used to perform the clustering.
#... Identify already the pair with the first maximum to initialize the clustering

foreach $line(@INTMP1){ #  For each line "nodeA nodeB networkSimilarity"
    #print join(' ',' ... Reading: ',$line),"\n";  # DEBUG
    if(substr($line,0,1) eq '#'){ # skip header and col names
	print join(' ','..Skip header: ',chomp($line)),"\n";
	next;
    }
    chomp($line);
    @fields=split(/\s+/,$line);
    $nodeA=$fields[$fieldNodeSimA];
    $nodeB=$fields[$fieldNodeSimB];
    $sim=$fields[$fieldSim];
    $value{$nodeA}{$nodeB}=$sim; # Store their similarity
    $value{$nodeB}{$nodeA}=$sim;
    $value{$nodeA}{$nodeA}=0; # This will work as a control
    $value{$nodeB}{$nodeB}=0;
    if((!defined($maxValue{$nodeA}))||($sim > $maxValue{$nodeA})){ # Recover maximum values for A
	$maxValue{$nodeA}=$sim;
	$maxPartner{$nodeA}=$nodeB;
    }
    if((!defined($maxValue{$nodeB}))||($sim > $maxValue{$nodeB})){ # for B
	$maxValue{$nodeB}=$sim;
	$maxPartner{$nodeB}=$nodeA;
    }
    if((!defined($maxGlobalValue))||($sim > $maxGlobalValue)){ # and for the pair
	$maxGlobalValue=$sim;
	$MaxPair="${nodeA}XXXX$nodeB";
    }

} # End foreach reading file

print join(' ','~~~ The maximum value found is:',$maxGlobalValue,'for pair',$MaxPair),"\n";
print "  \n";

# --- Start the clustering

print '*************************',"\n";
print ' STARTING CLUSTERING...',"\n";
print '*************************',"\n";
print "  \n";
$CtrlClust=0;
$Density=0;
$DensityInt=0;
$DensityExt=0;
$NcumInt=0;
$NcumExt=$TotalEdges;
$Ncum=$TotalEdges;
$Step=0;
while($CtrlClust==0){
    #... Recover the nodes from the pair
    $Step+=1;
    @pair=split("XXXX",$MaxPair);
    $nodeA=$pair[0];
    $nodeB=$pair[1];
    print join(' ','** JOINING: *************************'),"\n"; #DEBUG
    print join(' ',$nodeA,$nodeB),"\n"; #DEBUG
    print join(' ','************************************'),"\n"; #DEBUG
    print " ","\n";
    #... Count the number of nodes and edges before joining 
    #... it is neded for Average Linkage and the partition density
   
    print join(" ",'~~~ Results for A'),"\n"; #DEBUG
    #print join(' ',  @{$cluster2edges{$nodeA}},'KKK',@{$cluster2nodes{$nodeA}}),"\n"; #DEBUG
    &splitIntExtEdges(\@{$cluster2edges{$nodeA}},\@{$cluster2nodes{$nodeA}}); 
    $NumNodesA=$NumNodes;
    $NumEdgesA=$NumEdges;
    $NumIntNodesA=$NumIntNodes;
    $NumIntEdgesA=$NumIntEdges;
    $NumExtNodesA=$NumExtNodes;
    $NumExtEdgesA=$NumExtEdges;
    print join(" ",'~~~ Results for B'),"\n"; #DEBUG
    #print join(' ',  @{$cluster2edges{$nodeB}},'KKK',@{$cluster2nodes{$nodeB}}),"\n"; #DEBUG
    &splitIntExtEdges(\@{$cluster2edges{$nodeB}},\@{$cluster2nodes{$nodeB}});   
    $NumNodesB=$NumNodes;
    $NumEdgesB=$NumEdges;
    $NumIntNodesB=$NumIntNodes;
    $NumIntEdgesB=$NumIntEdges;
    $NumExtNodesB=$NumExtNodes;
    $NumExtEdgesB=$NumExtEdges;

    $maxGlobalTmp=$maxGlobalValue; # Store the maximum, it might be helpful for some measures

    #... Print some info to monitor the clustering
    print join(" ",'~~~ Clustering at Step/Thr ',$Step,$maxGlobalTmp),"\n";
    print join(" ",' ::: Nodes ',$nodeA,$nodeB,'which contain:'),"\n";
    print join(" ",' ::: Number of Nodes A/B:',$NumNodesA,'/',$NumNodesB),"\n";
    print join(" ",' ::: Number of Edges A/B:',$NumEdgesA,'/',$NumEdgesB),"\n";
    print join(" ",' ::: Number of Internal Nodes A/B:',$NumIntNodesA,'/',$NumIntNodesB),"\n";
    print join(" ",' ::: Number of External Nodes A/B:',$NumExtNodesA,'/',$NumExtNodesB),"\n";
    print join(" ",' ::: Number of Internal Edges A/B:',$NumIntEdgesA,'/',$NumIntEdgesB),"\n";
    print join(" ",' ::: Number of External Edges A/B:',$NumExtEdgesA,'/',$NumExtEdgesB),"\n";

    #... Print information of A and B  in the extended history file before joining them
    print HISTextend join(" ",'VALUES',$Step,$maxGlobalTmp,$nodeA,$nodeB,$NumIntNodesA,$NumIntNodesB,$NumExtNodesA,$NumExtNodesB,$NumIntEdgesA,$NumIntEdgesB,$NumExtEdgesA,$NumExtEdgesB),"\n";
    print HISTextend join(" ",'NODES_A',@{$cluster2nodes{$nodeA}}),"\n";
    print HISTextend join(" ",'NODES_B',@{$cluster2nodes{$nodeB}}),"\n";
    print HISTextend join(" ",'EDGES_A',@{$cluster2edges{$nodeA}}),"\n";
    print HISTextend join(" ",'EDGES_B',@{$cluster2edges{$nodeB}}),"\n";

    #... Start Joining the clusters. First, send the nodes from B to A
    push(@{$cluster2nodes{$nodeA}},@{$cluster2nodes{$nodeB}});

    #... Send edges from B to A 
    push(@{$cluster2edges{$nodeA}},@{$cluster2edges{$nodeB}});
    foreach $EdgeTmp(@{$cluster2edges{$nodeA}}){ # COMM1: We need to eliminate redundant interactions here
	$Edge2key{$EdgeTmp}=1;
    }
    @{$cluster2edges{$nodeA}}=keys(%Edge2key);
    undef(%Edge2key);

    # if(exists($nodes2netWeight{$nodeA}{$nodeB})){ # The link will be repeated
    # 	$edgeTmp="${nodeB}XXXX$nodeA"; # COMM1: Here is where we need to differentiate them
    # 	@{$cluster2edges{$nodeA}} = grep { $_ ne $edgeTmp  } @{$cluster2edges{$nodeA}};
    # }

    #... Determine internal and external nodes and edges of the new cluster
    #print join(" ",'~~~ Results for AB'),"\n"; #DEBUG
    &splitIntExtEdges(\@{$cluster2edges{$nodeA}},\@{$cluster2nodes{$nodeA}});    
    $NumNodesAB=$NumNodes;
    $NumEdgesAB=$NumEdges;
    $NumIntNodesAB=$NumIntNodes;
    $NumIntEdgesAB=$NumIntEdges;
    $NumExtNodesAB=$NumExtNodes;
    $NumExtEdgesAB=$NumExtEdges;

    #... Print the information after joining them
    print join(" ",' ::: Number of Nodes A+B:',$NumNodesAB),"\n";
    print join(" ",' ::: Number of Edges A+B:',$NumEdgesAB),"\n";
    print join(" ",' ::: Number of internal Nodes A+B:',$NumIntNodesAB),"\n";
    print join(" ",' ::: Number of internal Edges A+B:',$NumIntEdgesAB),"\n";
    print join(" ",' ::: Number of external Nodes A+B:',$NumExtNodesAB),"\n";
    print join(" ",' ::: Number of external Edges A+B:',$NumExtEdgesAB),"\n";

    #... extended file as well
    print HISTextend join(" ",'NODES_AB',@{$cluster2nodes{$nodeA}}),"\n";
    print HISTextend join(" ",'EDGES_AB',@{$cluster2edges{$nodeA}}),"\n";

    #... Update the relation values 

    @Nodes = grep { $_ ne $nodeB } @Nodes; # first, delete B from the list
    $maxValue{$nodeA}=0;
    foreach $ThirdTmp(@Nodes){ 
	if($ThirdTmp eq $nodeA){
	    next;
	}

	#print '   ~ Third / maxPartner / maxValue ',"\n";
	#print join(' ',' ::: ',$ThirdTmp,$maxPartner{$ThirdTmp},$maxValue{$ThirdTmp}),"\n";	
	#print join(' ',' :::  VALS(A-C)(B-C)',$value{$nodeA}{$ThirdTmp},$value{$nodeB}{$ThirdTmp}),"\n";
	#... update the values according with the different linkage criteria
	$valueA=$value{$nodeA}{$ThirdTmp};
	$valueB=$value{$nodeB}{$ThirdTmp};
	if(($valueA==0)&&($valueB==0)){
	    next;
	}
	if($Linkage[$LinkCode] eq "Single"){ # Update according with the type of clustering
	    $valueTmp=max($valueA,$valueB);
	    #print join(" ",$valueA,$valueB,$valueTmp),"\n"; # DEBUG
	    #print join(" ",$nodeA,$nodeB,$ThirdTmpA),"\n"; # DEBUG
	}elsif($Linkage[$LinkCode] eq "Complete"){
	    $valueTmp=min($valueA,$valueB);
	}else{                             # Average linkage, we weight by the size of the cluster
	    $valueTmp=($NumNodesA*$valueA+$NumNodesB*$valueB)/($NumNodesA+$NumNodesB);
	}
	if($valueTmp > 0){ # This should be just needed for Complete Linkage
	    $thirdListTmp{$ThirdTmp}=1; # We monitor if any third ends with a zero relation
	}
	$value{$nodeA}{$ThirdTmp}=$valueTmp; # Here is where the update actually takes place
	$value{$ThirdTmp}{$nodeA}=$valueTmp;
	$value{$nodeB}{$ThirdTmp}=0; # Fix to zero the relations with B
	$value{$ThirdTmp}{$nodeB}=0;
	# ... look for the new maximum in the new cluster
	if($valueTmp > $maxValue{$nodeA}){ 
	    $maxValue{$nodeA}=$valueTmp;
	    $maxPartner{$nodeA}=$ThirdTmp;
	}
	#... Verify that A and B were not the partners for which Third had a maximum
	#.... Not needed if nodeA AND single linkage or nodeB AND complete, 
	#.... but it is easy to have a mistake so I redo it anyway
	if(($maxPartner{$ThirdTmp} eq $nodeA)||($maxPartner{$ThirdTmp} eq $nodeB)){ # And check if they were maximum partners	  
	    $ValueTmp=$maxValue{$ThirdTmp};
	    $PartTmp=$maxPartner{$ThirdTmp};
	    $maxValue{$ThirdTmp}=0; # look for a new maximum
	    foreach $NewPartner (@Nodes){
		if($NewPartner eq $ThirdTmp){
		    next;
		}
		if($value{$ThirdTmp}{$NewPartner}>$maxValue{$ThirdTmp}){
		    $maxValue{$ThirdTmp}=$value{$ThirdTmp}{$NewPartner};
		    $maxPartner{$ThirdTmp}=$NewPartner;			
		}
	    }
	    #if($maxValue{$ThirdTmp}==0){ # It would be a disjoint cluster, why should I stop? Likely to happen with CompleteLinkage
	    #	print '~~~ I did not find a new maximum for ',$ThirdTmp,"\n";
	    #	&abort();
	    #}	    
	} #endif( Third partner was A or B)

    } #endforeach thirds

    #... look for the new maximum
    $maxGlobalValue=0;
    foreach $nodeTmp(@Nodes){
	#print join(" ",'>>>',$nodeTmp,$maxValue{$nodeTmp}),"\n"; # DEBUG
	if($maxValue{$nodeTmp}>$maxGlobalValue){
	    $maxGlobalValue=$maxValue{$nodeTmp}; #
	    $MaxPair="${nodeTmp}XXXX$maxPartner{$nodeTmp}";
	}
    }

    #... compute the partition density
    
    $Da=&clusterDensity($NumNodesA,$NumEdgesA,$NumIntNodesA,$NumIntEdgesA,$NumExtNodesA,$NumExtEdgesA);
    $Da=$D;
    $DintA=$D1;
    $DextA=$D2;
    $Db=&clusterDensity($NumNodesB,$NumEdgesB,$NumIntNodesB,$NumIntEdgesB,$NumExtNodesB,$NumExtEdgesB);
    $Db=$D;
    $DintB=$D1;
    $DextB=$D2;
    $Dab=&clusterDensity($NumNodesAB,$NumEdgesAB,$NumIntNodesAB,$NumIntEdgesAB,$NumExtNodesAB,$NumExtEdgesAB);
    $Dab=$D;
    $DintAB=$D1;
    $DextAB=$D2;
    $Density=$TotalEdges*$Density; # "Undo" the average
    $Density=$Density-$Da-$Db+$Dab; # Extract the contribution from the old clusters and add up the new one
    $Density=(1/$TotalEdges)*$Density; # "Redo" the average
    $DensityInt=$TotalEdges*$DensityInt; # Same for the internal contribution
    $DensityInt=$DensityInt-$DintA-$DintB+$DintAB; # 	 $D2=($NextEdges-$NextNodes)/($NextNodes*($NintNodes-1));

    $DensityInt=(1/$TotalEdges)*$DensityInt; # 
    $DensityExt=$TotalEdges*$DensityExt; # And external
    $DensityExt=$DensityExt-$DextA-$DextB+$DextAB; # 
    $DensityExt=(1/$TotalEdges)*$DensityExt; #
    $NcumInt=$NcumInt-$NumIntEdgesA-$NumIntEdgesB+$NumIntEdgesAB;
    $NcumExt=$NcumExt-($NumExtEdgesA+$NumExtEdgesB-$NumExtEdgesAB)/2;
    $Ncum=$NcumInt+$NcumExt;
    if($Density <0){
	print '~~~ Total density became negative',"\n";
	print join(" ",'TotalEdges', 'Da','Db','Dab','Density'),"\n";
	print join(" ",$TotalEdges,$Da,$Db,$Dab,$Density),"\n";
	&abort();
    }
    #... Control the clustering stopping
    if($Stop==1){ # We stop if we reach either
	if($StopStep>0){ 
	    if($Step>=$StopStep){ # a certain step
		$CtrlClust=1;
	    }
	}elsif($StopThr>0){
	    if($maxGlobalValue>=$StopThr){ # a certaing threshold
		$CtrlClust=1;
	    }
	}
    }elsif($Step==$TotalEdges){ # or we have a single cluster (last step)
 	$CtrlClust=1;	
     }elsif($maxGlobalValue==0){ # no significant links to join
 	$CtrlClust=1;
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',"\n";
	print "  \n";
	print 'WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~',"\n";
	print '::: Stopping clustering at Step: ',$Step,"  \n";
	print '::: No more significant relations were found ',"  \n";	
     }
  

    #... Print the compact history file, which includes the density
    print HISTcompact join(" ",$Step,$maxGlobalTmp,$Density,$DensityInt,$DensityExt,$NumNodesA,$NumEdgesA,$NumNodesB,$NumEdgesB,$NumNodesAB,$NumEdgesAB,$NumIntNodesA,$NumIntNodesB,$NumExtNodesA,$NumExtNodesB,$NumIntNodesAB,$NumExtNodesAB,$NumIntEdgesA,$NumIntEdgesB,$NumExtEdgesA,$NumExtEdgesB,$NumIntEdgesAB,$NumExtEdgesAB,$nodeA,$nodeB,$NcumInt,$NcumExt,$Ncum),"\n";
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',"\n";
    print "  \n";
} #endwhile clustering

if($Stop==1){ # If there is a classification and partition to print
    print CLU '# Clusters at Step/Threshold: ',$Step,' ',$maxGlobalTmp,"\n";
    print PAR '# Partitions at Step/Threshold: ',$Step,' ',$maxGlobalTmp,"\n";
    $Nclust=0;
    foreach $NodeTmp(@Nodes){
	$Nclust+=1;
	$NumNodes=$#{$cluster2nodes{$NodeTmp}}+1;
	$NumEdges=$#{$cluster2edges{$NodeTmp}}+1;
	print CLU 'CLUS_',$Nclust,'_INFO ',join(" ",'NumberNodes=',$NumNodes,'NumberEdges=',$NumEdges),"\n";
	print CLU 'CLUS_',$Nclust,'_NODES ',join(" ",@{$cluster2nodes{$NodeTmp}}),"\n";	
	print CLU 'CLUS_',$Nclust,'_EDGES ',join(" ",@{$cluster2edges{$NodeTmp}}),"\n";
	@NodesInClus=@{$cluster2nodes{$NodeTmp}};
	foreach $NodeInClusTmp(@NodesInClus){
	    print PAR join("\t",$NodeInClusTmp,$Nclust),"\n";
	}	
    } #endforeach @Nodes total
}


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
    my ($fieldNodeSimA,$fieldNodeSimB,$fieldSim,$fieldNodeNetA,$fieldNodeNetB,$fieldWeight,$StopStep,$StopThr,$LinkCode,@Linkage)=@_;

    print '~~ FIELDS for the FIRST FILE',"\n";
    print '~~~ Reading Node A from column ',$fieldNodeSimA+1,"\n";
    print '~~~ Reading Node B from column ',$fieldNodeSimB+1,"\n";
    print '~~~ Reading Similarities from column ',$fieldSim+1,"\n";
    print '~~ FIELDS for the SECOND FILE',"\n";
    print '~~~ Reading Node A from column ',$fieldNodeNetA+1,"\n";
    print '~~~ Reading Node B from column ',$fieldNodeNetB+1,"\n";
    print '~~~ Reading Similarities from column ',$fieldWeight+1,"\n";
    print '~~ CLUSTERING parameters: ',"\n";
    print '~~~ Performing a clustering with ',$Linkage[$LinkCode],' Linkage method',"\n";
    
    if(($StopStep==0)&&($StopThr==0)){
	$Stop=0;
	print '~~~ Clustering with no stopping point',"\n";
	$StopLabel="NoStop";
    }elsif(($StopStep!=0)&&($StopThr!=0)){
	print '~~~ Ambiguous stopping criteria...',"\n";
	&abort();
    }elsif($StopStep>0){
	$Stop=1;
	print '~~~ Clustering and recovering classification at step: ',$StopStep,"\n";
	$StopLabel="StopStep-".$StopStep;
    }elsif($StopThr>0){
	$Stop=1;
	print '~~~ Clustering and recovering classification at similarity threshold: ',$StopThr,"\n";
	$StopLabel="StopThr-".$StopThr;
    }    
    print '  ',"\n";
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',"\n";
    print '  ',"\n";
}

######################
#         filesInSetUp
######################
# Open a file obtained from standard input variable, process the label
# and return the content in an array

sub filesInSetUp{

    my $pathIn=$_[0];
    print join(' ','~~~ Input path: ',$pathIn),"\n";
    open(INTMP,$pathIn)  or die "Can't find file $pathIn\n";
    @fields=split("/",$pathIn);
    $fileIn=$fields[$#fields]; # Take the last field of the path as the name of the file for outputs
    print join(' ','~~~ Input file: ',$fileIn),"\n";
    @fields=split("Similarities",$fileIn);
    if($fields[1]){
	$labelIn=$fields[1];
    }else{
	$labelIn='.out'
    }
    print join(' ','~~~ Label for outputs: ',$labelIn),"\n";
    print '  ',"\n";
    @INTMP = <INTMP>;
    close(INTMP);
    return @INTMP,$fileIn,$labelIn;
}

######################
#        filesOutSetUp
######################
# Define handlers for output files, open them and
# creates a header with information of parameters and date. It does
# return nothing, because the array of handlers is not used, I've
# used the names directly (there is no loop and not many multiple prints)


sub filesOutSetUp(){

    # -- Define here file handlers
    $Handle[0]=*HISTcompact; $Handle[1]=*HISTextend; $Handle[2]=*CLU; $Handle[3]=*PAR;

    # -- Invoke the time to build headers
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $year = 1900 + $yearOffset;
    $theTime = " (mm,dd,yy) $month, $dayOfMonth, $year, and time: $hour:$minute:$second,";

    # -- Define here the file names
    $fileOut[0]='HistCompact-NL_'.$Linkage[$LinkCode].'_'.$StopLabel.$labelIn1;
    $fileOut[1]='HistExtend-NL_'.$Linkage[$LinkCode].'_'.$StopLabel.$labelIn1;
    $fileOut[2]='Clusters-NL_'.$Linkage[$LinkCode].'_'.$StopLabel.$labelIn1;
    $fileOut[3]='Partition-NL_'.$Linkage[$LinkCode].'_'.$StopLabel.$labelIn1;

    for($u=0; $u<=3; $u++){
	if($Stop==0){
	    if($u==2){
		last;
	    }
	}
	$ASA=$Handle[$u];
	$file=$fileOut[$u];
	open($ASA, ">$file") || die "Couldn't open file $file"; 
	print $ASA '# >> Output from NodeLinkage.pl <<',"\n";
	print $ASA '# >> Input file 1: "',$fileIn1,"\n";
	print $ASA '# >> Input file 2: "',$fileIn2,"\n";
	print $ASA '# >> Clustering of nodes with the ',$Linkage[$LinkCode],' Linkage method',"\n";
	if($StopStep>0){
	    print  $ASA '# >> Recovering classification at step: ',$StopStep,"\n";
	}elsif($StopThr>0){
	    print  $ASA '# >> Recovering classification at similarity threshold: ',$StopThr,"\n";
	}
	print $ASA '# >> Running at date: ',$theTime,"\n";
	print $ASA '# >> Silwood Park, Imperial College London, A.P-G.', "\n";
	print $ASA '# ', "\n";
	if($u==0){
	    print $ASA '#1Step, 2maxGlobalTmp, 3Density, 4DensityInt, 5DensityExt, 6NumNodesA, 7NumEdgesA, 8NumNodesB, 9NumEdgesB, 10NumNodesAB, 11NumEdgesAB, 12NumIntNodesA, 13NumIntNodesB, 14NumExtNodesA, 15NumExtNodesB, 16NumIntNodesAB, 17NumExtNodesAB, 18NumIntEdgesA, 19NumIntEdgesB, 20NumExtEdgesA, 21NumExtEdgesB, 22NumIntEdgesAB, 23NumExtEdgesAB, 20nodeA, 21nodeB, 22NcumInt, 22NcumExt, 23Ncum',"\n";
	}elsif($u==1){
	    print $ASA '#CODE4VALUES_1Step, 2maxGlobalTmp, 3nodeA, 4nodeB, 5NumIntNodesA, 6NumIntNodesB, 7NumExtNodesA, 8NumExtNodesB, 9NumIntEdgesA, 10NumIntEdgesB, 11NumExtEdgesA, 12NumExtEdgesB',"\n";
	}elsif($u==3){
	    print $ASA '# >>1Node,  2ClusterId',"\n";
	}
    }
}


###########################
#   splitIntExtEdges
###########################
# Given a cluster with N nodes and m links, we identify how many
# links connect the nodes within the cluster, how many connect with nodes
# outside the cluster, hoy many nodes have connections with other nodes 
# in the cluster, and how many external nodes are involved in these interactions

sub splitIntExtEdges{
    my ($arrayRef1 , $arrayRef2) = @_; 
    my @subEdges = @{$arrayRef1};  
    my @subNodes = @{$arrayRef2};
    print join(" ",'>>INTO> **EDGES',@subEdges,'**NODES:', @subNodes),"\n"; # DEBUG
    foreach $EdgeTmp(@subEdges){
	@fields=split("XXXX",$EdgeTmp);
	$NodeA=$fields[0];
	$NodeB=$fields[1];
	$keyA=0;
	$keyB=0;
	foreach $NodeTmp(@subNodes){
	    if($NodeTmp eq $NodeA){
		$keyA=1;
	    }elsif($NodeTmp eq $NodeB){
		$keyB=1;
	    }
	    if(($keyA==1)&&($keyB==1)){
		$IntNode2key{$NodeA}=1;
		$IntNode2key{$NodeB}=1;
		$IntEdge2key{$EdgeTmp}=1;
		last;
	    }	    
	}
	if(($keyA==1)&&($keyB==0)){
	    #$IntNode2key{$NodeA}=1; # a is internal node
	    $ExtNode2key{$NodeB}=1; # b is external  node
	    $ExtEdge2key{$EdgeTmp}=1;
	}elsif(($keyA==0)&&($keyB==1)){
	    #$IntNode2key{$NodeB}=1; # b is internal node
	    $ExtNode2key{$NodeA}=1; # a is the external  node
	    $ExtEdge2key{$EdgeTmp}=1;
	}elsif(($keyA==0)&&($keyB==0)){
	    print "~~~ Error in splitIntExtInteractions function","\n";
	    &abort();
	}	
    }
    @IntNodes=keys(%IntNode2key); # These hashes are not returned, but they might be useful in the future
    @IntEdges=keys(%IntEdge2key);
    @ExtNodes=keys(%ExtNode2key); # These hashes are not returned, but they might be useful in the future
    @ExtEdges=keys(%ExtEdge2key);
    undef(%IntNode2key);
    undef(%ExtNode2key);
    undef(%IntEdge2key);
    undef(%ExtEdge2key);
    #print join(' ','>>',@IntNodes,@ExtNodes),"\n"; # DEBUG
    $NumNodes=$#subNodes+1;
    $NumEdges=$#subEdges+1;
    $NumIntNodes=$#IntNodes+1;
    $NumIntEdges=$#IntEdges+1;
    $NumExtNodes=$#ExtNodes+1;
    $NumExtEdges=$#ExtEdges+1;
    #print join(' ',$NumNodes,$NumEdges,$NumIntNodes,$NumIntEdges,$NumExtNodes,$NumExtEdges),"\n"; # DEBUG
    return $NumNodes,$NumEdges,$NumIntNodes,$NumIntEdges,$NumExtNodes,$NumExtEdges;
}


######################
#      clusterDensity
######################
# For a given cluster, it computes its density. The partition density
# will be the average of the density of the clusters belonging to the
# partition. 
######## UPDATE OCTOBER 2016
# Note that NintNodes counts the number of nodes within the cluster
# having an interaction with another node in the cluster. This variable
# was used before but it is not used any more under this definition, 
# because for the partition density it is needed the number of nodes
# in the cluster instead (independently of having or not a link). This 
# is why the "if" control structures are commented as well. In addition, 
# the variable Nedges is no longer used neither.
#

sub clusterDensity{
     my ($Nnodes,$Nedges,$NintNodes,$NintEdges,$NextNodes,$NextEdges) = @_;
     # print join(" ",'$Nnodes,$Nedges,$NintNodes,$NintEdges,$NextNodes,$NextEdges'),"\n"; # DEBUG
     # print join(" ",$Nnodes,$Nedges,$NintNodes,$NintEdges,$NextNodes,$NextEdges),"\n"; # DEBUG

     if($Nnodes==1){
	 $D1=0;
	 $D2=0;
     }elsif($NextNodes==0){
	 #if($NintNodes > 1){ # For August definition I made ($NintNodes > 2)
	     #$D1=($NintEdges-($NintNodes-1))/(($NintNodes-2)*($NintNodes-1)); # August 2016
	     $D1=2*$NintEdges/($Nnodes*($Nnodes-1)); # October 2016
	 # }else{
	 #     $D1=0;
	 # }
     	 $D2=0;
     	 #$D1=2*$NintEdges/($Nnodes*($Nnodes-1)); # Old definition
     }else{
	# if($NintNodes > 1){  # For August definition I made ($NintNodes > 2)
	     #$D1=($NintEdges-($NintNodes-1))/(($NintNodes-2)*($NintNodes-1)); # August 2016
	     $D1=2*$NintEdges/($Nnodes*($Nnodes-1)); # October 2016

	 # }else{
	 #     $D1=0;
	 # }
	 $D2=($NextEdges-$NextNodes)/($NextNodes*($Nnodes-1));
     	 #$D2=$NextEdges/(2*$Nnodes*$NextNodes); # Old definition
     	 #$D1=2*$NintEdges/($Nnodes*($Nnodes-1)); # Old definition
     } 

     #$D=$D1+$D2; # old definition
     #$D=$D*($NextEdges+2*$NintEdges)/2; # old definition
     $D1=$NintEdges*$D1;
     $D2=$NextEdges/2*$D2;
     $D=$D1+$D2;
     if($D <0){
	 print '~~~ Total density became negative',"\n";
	 print join(" ",'D', 'D1','D2','NintNodes','NextNodes','NintEdges','NextEdges'),"\n";
	 print join(" ",$D,$D1,$D2,$NintNodes,$NextNodes,$NintEdges,$NextEdges),"\n";
	 &abort();
    }
     return $D,$D1,$D2;
}

######################
#               abort
######################
# prints a warning message and aborts script execution 

sub abort{
    print '~~~ I abort the execution...',"\n";
    print '~~~ Exit!',"\n";
    print '  ',"\n";
    exit
}
