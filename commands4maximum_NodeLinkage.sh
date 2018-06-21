#!/bin/bash
# This script extracts the "number" maximum values 
# of the partition density computed in a "HistoryCompact" command 
# obtained from NodeLinkage.pl
###########################
# Silwood Park (ICL)
# July 7th, 2016
# Alberto Pascual-García
##########################
#
# in the awk command: 
# $3 is the whole partition density
# $4 the internal partition density
# $5 the external partition density

# --- generic labels for files
number=20 # Number of maximum values you may want to extract
method=Average
# --- specific label for functional groups
function=ATP7 #
# --- define your files here
#file=HistCompact-NL_"$method"_NoStop.out
#file=HistCompact-NL_"$method"_NoStop_FunctGroups_SingleFunctions_"$function"_Zeros-Intercept.txt
file=HistCompact-NL_Average_NoStop_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt

echo
echo Extracting the $number maximum values for function $function:
echo
echo Step / Threshold / Partition Density
echo

grep -v '#' $file | awk '{print $1,$2,$3}' | sort -k 3,3 | tail -n $number

echo

