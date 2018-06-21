#!/bin/bash
# This script takes the Extended files generated with NodeLink.pl and extracts
# the last two clusters joined (note that the number used in the first line
# with "VALUES" should be modify, this number is the number of elements minus one) and
# it prints it in a file including the name of the file where it was found. Then
# it does some cleaning in the file to leave just the parameters. This script
# should be improved defining variables.
##############
# March 30th,2017; Silwood Park
# Imperial College London
# Alberto Pascual-García.
#######################

find . -iname '*Extend*' | xargs -i grep -H -A 2 "VALUES 92" {} > ClustersLastStep_AvLink.dat
sed 's/\.\/HistExtend-NL_Average_NoStop_Rmrz89_Exhaust_//g' ClustersLastStep_AvLink.dat | sed 's/.long.txt-/: /' | sort | grep -v VALUES  > ClustersLastStep_AvLink.tmp
mv -f ClustersLastStep_AvLink.tmp ClustersLastStep_AvLink.dat
