############################
#   nodeLinkage_analysis.R
############################
#
# This script is an example of how to use the function extractPartDensity.R to
# extract the maximum values of the partition densities measures. It is assumed
# that the user is following the Vignette described in the documentation, and
# that this script is executed with Rstudio.
############################
# USAGE: Provide path and name of the history compact file as indicated below.
# OUTPUT: The step at which the partition densities are maximum, and a summary plot of the partition densities.
############################

library(ggplot2)
library(here)

#### START EDITING
# --- Path to dir of history compact file
dir="." #  path to dir of history compact file relative to the root of the repo (e.g. fix to "." if it is located in the root of the directory)

# --- Name of the history file
file.hist="HistCompact-NL_Average_NoStop_Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.txt" #"history file name"

# ---- STOP EDITING HERE 
src.dir=here("scripts","analysis_R")
root.dir=here()
setwd(src.dir)
source("extractPartDensity.R")
if(dir == "."){
  setwd(root.dir)
}else{
  setwd(here(dir))
}

hist.comp=read.table(file=file.hist,sep="\t",header = TRUE) # for current NodeLink.pl version (Dec 2018)

part_density=extractPartDensity(hist.comp)

part_density$total_dens # maximum of the total partition density
part_density$int_dens # maximum of the internal partition density
part_density$ext_dens # maximum of the external partition density
part_density$total_dens_step # step of the clustering in which the maximum of the total partition density was found
part_density$int_dens_step # step of the clustering in which the maximum of the internal partition density was found
part_density$ext_dens_step # step of the clustering in which the maximum of the external partition density was found

