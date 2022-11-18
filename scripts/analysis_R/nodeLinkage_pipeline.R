############################
#   nodeLinkage_pipeline.R
############################
#
# This script is an example of how to use the function run_pipeline.R to
# compute the whole functionink pipeline. It is assumed
# that the user is following the Vignette described in the documentation and
# hence that the repository was cloned, and that this script is executed with Rstudio.
# The example is set to be executed in its location in the repository and hence
# the paths are automatically set, but you should accomodate these paths to your
# environment.
############################
# USAGE: Provide path and name of the network and of the repository, and set the options desired.
# OUTPUT: A directory located in the same path as the network, containing all output files
############################

#library(here) # you can skip if you use absolute paths

#### START EDITING
# --- Path to repository
pathRepo="." #  absolute path to the root of the repo, or "." if your working directory already is the repo.

# --- Name and path of the file containing the network
pathNet="data" # absolute path, or one relative to the root of the repo
fileNet="Rmrz89_Exhaust_Nest0.6_Conn0.08_CompConn0.15.long.format2.txt"

# ---- STOP EDITING HERE 

# --- Set up the paths
if(pathRepo != "."){
  setwd(pathRepo)
}
src.dir=paste("scripts","analysis_R",sep="/")
setwd(src.dir)
source("extractPartDensity.R")
source("run_pipeline.R")


# --- Run the pipeline
run_pipeline(fileNet,pathNet,pathRepo, # mandatory
weighted=TRUE,directed=FALSE, # options set as in the vignette
types=TRUE,method="Average",mode="all")
