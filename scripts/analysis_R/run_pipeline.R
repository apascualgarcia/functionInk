############################
#   run_pipeline.R
############################
#
# This script is a wrapper to run the whole functionink pipeline. Note that
# it simply perform calls to the system to run the Perl scripts, so the
# installation of Perl should be ready. The data generated during the
# process will be stored in files in a temporal directory. As a default, it will
# detect communities at the maximum of all partition densities.
#### INPUT:
# pathNet: An absolute path to the file containing the network, or relative
#         to the root path of the repository.
# fileNet: The name of the file containing the network
# pathRepo: An absolute path to the root of functionink repository
#### OPTIONS:
# weighted: logical, is the network weighted? defaults to FALSE
# directed: logical, is the network directed? defaults to FALSE
# types: logical, has the network different types of links? defaults to FALSE
# method: character, method of clustering, one of  "Average", "Single" or "Complete"
#         linkage. Defaults to average linkage.
# mode: character, criteria to extract the communities depending on the values of
#      the partition density. It will select the step where the maximum of the
#      "total", "internal", or "external", partition densities peak. Select "all"
#       for the maximum of the three partition densities (default) or "none" if you
#      are not interested in extracting the communities.
#
############################
# USAGE: Provide path and name of your network and path for the repository and set options
# OUTPUT: A directory with all the functionink analysis
############################
# Zürich, November 2022
# Theoretical Biology, ETH
# apascualgarcia.github.io
###################################################

run_pipeline = function(fileNet,pathNet,pathRepo, # mandatory
                        weighted=FALSE,directed=FALSE,
                        types=FALSE,method="Average",mode="all"){
  # ... set up the environment and source dependencies
  setwd(pathRepo)
  src.dir=paste("scripts","analysis_R",sep="/")
  setwd(src.dir)
  source("extractPartDensity.R")
  setwd(pathRepo)
  setwd(pathNet)
  dir.create("functionink_tmp")
  setwd("functionink_tmp")
  fileNetPath=paste0("../",fileNet)
  
  # ... process arguments
  if(weighted == FALSE){par_w = 0}else{par_w = 1}
  if(directed == FALSE){par_d = 0}else{par_d = 1}
  if(types == FALSE){par_t = 0}else{par_t = 1}

  # ... build commands basic run
  # ...... Node similarity
  script=paste(pathRepo,"NodeSimilarity.pl",sep="/")
  options=paste("-w",par_w,"-d",par_d,"-t",par_t,"-f",fileNetPath)
  comm_sim=paste(script,options)
  # ..... Node linkage
  fileSim=paste0("Nodes-Similarities_",fileNet)
  script=paste(pathRepo,"NodeLinkage.pl",sep="/")
  options=paste("-fn",fileNetPath,"-fs",fileSim,"-a",method)
  comm_link_base=paste(script,options)
  file.hist=paste0("HistCompact-NL_",method,"_NoStop_",fileNet) # expected history file output
  
  # --- Run the first analysis
  system(comm_sim) # compute nodes' similarities
  system(comm_link_base) # cluster nodes
  hist.comp=read.table(file=file.hist,sep="\t",header = TRUE) # read history file
  part_density=extractPartDensity(hist.comp) # extract partition densities
  
  # --- Run the extraction of the communities
  if(mode != "none"){ # if the user wants to retrieve the communities
    # .... determine the criteria to  be used
    labels=c("total_dens_step","int_dens_step","ext_dens_step")
    if(mode == "total"){
      idx=1
    }else if(mode == "internal"){
      idx=2
    }else if(mode == "external"){
      idx=3
    }else{ # mode "all", we identify the maximum among modes
      idx=which.max(c(part_density$total_dens,
                  part_density$int_dens,
                  part_density$ext_dens))
    }
    # ... find the step of the peak and create a new command
    value=part_density[labels[idx]]
    comm_link_spec=paste(comm_link_base,"-s step -v",value)
    # ... finally run
    system(comm_link_spec) # cluster nodes and extract partition
  }
  return(part_density)
}
