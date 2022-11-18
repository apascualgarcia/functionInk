extractPartDensity = function(hist.comp,plot = TRUE){
  # Function to retrieve the maximum of the difference partition density 
  # metrics.
  # INPUT:
  #  hist.comp = a data.frame with the contents of the file "HistCompact" obtained
  #     running the script NodeLinkage.pl with no optional arguments.
  #  plot = logical. If true it will create a plot with the partition density measures.
  # OUTPUT:
  #   A list with the maxima of each partition density (total, internal and external)
  ####################
  #
  require(ggplot2)
  
  # --- Retrieve maximum and report
  maxDens=hist.comp$Density[which.max(hist.comp$Density)]
  maxDensExt=hist.comp$DensityExt[which.max(hist.comp$DensityExt)]
  maxDensInt=hist.comp$DensityInt[which.max(hist.comp$DensityInt)]
  maxDensStep=hist.comp$Step[which.max(hist.comp$Density)]
  maxDensExtStep=hist.comp$Step[which.max(hist.comp$DensityExt)]
  maxDensIntStep=hist.comp$Step[which.max(hist.comp$DensityInt)]
  
  mes=paste("-- The maximum value of the total partition density is",round(maxDens,digits = 4),"found at step =",maxDensStep)
  print(mes)
  mes=paste("-- The maximum value of the internal partition density is",round(maxDensInt,digits=4),"found at step =",maxDensIntStep)
  print(mes)
  mes=paste("-- The maximum value of the external partition density is",round(maxDensExt,digits=4),"found at step =",maxDensExtStep)
  print(mes)

  
  if(plot == TRUE){
    system(command="mkdir -p figures")
    this.dir=getwd()
    setwd("figures")
    fileOut=paste("Plot_PartitionDensityVsStep.pdf",sep="")
    pdf(file=fileOut,width=10,height=8)
    print(ggplot(data = hist.comp) + 
            geom_point(mapping = aes(x=Step,y=Density,color="Total"))+
            geom_point(mapping = aes(x=Step,y=DensityInt,color="Internal"),shape=2)+
            geom_point(mapping = aes(x=Step,y=DensityExt,color="External"),shape=6) + 
            geom_line(mapping = aes(x=Step,y=Density,color="Total"))+
            geom_line(mapping = aes(x=Step,y=DensityInt,color="Internal"))+
            geom_line(mapping = aes(x=Step,y=DensityExt,color="External"))+
            xlab("Clustering Step")+
            theme(legend.title =  element_blank(),
                  axis.text = element_text(size=20),axis.title = element_text(size=24),
                  legend.text =element_text(size=18)))
    dev.off()
    setwd(this.dir)
  }
  
  return(list("total_dens" = maxDens,"int_dens"=maxDensInt,"ext_dens"=maxDensExt,
              "total_dens_step" = maxDensStep,"int_dens_step"=maxDensIntStep,"ext_dens_step"=maxDensExtStep))
}