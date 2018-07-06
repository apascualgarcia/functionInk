

library(ggplot2)

dir="/home/apascual/Nextcloud/Research/Projects/FunctionalLinkage/SoniaKefi"
setwd(dir)
file.hist="HistCompact-NL_Average_NoStop_Network_chilean_merged.tsv"
hist.comp=read.table(file=file.hist,sep="\t",header = TRUE)

ggplot(data = hist.comp) + 
  geom_point(mapping = aes(x=Step,y=Density,color="Total"))+
  geom_point(mapping = aes(x=Step,y=DensityInt,color="Internal"),shape=2)+
  geom_point(mapping = aes(x=Step,y=DensityExt,color="External"),shape=6) + 
  geom_line(mapping = aes(x=Step,y=Density,color="Total"))+
  geom_line(mapping = aes(x=Step,y=DensityInt,color="Internal"))+
  geom_line(mapping = aes(x=Step,y=DensityExt,color="External")) 

maxDens=hist.comp$Step[which.max(hist.comp$Density)]
maxDensExt=hist.comp$Step[which.max(hist.comp$DensityExt)]
maxDensInt=hist.comp$Step[which.max(hist.comp$DensityInt)]


print("Maximum total partition density found at step =")
print(maxDens,digits=3)
print("Maximum external partition density found at step =")
print(maxDensExt,digits=3)
print("Maximum internal partition density found at step =")
print(maxDensInt,digits=3)
