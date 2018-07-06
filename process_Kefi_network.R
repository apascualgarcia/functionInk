#####################################
# process_Kefi_network.R
####################################
#
# With this script we aim to process the original data
# and results published by Kefi et al. Plos Biology 2016.

library(reshape2)
# Load data -----------

dir="/home/apascual/Nextcloud/Research/Projects/FunctionalLinkage/SoniaKefi"
setwd(dir)
file=matrix(nrow=3,ncol=1)
Type=matrix(nrow=3,ncol=1)
file[1]="chilean_Trophic.txt"
Type[1]="Trophic"
file[2]="chilean_NonTrophic_pos.txt"
Type[2]="Pos"
file[3]="chilean_NonTrophic_neg.txt"
Type[3]="Neg"
all.long=list()
for(i in 1:3){
  mat=read.table(file=file[i],sep="\t",skip=1)
  rownames(mat)=mat$V2
  mat=subset(mat,select=-c(V1))
  colnames(mat)[1]="SpeciesA"
  colnames(mat)[2:dim(mat)[2]]=as.character(mat$Species)
  long=reshape(mat,direction="long",idvar="SpeciesA",varying=list(colnames(mat)[2:dim(mat)[2]]),
               timevar="SpeciesB",times=colnames(mat)[2:dim(mat)[2]],v.names="Weight") # I had problems with melt
  nn=dim(long)[1]
  type.vec=replicate(nn,Type[i])
  long=cbind(long,type.vec)
  long=long[abs(long$Weight)>0,]
  all.long[[i]]=long
  if(i == 1){
    all.long.short=all.long[[i]]
  }else{
    all.long.short=rbind(all.long.short,all.long[[i]])
  }
}

file.out="Network_chilean_merged.tsv"
write.table(all.long.short, file=file.out, sep="\t",row.names = FALSE,quote = FALSE)
