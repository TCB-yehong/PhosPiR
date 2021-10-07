#######################download all libraries###################################
if(!require(reshape2)){install.packages("reshape2")}
library(reshape2)
################################################################################

#two types of special input:
#suspected falsely-classified-contamination included in input
#expanded phosphorylation count entries in input 
inputselec=dlgList(c("Expanded phosphorylation count","Suspected falsely-classified-contamination included"),multiple=T,title="Select special input type")$res
dir.create("Raw File",showWarning=F)
if ("Expanded phosphorylation count" %in% inputselec) {
raw=read.delim(fname)
#intensity file:
#extract data info column number
didx1=match(c("Protein","Protein.names","Gene.names","Amino.acid","Position",
"Sequence.window","Reverse","Potential.contaminant"),colnames(raw))
#extract all intensity column number
intidx=grep("Intensity",colnames(raw))
#find out which columns are the expanded phospho count intensities for each sample
intidx2=which(diff(intidx)>1)
#extract intensity column number of expanded columns for all samples
didx2=intidx[(intidx2[2]+1):length(intidx)]
#combine with label columns 
intsep=raw[,c(didx1,didx2)]
#check how many phosphorylation counts the intensity is separated into
intsepgrp=length(table(unlist(lapply(colnames(intsep)[(length(didx1)+1):ncol(intsep)],FUN=function(x) strsplit(x,"__")[[1]][2]))))
n=1
#for each start of the sample columns, extract target sample columns and IDs, change to long format, then merge with other sample's data
for (i in seq(length(didx1)+1,ncol(intsep),intsepgrp)) {
tempsepdf=intsep[,c(1:length(didx1),i:(i+intsepgrp-1))]
tempsepdfm=melt(tempsepdf,id.vars=colnames(tempsepdf)[1:length(didx1)],variable.name="expandID",value.name="Intensity")
#tempsepmo=tempsepdfm[with(tempsepdfm,order(Protein,Amino.acid,Position,expandID)),]
tempsepmo=tempsepdfm
#add phospho-count information to protein description
tempsepmo$Protein.names=paste0(tempsepmo$Protein.names,"; phospho-count=",unlist(lapply(as.character(tempsepmo$expandID), FUN=function(x) strsplit(x,"___")[[1]][2])))
#change column "Intensty" to sample name
colnames(tempsepmo)[colnames(tempsepmo)=="Intensity"]=strsplit(names(table(tempsepmo$expandID))[1],"___")[[1]][1]
#remove column of expanded columns' names
tempsepmo=tempsepmo[,-which(colnames(tempsepmo)=="expandID")]
#merge with previous samples' data, if it's the first sample, simply record the data
if (n!=1) {
finalsepdf=merge(finalsepdf,tempsepmo,by=colnames(tempsepdf)[1:length(didx1)],all=T)} else {
finalsepdf=tempsepmo
}
n=n+1
}
#sort by IDs
finalsepdf=finalsepdf[with(finalsepdf,order(Protein,Amino.acid,Position,Protein.names)),]
rawdatafile=finalsepdf

#remove reverse and contaminant entries
rawdatafile$Reverse[is.na(rawdatafile$Reverse)]=""
rawdatafile$Potential.contaminant[is.na(rawdatafile$Potential.contaminant)]=""
rawdatafile=rawdatafile[which(rawdatafile$Reverse==""&rawdatafile$Potential.contaminant==""),-c(7,8)]
#get original dataset name
ofname=strsplit(fname,"/|\\.")[[1]]
ofnames=ofname[length(ofname)-1]
#save generated input file 
write.csv(rawdatafile,paste0("Raw File/",ofnames,"_expanded_phospho_count.csv"),row.names=F)
}

if ("Suspected falsely-classified-contamination included" %in% inputselec) {
raw=read.delim(fname)

#intensity file:
#extract data info column number
didx1=match(c("Protein","Protein.names","Gene.names","Amino.acid","Position",
"Sequence.window","Reverse","Potential.contaminant"),colnames(raw))
#extract all intensity column number
intidx=grep("Intensity",colnames(raw))
#find out which columns are the combined intensities for each sample
intidx2=which(diff(intidx)>1)
#extract intensity column number of combined intensities for each sample
didx2=intidx[(intidx2[1]+1):intidx2[2]]
rawdatafile=raw[,c(didx1,didx2)]
rawdatafile$Reverse[is.na(rawdatafile$Reverse)]=""
rawdatafile$Potential.contaminant[is.na(rawdatafile$Potential.contaminant)]=""
#extract potential contaminant entries 
pcontam=rawdatafile[which(rawdatafile$Potential.contaminant=="+"),]
#if no contaminants, print warning 
if (nrow(pcontam)==0) {
dlgMessage("There are no marked contaminants in the data!")
#save generated input file 
write.csv(rawdatafile[,-c(7,8)],paste0("Raw File/",ofnames,"_no_potential_contaminant.csv"),row.names=F)
} else {
#mark possible contaminant in protein name 
pcontam$Protein.names=paste0(pcontam$Protein.names,"; Possible contaminant")
#keep only non-reverse entries 
pcontam=pcontam[which(pcontam$Reverse==""),]
#remove "CON__" from protein name
pcomprot=gsub("CON__","",pcontam$Protein)
#get non-reverse non-contaminant protein names from dataset 
rawprot=rawdatafile[which(rawdatafile$Reverse==""&rawdatafile$Potential.contaminant==""),which(colnames(rawdatafile)=="Protein")]
#keep any potential contaminant entries with overlapping protein ID with non-contaminant
pcontamf=pcontam[which(pcomprot %in% rawprot),-c(7,8)]

#remove reverse or contaminant entries from main dataset 
rawdatafile=rawdatafile[which(rawdatafile$Reverse==""&rawdatafile$Potential.contaminant==""),-c(7,8)]
#add kept potential contaminant entries
rawdatafile=as.data.frame(rbind(rawdatafile,pcontamf))
#get original dataset name
ofname=strsplit(fname,"/|\\.")[[1]]
ofnames=ofname[length(ofname)-1]
#save generated input file 
write.csv(rawdatafile,paste0("Raw File/",ofnames,"_potential_false_contaminant_added.csv"),row.names=F)
#report number of entries added to the dataset 
dlgMessage(paste0(nrow(pcontamf)," potential false contaminant-labeled entries are added to dataset."))}
}
