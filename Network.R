#######################download all libraries###################################
if(!require(circlize)){install.packages("circlize")}
library(circlize)
if(!require(STRINGdb)){BiocManager::install("STRINGdb",update=T,ask=F)}
library(STRINGdb)
if(!require(gplots)){install.packages("gplots")}
library(gplots)
if(!require(plyr)){install.packages("plyr")}
library(plyr)
if(!require(tiff)){install.packages("tiff")}
library(tiff)
################################################################################

#######################protein network##########################################
strspecies=read.delim("https://stringdb-static.org/download/species.v11.5.txt")
#check for species, download respective species data from STRING database
if (phylo %in% strspecies$official_name_NCBI) {
dir.create("Network",showWarning=F)
msgBox("Now proceeding to Network Analysis.")
spcid=strspecies$X.taxon_id[match(phylo,strspecies$official_name_NCBI)]
string_db <- STRINGdb$new( version="11.5", species=spcid ,score_threshold=0, input_directory="" )
#statsRes$Protein=pepalllabcleaned
#cleanedids=data.frame(Protein=pepalllabcleaned,Gene.names=statsRes$Gene.names)
#globmapping=string_db$map(cleanedids, "Protein", removeUnmappedRows = TRUE)

#get STRING ID for all proteins in the dataset
globmapping=string_db$map(statsRes[c("Protein","Gene.names")], "Protein", removeUnmappedRows = TRUE)
globmapping=ddply(globmapping,.(Protein),summarize,Gene.names=Gene.names[1],STRING_id=STRING_id[1])

#option to include hub analysis and hub enrichment analysis
hubyes=dlgMessage("Find hub nodes?",type="yesno")$res
if (hubyes=="yes") {
hubcof=dlgList(c("1 standard deviation above mean","top 10% highest interaction"),title="Hub selection cutoff")$res
}
hubenrichyes=dlgMessage("Check hub interaction enrichment for your specific protein set against background sets?",type="yesno")$res

#function for cleaning up accession IDs
for (i in 1:length(nwgenes)) {
cleanedsigacc=unlist(lapply(nwgenes[[i]]$Protein ,FUN=function(x) {
cleaned=c()
if (!is.na(x)) {
temp=strsplit(x,";")
for (k in 1:length(temp)) {
	if (sum(nchar(temp[[k]])<7)>0) {
cleaned=c(cleaned,temp[[k]][which(nchar(temp[[k]])<7)[1]])} 
	else if (grepl("-",temp[[k]][1])){
cleaned=c(cleaned,substr(temp[[k]][1],1,(gregexpr("-",temp[[k]][1])[[1]][1]-1)))}
	else {cleaned=c(cleaned,temp[[k]][1])}
#criteria for cleaning: choose the first one that's uniprot 6-symbol acc., if none, return first accession, if first has sub category, meaning dash, remove sub
} } else {cleaned=c(cleaned,NA)}
return(cleaned)}))
nwgenes[[i]]$Protein=cleanedsigacc

#map query accessions to STRING database to get interaction map
tempmapping=string_db$map(nwgenes[[i]][,c(1,3,7:ncol(nwgenes[[i]]))], "Protein", removeUnmappedRows = FALSE)
#record accessions that cannot be found in STRING database
if (sum(is.na(tempmapping$STRING_id))>0) {
write.csv(tempmapping[which(is.na(tempmapping$STRING_id)),],
paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_unmappedEntries.csv"))
}
#as long as some accessions are found, create tiff figure of the network for these accessions
if (sum(is.na(tempmapping$STRING_id))<length(tempmapping$STRING_id)) {
tempmapping2=tempmapping[which(!is.na(tempmapping$STRING_id)),]
tiff(paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_STRINGnetwork.tiff"), 
units="in", width=15, height=15, res=300)
tryCatch(string_db$plot_network(tempmapping2$STRING_id,required_score=400),error=function(e) {textplot("Connection Error!")})
dev.off()
#set confidence score cutoff to 400, if any interaction passes the threshold, it is recorded in a csv file 
tempinteraction=string_db$get_interactions(tempmapping2$STRING_id)
tempconfint=tempinteraction[tempinteraction$combined_score>400,]
if (nrow(tempconfint)>0) {
tempconfint$from_gene=tempmapping2$Gene.names[match(tempconfint$from,tempmapping2$STRING_id)]
tempconfint$to_gene=tempmapping2$Gene.names[match(tempconfint$to,tempmapping2$STRING_id)]
write.csv(tempconfint,paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_STRINGinteractions.csv"),row.names=F)
#tidy up interaction list, and, for two proteins with more than 1 type of interactions, create a variable holds the combined interaction, and highest confidence interaction is recorded
tempconfint2=tempconfint[which(!(tempconfint$from_gene==""|tempconfint$to_gene=="")),]
tempconfint2=ddply(tempconfint2,.(from, to), summarize, highest_score=max(combined_score),
combined_score=sum(combined_score), from_gene=from_gene[1],to_gene=to_gene[1])

#implementing background colors for nodes based on fold change, and creating tiff figures for these networks
temppvalidx=grep("Pvalue|FDR|PFP",colnames(tempmapping2))
tempfcidx=grep("FoldChange",colnames(tempmapping2))
if (length(temppvalidx)>0&length(tempfcidx)>0) {
nwsubcolor=string_db$add_diff_exp_color(tempmapping2[which(tempmapping2[,temppvalidx[1]]<0.05),],logFcColStr=colnames(tempmapping2)[tempfcidx[1]])
payload_id = string_db$post_payload(nwsubcolor$STRING_id,colors=nwsubcolor$color)
tiff(paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_STRINGnetwork_FCcolored.tiff"), 
units="in", width=15, height=15, res=300)
tryCatch(string_db$plot_network(tempmapping2$STRING_id,required_score=400,payload_id=payload_id),error=function(e) {textplot("Connection Error!")})
dev.off()
} else if (length(tempfcidx)>0) {
nwsubcolor=string_db$add_diff_exp_color(tempmapping2,logFcColStr=colnames(tempmapping2)[tempfcidx[1]])
payload_id = string_db$post_payload(nwsubcolor$STRING_id,colors=nwsubcolor$color)
tiff(paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_STRINGnetwork_FCcolored.tiff"), 
units="in", width=15, height=15, res=300)
tryCatch(string_db$plot_network(tempmapping2$STRING_id,required_score=400,payload_id=payload_id),error=function(e) {textplot("Connection Error!")})
dev.off()
} 

#check for interaction clusters in the network, if any cluster includes 10 or above proteins, a tiff figure is created for the clusters (4 clusters per figure)
clustersList = string_db$get_clusters(tempmapping2$STRING_id)
clustersize=unlist(lapply(clustersList,length))
if (sum(clustersize>10)>0) {
clustersList=clustersList[which(clustersize>10)]
for (j in 1:ceiling(sum(clustersize>10)/4)) {
tiff(paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_STRINGclusters",j,".tiff"), 
units="in", width=15, height=15, res=300)
par(mfrow=c(2,2))
for (k in 1:4) {
if ((j-1)*4+k<=length(clustersList)) {
tryCatch(string_db$plot_network(clustersList[[(j-1)*4+k]],required_score=400),error=function(e) {textplot("Connection Error!")})
}}
dev.off()
}
}

#hub analysis starts if selected yes
if (hubyes=="yes") {
#get interaction count for all accessions in the network and order them
tempintcount=table(c(tempconfint2$from_gene,tempconfint2$to_gene))
tempintcount=tempintcount[order(tempintcount,decreasing=T)]
#select hub based on 1 SD above mean when this option is chosen
if (match(hubcof,c("1 standard deviation above mean","top 10% highest interaction"))==1) {
temphubs=tempintcount[which(tempintcount>(mean(tempintcount)+sd(tempintcount))&(tempintcount>1))]
} 
#get top 10% protein with highest interaction counts, excluding proteins with only 1 interaction, when this option is chosen
else if (match(hubcof,c("1 standard deviation above mean","top 10% highest interaction"))==2) {
if (sum(tempintcount>1)>9) {
temphubs=tempintcount[1:floor(length(tempintcount)/10)] } else {temphubs=c()}
} else {temphubs=c()}

#if hub candidates are present from the chosen cutoff, a csv file is generated to include the following for each hub:
#interaction count for standard confidence (>400) interactions, proteins it interact with;
#interaction count for high confidence (>700) interactions, proteins it interact with
if (length(temphubs)>0) {
hubintprots=unlist(lapply(names(temphubs),FUN=function(x) paste(c(tempconfint2$to_gene[tempconfint2$from_gene %in% x],
tempconfint2$from_gene[tempconfint2$to_gene %in% x]),collapse=',')))
temphighint=tempconfint2[tempconfint2$highest_score>700,]
temphighintcount=table(c(temphighint$from_gene,temphighint$to_gene))
hhint=temphighintcount[match(names(temphubs),names(temphighintcount))]
hhint[is.na(hhint)]=0
hhintprots=unlist(lapply(names(temphubs),FUN=function(x) paste(c(temphighint$to_gene[temphighint$from_gene %in% x],
temphighint$from_gene[temphighint$to_gene %in% x]),collapse=',')))
temphubinfo=data.frame(Hub=names(temphubs),Interaction_count=as.vector(temphubs),Interactions=hubintprots,
High_confidence_interaction_count=as.vector(hhint),High_confidence_interaction=hhintprots)
#high confidence interaction cutoff=700,interaction cutoff=400
write.csv(temphubinfo,paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_Hubs_Info.csv"),row.names=F)

#starts hub enrichment analysis if yes is selected  
if (hubenrichyes=="yes") {
set.seed(123)
#extract mappings for all hubs and get their STRING ID 
tempglob=globmapping[-na.omit(which(globmapping$Gene.names %in% names(temphubs))),]
temphubid=globmapping$STRING_id[match(names(temphubs),globmapping$Gene.names)]

#generate 1000 background networks for each hub
bgintlist=c()
bgpval=c()
bgfdr=c()
for (j in 1: length(temphubs)) {
tempintlist=c(as.vector(temphubs)[j])
for (k in 1:1000) {
#for each of the 1000 background network, randomly select proteins from the dataset, until there are the same number of proteins as the query network input when adding the target hub
bgid=c(sample(tempglob$STRING_id,length(unique(tempmapping2$STRING_id))-1),temphubid[j])
#get STRING interaction for the background network, organize result data, and count standard confidence interactions for target hub 
tempintk=string_db$get_interactions(bgid)
tempconfintk=tempintk[tempintk$combined_score>400,]
tempconfintk$from_gene=globmapping$Gene.names[match(tempconfintk$from,globmapping$STRING_id)]
tempconfintk$to_gene=globmapping$Gene.names[match(tempconfintk$to,globmapping$STRING_id)]
tempconfintk=tempconfintk[which(!(tempconfintk$from_gene==""|tempconfintk$to_gene=="")),]
tempconfintk2=ddply(tempconfintk,.(from, to), summarize, highest_score=max(combined_score))

tempbgint=table(c(tempconfintk2$from,tempconfintk2$to))
#score 1000 background network interaction count for target hub 
tempintlist=c(tempintlist, as.vector(tempbgint[match(temphubid[j],names(tempbgint))]))
}
#if for some reason all interaction counts are 0 or NA, enrichment for this hub is not calculated
#if not enrichment pvalue and adjusted pvalue (by BH method) is calculated 
tempintlist[is.na(tempintlist)]=0
if (sum(tempintlist==0)>=1000) {
tempbgpval=NA
tempbgfdr=NA
} else {
tempbgpval=sum(tempintlist[2:1001]>tempintlist[1],na.rm=T)/(1000-sum(is.na(tempintlist)))
tempbgfdr=p.adjust(tempbgpval, method = "BH", n = length(temphubs))
}
#record 1001 interaction counts, pvalue, and adjusted pvalue for all hubs, and output a csv file that includes this information
bgintlist=rbind(bgintlist,tempintlist)
bgpval=rbind(bgpval,tempbgpval)
bgfdr=rbind(bgfdr,tempbgfdr)
}
bghubenrich=data.frame(bgpval,bgfdr,bgintlist)
colnames(bghubenrich)=c("P-value","FDR","query_network_interaction_count",paste0("background_network",1:1000,"_interaction_count"))
rownames(bghubenrich)=names(temphubs)
write.csv(bghubenrich,paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_Hubs_Interaction_Enrichment.csv"))

#organize data and plot hub enrichment with ggplot boxplot, and create a tiff and an EPS figure 
bghubenrich=bghubenrich[,-which(colnames(bghubenrich)=="P-value")]
bghubenrich$hubs=rownames(bghubenrich)

bghubenrichip=melt(bghubenrich,id.vars=c("hubs","FDR"))
colnames(bghubenrichip)[colnames(bghubenrichip)=="FDR"]="FDR"
labdat = ddply(bghubenrichip, .(hubs), summarize, ypos= max(value)*1.1, FDR=FDR[1])
	 
hubgg=ggplot(data=bghubenrichip,aes(x=hubs,y=value)) +
  geom_boxplot(notch=T) +
  geom_point(data=subset(bghubenrichip,variable=="query_network_interaction_count"),shape=8,color="red",size=1.5) +
  geom_text(data=labdat,aes(label = FDR, y = ypos), position = position_dodge(width = .75)) +
  labs(title="Hub intereaction enrichment comparing to 1000 background networks",x="Hubs", y = "Interaction Count") +
  theme(axis.text.x = element_text(angle = 90,hjust=0.98,vjust=0.2))

#automatically adjusting figure size based on number of hubs  
wd=10+(length(temphubs)-10)*0.18
ht=7
fname=paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_Hubs_Interaction_Enrichment.tiff")
tiff(filename =fname, units="in", width=wd, height=ht, res=300)
print(hubgg)
dev.off()
setEPS()
postscript(paste0("Network/",substring(names(nwgenes)[i],1,gregexpr("\\_",names(nwgenes)[i])[[1]][2]-1),"_Hubs_Interaction_Enrichment.eps"), width=wd, height=ht)
print(hubgg)
dev.off()

}
}
}
}
}
}

} else {msgBox("Organism not supported for network analysis.")}
################################################################################

#######################kinase network image#####################################
#available for phosphoproteomics data when there are significant networks based on KinSwingR analysis 
if (dir.exists("Kinase Analysis")) {
if (sum(grepl("significant_kinaseNetwork",list.files(path="Kinase Analysis",pattern="*.csv")))>0) {

#search for all significant kinase network files and read content
allf=list.files(path="Kinase Analysis",pattern="*.csv")

for (i in 1:length(allf[grep("significant_kinaseNetwork",allf)])) {
tgtname=allf[grep("significant_kinaseNetwork",allf)][i]
tempkinnw=read.csv(paste0("Kinase Analysis\\",tgtname),stringsAsFactors=F)
tempkinnw=tempkinnw[which(nchar(unlist(lapply(tempkinnw[,2],FUN=function(x) strsplit(x,"\\|")[[1]][4])))!=2),]

#if there are more than 250 interactions, due to image size, the top 250 is kept 
if (nrow(tempkinnw)>250) {
tempkinnw=tempkinnw[order(tempkinnw$Pvalue,decreasing=F),][1:250,]
}
#circosplots are made for any file with more than 10 significant interactions 
if (nrow(tempkinnw)>10&nrow(tempkinnw)<251) {
#input data and circosplot parameters are set 
temppeplab=data.frame(matrix(unlist(strsplit(tempkinnw[,2],"\\||\\::")),ncol=5,byrow=T),stringsAsFactors=F)
colnames(temppeplab)=c("accession","gene","position","seq","seq2")
temppeplab$phoseq=substring(temppeplab$seq,ceiling(nchar(temppeplab$seq)/2),ceiling(nchar(temppeplab$seq)/2))
temppeplab$site=paste0(temppeplab$phoseq,temppeplab$position)
knwgrp=c(rep("Kinase",nrow(temppeplab)),temppeplab$gene)
names(knwgrp)=c(tempkinnw$Kinase,temppeplab$site)
tempknwinfo=data.frame(kinase=tempkinnw$Kinase,site=temppeplab$site)
gap=10*(3/length(unique(knwgrp)))+0.3
choralpha=0.5

#figure size adjusted automatically 
wd=10+(length(knwgrp)-20)*0.1
if (wd>75) {wd=75}
ht=wd

#tiff figure is created for circosplot 
fname=paste0("Network/",substring(tgtname,1,gregexpr("\\.",tgtname)[[1]][1]-1),"_ChordDiagram.tiff")
tiff(filename =fname, units="in", width=wd, height=ht, compression="lzw", res=300)
#making the actual circosplot 
chordDiagram(tempknwinfo,group = knwgrp,transparency = choralpha,big.gap = gap,
    annotationTrack = c("grid", "axis"),
    preAllocateTracks = list(
        track.height = mm_h(4),
        track.margin = c(mm_h(4), 0)
))
circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, col="cyan", niceFacing = TRUE)
}, bg.border = NA)
for (j in 1:length(unique(knwgrp))) {
if (unique(knwgrp)[j]=="Kinase"){
highlight.sector(names(knwgrp[which(knwgrp==unique(knwgrp)[j])]), track.index = 1, cex = 0.8,col = "darkgoldenrod4",
    text = unique(knwgrp)[j],  text.col = "white", niceFacing = TRUE)
} else {
highlight.sector(names(knwgrp[which(knwgrp==unique(knwgrp)[j])]), track.index = 1, cex = 0.8,col = "aquamarine4",
    text = unique(knwgrp)[j],  text.col = "white", niceFacing = TRUE)
}
}
circos.clear()
dev.off()



#a second style circosplot with larger labels is created and saved in tiff figure file 
fname=paste0("Network/",substring(tgtname,1,gregexpr("\\.",tgtname)[[1]][1]-1),"_ChordDiagram_style2.tiff")
tiff(filename =fname, units="in", width=wd, height=ht, compression="lzw",res=300)
chordDiagram(tempknwinfo,group = knwgrp,transparency = choralpha,big.gap = gap,
    annotationTrack = c("grid", "axis"),
    preAllocateTracks = list(
        list(track.height = max(strwidth(c(tempknwinfo$kinase,tempknwinfo$site)))*5),
        list(track.height =max(strwidth(names(knwgrp)))*3, track.margin = c(0, 0))
))
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, cex=4.5,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
circos.track(track.index = 3, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, col="cyan", niceFacing = TRUE)
}, bg.border = NA)
for (j in 1:length(unique(knwgrp))) {
if (unique(knwgrp)[j]=="Kinase"){
highlight.sector(names(knwgrp[which(knwgrp==unique(knwgrp)[j])]), track.index = 2, cex = 4,col = "darkgoldenrod4",
    text = unique(knwgrp)[j],  text.col = "white", niceFacing = TRUE
	)
} else {
highlight.sector(names(knwgrp[which(knwgrp==unique(knwgrp)[j])]), track.index = 2, cex = 2.6,col = "aquamarine4",
    text = unique(knwgrp)[j],  text.col = "white", niceFacing = TRUE,facing = "clockwise")
}
}
circos.clear()
dev.off()

}	
}
}}
################################################################################
