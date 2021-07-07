#######################download all libraries###################################
if(!require(reshape2)){install.packages("reshape2")}
library(reshape2)
if(!require(pheatmap)){install.packages("pheatmap")}
library(pheatmap)
if(!require(rcdk)){install.packages("rcdk")}
library(rcdk)
if(!require(fingerprint)){install.packages("fingerprint")}
library(fingerprint)
if(!require(vegan)){install.packages("vegan")}
library(vegan)
if(!require(rgl)){install.packages("rgl")}
library(rgl)
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
if(!require(plot3D)){install.packages("plot3D")}
library(plot3D)
if(!require(magick)){install.packages("magick")}
library(magick)
if(!require(plot3Drgl)){install.packages("plot3Drgl")}
library(plot3Drgl)
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
library(RColorBrewer)
if(!require(FactoMineR)){install.packages("FactoMineR")}
library(FactoMineR)
if(!require(factoextra)){install.packages("factoextra")}
library(factoextra)
################################################################################

msgBox("Now proceeding to Overview Figure step.")
dir.create("Overview Figure",showWarnings=F)
colrgrp=dlgList(names(group),title="Choose color grouping for figures")$res
gcolor=group[[match(colrgrp,names(group))]]

grapov=function(x,type,group) {
#x=dataset, type=graph type
#type=c("heatmap","PCA","PCA3d","histogram","boxplot")
x[is.na(x)]=0
hiss=list()
a=x
a$name=rownames(x)
forgraph=melt(a,id.vars="name")
l2=forgraph$value
forgraph$Log2values=log2(l2)

if (type=="heatmap") {
if (nrow(x)>10000) {
kmeansk=floor(nrow(x)/10)} else {kmeansk=NA}
fname="Overview Figure/Overview_heatmap.tiff"
hmp=pheatmap(log2(x),show_rownames=F,cluster_rows=F,kmeans_k=kmeansk,main="Data Overview - Log 2 Intensity",filename=fname,)
}

else if (type=="PCA") {
kmc=length(unique(group))+1
tdata=t(x)
km <-kmeans(tdata,kmc,10000)
pc<-prcomp(tdata)
fname="Overview Figure/PCA-kmeans_plot.tiff"
tiff(filename =fname, units="in", width=12, height=10, res=300)
plot(pc$x[,1], pc$x[,2],col=group,pch=16,main="K-means Clustering of PCA", xlab="Principal Component 1", ylab="Principal Component 2")
text(pc$x[,1], pc$x[,2],labels=rownames(pc$x),pos=3,offset=0.4,cex=0.7)
pc<-cbind(pc$x[,1], pc$x[,2])
ordispider(pc, factor(km$cluster), label = TRUE)
ordihull(pc, factor(km$cluster), lty = "dotted")
dev.off()

accpca=PCA(tdata,graph=F)
tiff("Overview Figure/PCA_plot.tiff", units="in", width=12, height=10, res=300)
print(fviz_pca_ind(accpca,geom.ind="point",label="all",col.ind=factor(group),palette="Set1",addEllipses = TRUE,ellipse.type="confidence",ellipse.alpha=0.3,legend.title = "Groups"))
dev.off()
}

else if (type=="PCA3d") {
tdata=t(x)
cmd3d<-cmdscale(dist(tdata),k=3)
pc<-prcomp(tdata)
pc3d<-cbind(pc$x[,1], pc$x[,2], pc$x[,3])
fname="Overview Figure/PCA3D.tiff"
#tiff(filename =fname, units="in", width=10, height=10, res=150)
#scatter3D(pc$x[,1], pc$x[,2], pc$x[,3], phi = 0, bty ="g")
#dev.off()

moviename="PCA3D"
if (length(unique(group))<8&length(unique(group))>2) {
colset=as.data.frame(cbind(brewer.pal(n = length(unique(group)), name = "Accent"),unique(group)))
color2=colset[match(group,colset[,2]),1]
 } else if (length(unique(group))<=2) {
 colset=as.data.frame(cbind(brewer.pal(n = 3, name = "Accent")[1:length(unique(group))],unique(group)))
 color2=colset[match(group,colset[,2]),1]
 } else {
 colset=as.data.frame(cbind(unique(group),unique(group)))
 color2=colset[match(group,colset[,2]),1]}

par3d(windowRect = c(0,50, 800, 800))
legend3d("topright",legend = colset[,2][order(colset[,2])],pch = 16, col =colset[,1][order(colset[,2])], cex=2, inset=c(0.03))
plot3d(pc3d,col=color2, type="s",size=2,scale=0.2,pch=20,cex.lab=1,resfac=3,xlab="PC1",ylab="PC2",zlab="PC3")
text3d(pc$x[,1], pc$x[,2], pc$x[,3],row.names(pc$x),adj=1.2,family="sans",cex=1)
movie3d(spin3d(axis=c(0,0,1)), duration=7, fps=10, movie = moviename,dir=paste0(getwd(),"/Overview Figure"),type = "gif",clean=T)
}

else if (type=="histogram") {
xmin=floor(min(forgraph$Log2values))-1
xmax=10*ceiling(max(forgraph$Log2values)/10)
ymax=0
for (i in 1:ncol(x)) {
v=forgraph[forgraph$variable==unique(forgraph$variable)[i],]$Log2values
if (max(hist(v,breaks=seq(xmin,xmax,by=0.1),plot=F)$counts)>ymax) {
ymax=max(hist(v,breaks=seq(xmin,xmax,by=0.1),plot=F)$counts)
}}
for (i in 1:ncol(x)) {
graphi=
ggplot(forgraph, aes(x=Log2values)) +
	labs(title=paste("Log 2 Peptide Intensities of ",colnames(x)[i],collapse="")) +
	xlim(xmin,xmax) +	
	ylim(0,ymax) +	
	geom_histogram(data=subset(forgraph,variable==unique(variable)[i]),fill=(i+1), alpha=1, binwidth = 0.1) 
hiss[[length(hiss)+1]]=graphi}
for (i in 0:floor((length(hiss)-1)/6)){
a=hiss[[(i*6+1)]]
b=tryCatch(hiss[[(i*6+2)]],error=function(err) "")
c=tryCatch(hiss[[(i*6+3)]],error=function(err) "")
d=tryCatch(hiss[[(i*6+4)]],error=function(err) "")
e=tryCatch(hiss[[(i*6+5)]],error=function(err) "")
f=tryCatch(hiss[[(i*6+6)]],error=function(err) "")
source("http://peterhaschke.com/Code/multiplot.R")
jpeg(filename = paste0("Overview Figure/histogram",i+1,".jpg"), width = 840, height = 630,pointsize =12, quality = 500, bg = "white", res = NA, restoreConsole = TRUE)
multiplot(a, b, c, d, e, f, cols=2)
dev.off()
}
}

else if (type=="boxplot") {
o <- function(y) {
  subset(y, y < quantile(y,probs=c(0.08))[1] | quantile(y,probs=c(0.92))[1] < y)
}
bpd=ggplot(forgraph, aes(x=variable, y=Log2values, colour=factor(variable))) 
bpd=bpd+geom_boxplot()
bpd=bpd+stat_summary(fun.y=o, geom="point", aes(colour=factor(variable)))+labs(title="Boxplot of Log 2 Intensity Distributions of All Samples", x="Sample Label", y="Intensity") 
bpd=bpd+theme(legend.title=element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
wd=12+(length(unique(forgraph$variable))-10)*0.2
ht=8
fname="Overview Figure/Overview_boxplot.tiff"
tiff(filename =fname, units="in", width=wd, height=ht, res=300)
print(bpd)
dev.off()

}

else {stop("invalid variable")}
}

grapov(rawdatafile[,7:ncol(rawdatafile)],type="heatmap",group=gcolor)
grapov(rawdatafile[,7:ncol(rawdatafile)],type="PCA",group=gcolor)
grapov(rawdatafile[,7:ncol(rawdatafile)],type="PCA3d",group=gcolor)
grapov(rawdatafile[,7:ncol(rawdatafile)],type="histogram",group=gcolor)
grapov(rawdatafile[,7:ncol(rawdatafile)],type="boxplot",group=gcolor)
