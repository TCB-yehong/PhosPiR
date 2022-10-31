#######################download all libraries###################################
if(!require(clusterProfiler)){BiocManager::install("clusterProfiler",update=F,ask=F)}
library(clusterProfiler)
options(clusterProfiler.download.method = "wininet")
if(!require(pathview)){BiocManager::install("pathview",update=F,ask=F)}
library(pathview)
if(!require(vroom)){install.packages("vroom")}
library(vroom)
if(!require(textreadr)){install.packages("textreadr")}
library(textreadr)
################################################################################

#######################enrichment run###########################################
#check if organism has annotation package available, if it does, proceed to enrichment analysis
dbinfo=data.frame(organism=c("Homo sapiens","Mus musculus","Rattus norvegicus","Drosophila melanogaster",
"Saccharomyces cerevisiae","Arabidopsis thaliana","Danio rerio","Caenorhabditis elegans","Bos taurus",
"Gallus gallus","Macaca mulatta","Sus scrofa","Canis familiaris","Escherichia coli","Xenopus laevis",
"Anopheles gambiae","Pan troglodytes","Plasmodium falciparum"),database=c("org.Hs.eg.db","org.Mm.eg.db",
"org.Rn.eg.db","org.Dm.eg.db","org.Sc.sgd.db","org.At.tair.db","org.Dr.eg.db","org.Ce.eg.db","org.Bt.eg.db",
"org.Gg.eg.db","org.Mmu.eg.db","org.Ss.eg.db","org.Cf.eg.db","org.EcK12.eg.db","org.Xl.eg.db","org.Ag.eg.db",
"org.Pt.eg.db","org.Pf.plasmo.db"),code3=c("hsa","mmu","rno","dme","sce","ath","dre","cel","bta","gga","mcc",
"ssc","cfa","eco","xla","aga","ptr","pfa"))
phylocheck=match(phylo,dbinfo$organism)
if (is.na(phylocheck)) {
msgBox("Organism not supported for enrichment analysis. Moving on!")
} else {
msgBox("Now proceeding to Enrichment Analysis.")
dir.create("Enrichment",showWarnings=F)
#download species database
genodb=dbinfo[phylocheck,2]
if(!require(genodb,character.only=T)){BiocManager::install(genodb,update=F,ask=F)}
library(genodb,character.only=T)

#set up background data & datalist data
allgeneid=statsRes$entrez_id
#read sig lists
nwgenes=list()
n=0
inlist=c()
for (i in siglistFiles) {
n=n+1
nwgenestemp=read.csv(paste0("Statistical Analysis\\Significant Lists\\",i),stringsAsFactors=F)
nwgenes[[n]]=nwgenestemp
names(nwgenes)[[n]]=i}
#create gene list and background list 
genesets=list()
allgeneset=list()
#for all sig lists:
#check if fold change is used as cutoff and thus recorded
#if yes, extract from the entire dataset the specific column in dataset this particular comparison's FC is located, combine with gene symbols
#if no, find any pvalues in the cutoff, and do the same as FC 
#clean out NA, get ENTREZ ID for all background genes 
for (i in 1:length(nwgenes)) {
if (sum(grepl("FoldChange",colnames(nwgenes[[i]])))>0) {
allgeneod=order(statsRes[,intersect(grep(strsplit(names(nwgenes)[[i]],"_")[[1]][1],colnames(statsRes)),grep("_FoldChange",colnames(statsRes)))[1]],decreasing=T)
tempallgeneset=statsRes[allgeneod,intersect(grep(strsplit(names(nwgenes)[[i]],"_")[[1]][1],colnames(statsRes)),grep("_FoldChange",colnames(statsRes)))[1]]
listfccol=grep("FoldChange",colnames(nwgenes[[i]]))} else {
allgeneod=order(statsRes[,intersect(grep(strsplit(names(nwgenes)[[i]],"_")[[1]][1],colnames(statsRes)),grep("Pvalue|FDR|PFP",colnames(statsRes)))[1]],decreasing=T)
tempallgeneset=statsRes[allgeneod,intersect(grep(strsplit(names(nwgenes)[[i]],"_")[[1]][1],colnames(statsRes)),grep("Pvalue|FDR|PFP",colnames(statsRes)))[1]]
listfccol=grep("Pvalue|FDR|PFP",colnames(nwgenes[[i]]))}

names(tempallgeneset)=allgeneid[allgeneod]
tempallgeneset=tempallgeneset[which(!is.na(names(tempallgeneset)))]
listgeneid=statsRes$entrez_id[match(nwgenes[[i]]$Protein,statsRes$Protein)]

#do the same with sig genes, this time extract column from sig list, not entire dataset 
if (length(listfccol)>1) {
listfccol=listfccol[1]}
if (length(listfccol)!=0){
inlist=c(inlist,i)
listgeneod=order(nwgenes[[i]][,listfccol],decreasing=T)
templistgene=nwgenes[[i]][,listfccol][listgeneod]
names(templistgene)=listgeneid[listgeneod]
templistgene=templistgene[which(!is.na(names(templistgene)))]
#record in corresponding lists 
genesets[[i]]=templistgene
allgeneset[[i]]=tempallgeneset
}
} 
#name each item in the list by their comparison and cutoff number to distinguish the sig/background list origin 
names(genesets)=paste(unlist(lapply(names(nwgenes)[inlist],FUN=function(x) substring(x,1,gregexpr("\\_",x)[[1]][2]-1))),
"significant","genes",sep="_")
names(allgeneset)=paste(unlist(lapply(names(nwgenes)[inlist],FUN=function(x) substring(x,1,gregexpr("\\_",x)[[1]][2]-1))),
"significant","genes",sep="_")

###KEGG enrichment###
dir.create("Enrichment/KEGG enrichment",showWarnings=F)
#wddflt=getwd()
setwd(paste0(wddflt,"/Enrichment/KEGG enrichment"))
#function for making an error notice plot for when connection is unstable and pathview plots are giving errors
textPlot <- function(plotname, string){
  par(mar=c(0,0,0,0))
  png(paste0(plotname, ".png"))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(string), cex = 1, col = "black", family="serif", font=2, adj=0.5)
  dev.off()
	dev.off()
}
#create excel file for organizing all pathway figure links
plink=createWorkbook("PathviewLinks")
addWorksheet(plink,"Sheet1")
n=0

#setup pathview color parameter for different range of FC/pval
for (i in 1:length(genesets)) {
if (max(genesets[[i]])<=1&min(genesets[[i]])>=0) {
midcolor=list(gene = "seashell", cpd = "seashell")
highcolor=list(gene = "grey", cpd = "grey")
lowcolor=list(gene = "red", cpd = "orange")
pvlimits=list(gene=c(0,0.1),cpd=c(0,1))
bothdir=list(gene=TRUE, cpd=TRUE)
pvbins = list(gene = 10, cpd= 10)
} else if (sum(abs(genesets[[i]])==Inf)>0) {
midcolor=list(gene = "white", cpd = "white")
highcolor=list(gene = "red", cpd = "orange")
lowcolor=list(gene = "blue", cpd = "purple")
pvlimits=list(gene=c(-1,1),cpd=c(-1,1))
bothdir=list(gene=TRUE, cpd=TRUE)
pvbins = list(gene = 10, cpd= 10)
} else {
midcolor=list(gene = "white", cpd = "white")
highcolor=list(gene = "red", cpd = "orange")
lowcolor=list(gene = "blue", cpd = "purple")
pvlimits=list(gene=c(min(genesets[[i]]),max(genesets[[i]])),cpd=c(-1,1))
bothdir=list(gene=TRUE, cpd=TRUE)
pvbins = list(gene = 10, cpd= 10)
}
#run enrichment analysis with entire dataset as background and with universal background 
kk2=enrichKEGG(gene = names(genesets[[i]]), 
			   universe=names(allgeneset[[i]]),
               organism = dbinfo$code3[phylocheck],
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.25,
               pAdjustMethod = "BH",
               minGSSize = 10,
               maxGSSize = 500)
#kk2=setReadable(kk2, eval(as.name(genodb)))
kk3=enrichKEGG(gene = names(genesets[[i]]), 
#			   universe=names(allgeneset[[i]]),
               organism = dbinfo$code3[phylocheck],
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.25,
               pAdjustMethod = "BH",
               minGSSize = 10,
               maxGSSize = 500)
#kk3=setReadable(kk3, eval(as.name(genodb)))

#if there are any enriched pathways in dataset background enrichment:
#create dotplot in tiff format, and record results in csv format
if (!is.null(nrow(kk2))) {	   
if (nrow(kk2)>0) {
n=n+1
#auto adjust figure size based on number of sig pathway hits (formula created through repeated testing specfically for this)
ht=5+log(nrow(kk2),base=1.4)
wd=2.5+max(nchar(kk2$Description))/4+(15/(5+log(nrow(kk2),base=1.22)))
#setup tiff parameters and create tiff 
fname=paste0(names(genesets)[[i]],"_KEGG_enrichment_dotplot_bgdata.tiff")
ftitle=paste0(names(genesets)[[i]],"_dataBackground")
catcount=nrow(kk2)
tiff(filename =fname, units="in", width=wd, height=ht, res=100)
print(dotplot(kk2,showCategory = catcount,title=ftitle))
dev.off()
#create csv file to record result 
fname2=paste0(names(genesets)[[i]],"_dataBackground_KEGGenrichmentResult.csv")
write.csv(kk2,file=fname2)

#create pathview plot for each sig pathway,then organize all the figure link in excel 
for (j in kk2$ID) {
#pathview function, with error function embedded. Upon error, an error notice figure will be made
dme <- tryCatch(pathview(gene.data=genesets[[i]], pathway.id=j, species = dbinfo$code3[phylocheck],
out.suffix=paste0("from_",names(genesets)[[i]],"_enrichment_dataBg"),both.dirs=bothdir,bins=pvbins,
limit=pvlimits,low = lowcolor,mid=midcolor,high=highcolor),
error=function(e) {textPlot(paste0(j,".from_",names(genesets)[[i]],"_enrichment_dataBg"),"Pathway error from KEGG database!")})
}
#create excel entries for the plotted pathways, with description as text and hyperlink to the plots 
kklink=paste0(kk2$ID,".from_",names(genesets)[[i]],"_enrichment_dataBg.png")
names(kklink)=kk2$Description
class(kklink)="hyperlink"
#saving in excel, each enrichment run takes 1 column, first row is descrition of the run, rest rows are pathway name and link to the pathview plot of the pathway 
c1=paste0(names(genesets)[[i]],"_KEGGenrichment_dataBackground")
writeData(plink,sheet=1,x=c1,startCol=n,startRow=1)
writeData(plink,sheet=1,x=kklink,startCol=n,startRow=2)
}}

#repeat of above functions, but for universal analysis enrichment 
if (!is.null(nrow(kk3))) {	  
if (nrow(kk3)>0) {
n=n+1
ht=3+log(nrow(kk3),base=1.4)
wd=2.5+max(nchar(kk3$Description))/4+(15/(5+log(nrow(kk3),base=1.22)))
fname=paste0(names(genesets)[[i]],"_KEGG_enrichment_dotplot_bgall.tiff")
ftitle=paste0(names(genesets)[[i]],"_universalBackground")
catcount=nrow(kk3)
tiff(filename =fname, units="in", width=wd, height=ht, res=100)
print(dotplot(kk3,showCategory = catcount,title=ftitle))
dev.off()
fname2=paste0(names(genesets)[[i]],"_universalBackground_KEGGenrichmentResult.csv")
write.csv(kk3,file=fname2)

for (j in kk3$ID) {
dme <- tryCatch(pathview(gene.data=genesets[[i]], pathway.id=j, species = dbinfo$code3[phylocheck],
out.suffix=paste0("from_",names(genesets)[[i]],"_enrichment_allBg"),both.dirs=bothdir,bins=pvbins,
limit=pvlimits,low = lowcolor,mid=midcolor,high=highcolor),
#limit=list(gene=max(abs(genesets[[i]])),cpd=1),low = list(gene = "blue", cpd = "purple")),
error=function(e) {textPlot(paste0(j,".from_",names(genesets)[[i]],"_enrichment_allBg"),"Pathway error from KEGG database!")})
}
kklink=paste0(kk3$ID,".from_",names(genesets)[[i]],"_enrichment_allBg.png")
names(kklink)=kk3$Description
class(kklink)="hyperlink"

c1=paste0(names(genesets)[[i]],"_KEGGenrichment_universalBackground")
writeData(plink,sheet=1,x=c1,startCol=n,startRow=1)
writeData(plink,sheet=1,x=kklink,startCol=n,startRow=2)
}}
}
saveWorkbook(plink,file="KEGG_enrichedPathways_pathview_link.xlsx",overwrite=T)

###GO enrichment###
setwd(wddflt)
dir.create("Enrichment/GO enrichment",showWarnings=F)
setwd(paste0(wddflt,"/Enrichment/GO enrichment"))
#performs GO enrichment and creates dotplot and csv record of results whenever there are results
#See KEGG enrichment section for specific code descriptions  
for (i in 1:length(genesets)) {
	for (j in c("CC","MF","BP")) {
	kk2 <- enrichGO(gene = names(genesets[[i]]), 
                universe      = names(allgeneset[[i]]),
                OrgDb         = eval(as.name(genodb)),
                ont           = j,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
	kk3 <- enrichGO(gene = names(genesets[[i]]), 
            #    universe      = names(allgeneset[[i]]),
                OrgDb         = eval(as.name(genodb)),
                ont           = j,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
		
if (!is.null(nrow(kk2))) {		
	if (nrow(kk2)>0) {
		ht=5+log(nrow(kk2),base=1.4)+1.016^nrow(kk3)
		wd=3.5+max(nchar(kk2$Description))/4+(15/(5+log(nrow(kk2),base=1.22)))
		fname=paste0(names(genesets)[[i]],"_GO_",j,"_enrichment_dotplot_bgdata.tiff")
		ftitle=paste0(names(genesets)[[i]],"_GO_",j,"_dataBackground")
		catcount=nrow(kk2)
		tiff(filename =fname, units="in", width=wd, height=ht, res=100)
		print(dotplot(kk2,showCategory = catcount,title=ftitle))
		dev.off()
		fname2=paste0(names(genesets)[[i]],"_dataBackground_GO_",j,"_enrichmentResult.csv")
		write.csv(kk2,file=fname2)}}
if (!is.null(nrow(kk3))) {			
	if (nrow(kk3)>0) {
		ht=3+log(nrow(kk3),base=1.4)+1.016^nrow(kk3)
		wd=3.5+max(nchar(kk3$Description))/4+(15/(5+log(nrow(kk3),base=1.22)))
		fname=paste0(names(genesets)[[i]],"_GO_",j,"_enrichment_dotplot_bgall.tiff")
		ftitle=paste0(names(genesets)[[i]],"_GO_",j,"_universalBackground")
		catcount=nrow(kk3)
		tiff(filename =fname, units="in", width=wd, height=ht, res=100)
		print(dotplot(kk3,showCategory = catcount,title=ftitle))
		dev.off()
		fname2=paste0(names(genesets)[[i]],"_universalBackground_GO_",j,"_enrichmentResult.csv")
		write.csv(kk3,file=fname2)}}
		}}

###DOSE enrichment###
#performs DOSE enrichment if organism is human, and creates dotplot and csv record of results whenever there are results
#See KEGG enrichment section for specific code descriptions 
if (phylo=="Homo sapiens") {
setwd(wddflt)
dir.create("Enrichment/Disease Association enrichment",showWarnings=F)
setwd(paste0(wddflt,"/Enrichment/Disease Association enrichment"))
library(DOSE)
for (i in 1:length(genesets)) {
kk2 <- enrichDO(gene          = names(genesets[[i]]), 
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(allgeneset[[i]]),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = TRUE)
kk3 <- enrichDO(gene          = names(genesets[[i]]), 
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
#              universe      = names(allgeneset[[i]]),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = TRUE)

if (!is.null(nrow(kk3))) {	
	if (nrow(kk2)>0) {
		ht=5+log(nrow(kk2),base=1.4)+1.016^nrow(kk3)
		wd=3+max(nchar(kk2$Description))/4+(15/(5+log(nrow(kk2),base=1.22)))
		fname=paste0(names(genesets)[[i]],"_DOSE_enrichment_dotplot_bgdata.tiff")
		ftitle=paste0(names(genesets)[[i]],"_DOSE_dataBackground")
		catcount=nrow(kk2)
		tiff(filename =fname, units="in", width=wd, height=ht, res=100)
		print(dotplot(kk2,showCategory = catcount,title=ftitle))
		dev.off()
		fname2=paste0(names(genesets)[[i]],"_dataBackground_DOSE_enrichmentResult.csv")
		write.csv(kk2,file=fname2)}}
if (!is.null(nrow(kk3))) {			
	if (nrow(kk3)>0) {
		ht=3+log(nrow(kk3),base=1.4)+1.016^nrow(kk3)
		wd=3+max(nchar(kk3$Description))/4+(15/(5+log(nrow(kk3),base=1.22)))
		fname=paste0(names(genesets)[[i]],"_DOSE_enrichment_dotplot_bgall.tiff")
		ftitle=paste0(names(genesets)[[i]],"_DOSE_universalBackground")
		catcount=nrow(kk3)
		tiff(filename =fname, units="in", width=wd, height=ht, res=100)
		print(dotplot(kk3,showCategory = catcount,title=ftitle))
		dev.off()
		fname2=paste0(names(genesets)[[i]],"_universalBackground_DOSE_enrichmentResult.csv")
		write.csv(kk3,file=fname2)}}
		}}

###Wikipathway enrichment###
#Wikipathway enrichment is available for a selected set of organisms, if there are no results, nothing will be recorded
#See KEGG enrichment section for specific code descriptions 
if (phylo %in% c("Sus scrofa","Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
"Caenorhabditis elegans","Canis familiaris","Danio rerio","Drosophila melanogaster","Gallus gallus",
"Pan troglodytes","Rattus norvegicus","Saccharomyces cerevisiae")) {
setwd(wddflt)
dir.create("Enrichment/Wikipathway enrichment",showWarnings=F)
setwd(paste0(wddflt,"/Enrichment/Wikipathway enrichment"))
URLwikidir="http://data.wikipathways.org/current/gmt/"
wikilist=read_html(URLwikidir)
wikidbname=wikilist[grep(gsub(" ","_",phylo),wikilist)]
download.file(paste0(URLwikidir,wikidbname),wikidbname)

wp2gene <- read.gmt(wikidbname)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) 
wpid2name <- wp2gene %>% dplyr::select(wpid, name) 

for (i in 1:length(genesets)) {
kk2 <- enricher(names(genesets[[i]]), TERM2GENE = wpid2gene, TERM2NAME = wpid2name,minGSSize = 5)
	if (!is.null(nrow(kk2))) {
		if (nrow(kk2)>0) {
		ht=5+log(nrow(kk2),base=1.4)+1.016^nrow(kk3)
		wd=3+max(nchar(kk2$Description))/4+(15/(5+log(nrow(kk2),base=1.22)))
		fname=paste0(names(genesets)[[i]],"_Wikipathway_enrichment_dotplot.tiff")
		ftitle=paste0(names(genesets)[[i]],"_Wikipathway Enrichment")
		catcount=nrow(kk2)
		tiff(filename =fname, units="in", width=wd, height=ht, res=100)
		print(dotplot(kk2,showCategory = catcount,title=ftitle))
		dev.off()
		fname2=paste0(names(genesets)[[i]],"_Wikipathway_enrichmentResult.csv")
		write.csv(kk2,file=fname2)}}
}
if(length(dir(all.files=TRUE))<=3) {
setwd(wddflt)
unlink("Enrichment/Wikipathway enrichment",recursive=T)}
}

###cellmarker enrichment###
#performs cell marker enrichment if organism is human or mouse, and creates dotplot and csv record of results whenever there are results
#See KEGG enrichment section for specific code descriptions 
if (phylo=="Homo sapiens"|phylo=="Mus musculus"){
setwd(wddflt)
dir.create("Enrichment/Cell Marker enrichment",showWarnings=F)
setwd(paste0(wddflt,"/Enrichment/Cell Marker enrichment"))
if (phylo=="Homo sapiens") {
URL="http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt"
} else {
URL="http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt"
}
cell_markers <- vroom::vroom(URL) %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
   dplyr::select(cellMarker, geneID) %>%
   dplyr::mutate(geneID = strsplit(geneID, ', '))
   
for (i in 1:length(genesets)) {
kk2 <- enricher(names(genesets[[i]]), TERM2GENE=cell_markers, minGSSize=1)
	if (!is.null(nrow(kk2))) {
		if (nrow(kk2)>0) {
		ht=5+log(nrow(kk2),base=1.4)+1.016^nrow(kk3)
		wd=3+max(nchar(kk2$Description))/4+(15/(5+log(nrow(kk2),base=1.22)))
		fname=paste0(names(genesets)[[i]],"_CellMarker_enrichment_dotplot.tiff")
		ftitle=paste0(names(genesets)[[i]],"_Cell Marker Enrichment")
		catcount=nrow(kk2)
		tiff(filename =fname, units="in", width=wd, height=ht, res=100)
		print(dotplot(kk2,showCategory = catcount,title=ftitle))
		dev.off()
		fname2=paste0(names(genesets)[[i]],"_CellMarker_enrichmentResult.csv")
		write.csv(kk2,file=fname2)}}
}
if(length(dir(all.files=TRUE))<=2) {
setwd(wddflt)
unlink("Enrichment/Cell Marker enrichment",recursive=T)}
}
setwd(wddflt)
}
################################################################################