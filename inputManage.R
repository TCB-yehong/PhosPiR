#######################download all libraries###################################
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#	BiocManager::install(version = BiocManager::version(),ask=F,update=F)
if(!require(htmltools)){BiocManager::install("htmltools",update=F,ask=F)}
library(htmltools)
if(!require(msImpute)){BiocManager::install("msImpute",update=F,ask=F)}
library(msImpute)
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
if(!require(openxlsx)){install.packages("openxlsx")}
library(openxlsx)
if(!require(remotes)){install.packages("remotes")}
library(remotes)
if(!require(proBatch)){remotes::install_github("symbioticMe/proBatch", upgrade='never')}
library(proBatch)
#bioc_deps <- c("GO.db", "impute", "preprocessCore", "pvca","sva" )
#cran_deps <- c("corrplot", "data.table", "ggplot2", "ggfortify","lazyeval", "pheatmap", "reshape2", "rlang", 
#               "tibble", "dplyr", "tidyr", "wesanderson","WGCNA") 
#for (i in bioc_deps) {
#if(!require(i, character.only = TRUE)){BiocManager::install(i,update=F,ask=F)}
#}
#for (i in cran_deps) {
#if(!require(i, character.only = TRUE)){install.packages(i,update=F,ask=F)}
#}
#download.file("https://github.com/symbioticMe/proBatch/archive/refs/heads/master.zip", 
#              destfile = "proBatch-master.zip")
#unzip("proBatch-master.zip")
#install.packages("proBatch-master", repos = NULL, type = "source")
################################################################################

#######################deposit location#########################################
setwd(dlg_dir(default = getwd(),title="Where should results be stored?")$res)
fdn=paste0("Phosphoproteomics_analysis_",Sys.Date())
#if location exists, then add suffix of dash and the next number for the new result folder
if (!dir.exists(fdn)) {
dir.create(fdn) } else {
n=1
while (dir.exists(paste0(fdn,"-",n))) {
n=n+1}
dir.create(paste0(fdn,"-",n))
fdn=paste0(fdn,"-",n)}
setwd(fdn)
wddflt=getwd()
################################################################################

#######################input type###############################################
fileType=dlgList(c("MaxQaunt","Other"),title="Choose input file source")$res
if (fileType=="MaxQaunt") {
fname=dlg_open(title = "Select MaxQuant result file \'Phospho(STY)Sites\' for input",multiple=F)$res
annotpep=NULL
annotsamp=NULL
} else {
fname=dlg_open(title = "Select input file with correct format (.csv or .xlsx)",multiple=F)$res
#add annotation file input option if file is not from maxquant
annotFA=dlgList(c("No annotation files","Row (peptide) annotation available",
"Column (sample) annotation available","Both annotations available"),title="Add annotation files?")$res
if (annotFA=="Both annotations available") {
annotpep=dlg_open(title = "Select row (peptide) annotation file (.csv)",multiple=F)$res
annotsamp=dlg_open(title = "Select column (sample) annotation file (.csv)",multiple=F)$res
} else if (annotFA=="Row (peptide) annotation available"){
annotpep=dlg_open(title = "Select row (peptide) annotation file (.csv)",multiple=F)$res
annotsamp=NULL
} else if (annotFA=="Column (sample) annotation available"){
annotpep=NULL
annotsamp=dlg_open(title = "Select column (sample) annotation file (.csv)",multiple=F)$res
} else {
annotpep=NULL
annotsamp=NULL
}
}
################################################################################


#######################run option###############################################
#setup which step to run
if (fileType=="MaxQaunt") {
pipestep=dlgList(c("Generate Special Input File","Annotation","Overview Figures","Differential Analysis","All Steps"),
title="Select which analysis to run")$res
} else {
pipestep=dlgList(c("Annotation","Overview Figures","Differential Analysis","All steps"),
title="Select which analysis to run")$res}
################################################################################

#######################maxquant input formatting################################
if (pipestep!="Generate Special Input File") {
if (fileType=="MaxQaunt") {
raw=read.delim(fname)
#2 files are created from maxquant raw file, one for intensity and one for ID info

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

#ID info file:
#extract available id and info columns by column number
iidx=match(c("Protein","Protein.names","Gene.names","Amino.acid","Position",
"Sequence.window","Proteins","Positions.within.proteins","Fasta.headers",
"Localization.prob","Number.of.Phospho..STY.","Modification.window",
"Phospho..STY..Probabilities","Position.in.peptide","id","Protein.group.IDs",
"Positions","Peptide.IDs","Mod..peptide.IDs","Evidence.IDs","MS.MS.IDs"),
colnames(raw))
rawinfofile=raw[,iidx]
#create new columns for each sample recording % of each phosphosite count's
#intensity contributing to the combined intensity
#extract columns with separated intensity info (by # of phosphorylation)
intsep=raw[,intidx[(intidx2[2]+1):length(intidx)]]
#check how many phosphorylation counts the intensity is separated into
intsepgrp=length(table(unlist(lapply(colnames(intsep),FUN=function(x) strsplit(x,"__")[[1]][2]))))
#record column number of previously extracted columns
rawinfofilecoln=ncol(rawinfofile)
#for each group of intensities (grouped by sample), intensity is divided by
#sum intensity, then multiplied by 100 for percentage
#all percentage from 1 group are combined together separated by semicolon
#and added to the id info file
for (i in seq(1,ncol(intsep),intsepgrp)) {
temperc=intsep[,i:(i+intsepgrp-1)]/rowSums(intsep[,i:(i+intsepgrp-1)])*100
rawinfofile[,ncol(rawinfofile)+1]=apply(temperc,1,paste,collapse=";")
}
#add proper names for these columns
colnames(rawinfofile)[(rawinfofilecoln+1):ncol(rawinfofile)]=
paste0(colnames(raw)[didx2],"_percentIntensityFrom_phosphoCount",paste(names(
table(unlist(lapply(colnames(intsep),FUN=function(x) strsplit(x,"__")[[1]][2])))),collapse=""))

#save both files just in case
dir.create("Raw File")
write.csv(rawdatafile,"Raw File/maxquant_raw_data.csv",row.names=F)
write.csv(rawinfofile,"Raw File/maxquant_raw_info.csv",row.names=F)

################################################################################

#######################refine matrix############################################
#remove reverse and contaminant entries
rawdatafile$Reverse[is.na(rawdatafile$Reverse)]=""
rawdatafile$Potential.contaminant[is.na(rawdatafile$Potential.contaminant)]=""
row.filt <- which(rawdatafile$Reverse==""&rawdatafile$Potential.contaminant=="")                    
rawdatafile=rawdatafile[row.filt,-c(7,8)]
rawinfofile=rawinfofile[row.filt,]
annotpepdf=NULL
annotsampdf=NULL
} else {
#read input file from 'Other' based on whether it's a csv or xlsx file 
fnameparts=strsplit(fname,"\\.")[[1]]
if (fnameparts[length(fnameparts)]=="csv") {
rawdatafile=read.csv(fname,stringsAsFactors=F)
} else if (fnameparts[length(fnameparts)]=="xlsx") {
rawdatafile=read.xlsx(fname,sheet=1,check.names=T)
} else {
dlg_message("Only .csv and .xslx formats are supported!")
}

if (!is.null(annotpep)) {
annotpepparts=strsplit(annotpep,"\\.")[[1]]
if (annotpepparts[length(annotpepparts)]=="csv") {
annotpepdf=read.csv(annotpep,stringsAsFactors=F)
} else {
dlg_message("Only .csv format is supported for row annotation file!")
annotpepdf=NULL
} } else {annotpepdf=NULL}

if (!is.null(annotsamp)) {
annotsampparts=strsplit(annotsamp,"\\.")[[1]]
if (annotsampparts[length(annotsampparts)]=="csv") {
annotsampdf=read.csv(annotsamp,stringsAsFactors=F)
} else {
dlg_message("Only .csv format is supported for column annotation file!")
annotsampdf=NULL
} } else {annotsampdf=NULL}


colnames(rawdatafile)[1:6]=c("Protein","Protein.names","Gene.names","Amino.acid","Position","Sequence.window")
msg_box("Please make sure the file follows the described
format, and have reverse sequence and potential 
contaminants removed to ensure accurate analysis. 
Thank you!")
}

#clean gene names
genecleaned=unlist(lapply(rawdatafile$Gene.names ,FUN=function(x) {
cleaned=c()
if (!is.na(x)) {
temp=strsplit(x,";|,|, ")
for (k in 1:length(temp)) {
cleaned=c(cleaned,temp[[k]][1])
} } else {cleaned=c(cleaned,NA)}
return(cleaned)}))
rawdatafile$Gene.names=genecleaned

#internal check row uniqueness
ucheck=sum(table(paste(rawdatafile$Protein,rawdatafile$Position,sep="_"))!=1) 
if (ucheck>0) {
msg_box("Please note nonunique sites (rows) exist, 
merging nonunique sites is recommended before
running the pipeline.")}

#look for outliers 
rintfile=rawdatafile[,7:ncol(rawdatafile)]
rintfile[is.na(rintfile)]=0
#remove rows with too many missing value
rrgo=dlgMessage("Remove rows with too many missing values? (Yes is recommended)",
type="yesno")$res
if (rrgo=="yes") {
rcoff=dlgList(c("Remove rows with missing value >50%","Remove rows with missing value >30%",
"Remove rows with missing value >10%"),title="Choose missing value cutoff")$res
if (rcoff=="Remove rows with missing value >50%") {rcoffv=0.5
} else if (rcoff=="Remove rows with missing value >30%") {rcoffv=0.3
} else {rcoffv=0.1}
rpass=apply(rintfile,1,FUN=function(x) (sum(x==0)/length(x))<rcoffv)
rintfile=rintfile[which(rpass),]
rawdatafile=rawdatafile[which(rpass),]
if (!is.null(annotpepdf)) {
annotpepdf=annotpepdf[which(rpass),]
}
if (exists("rawinfofile")) {
rawinfofile=rawinfofile[which(rpass),]}
}
#remove columns with too many missing value
crgo=dlgMessage("Remove samples with too many missing values? (Not 
recommended unless outliers are expected, i.e. possible contaminations)",type="yesno")$res
if (crgo=="yes") {
cnasum=apply(rintfile,2,FUN=function(x) sum(x==0))
rintfile=rintfile[,which(cnasum<=mean(cnasum)+2*sd(cnasum))]
rawdatafile=rawdatafile[,c(1:6,which(cnasum<=mean(cnasum)+2*sd(cnasum))+6)]
if (exists("rawinfofile")) {
rawinfofile=rawinfofile[,c(1:rawinfofilecoln,which(cnasum<=mean(cnasum)+2*sd(cnasum))+rawinfofilecoln)]}
}
}
################################################################################

#######################universal format#########################################
#has more to do with non-maxquant input when it's not in the right format
#at the moment this is user resposibility
#will implement in newer versions some build-in data formatting to help user 
#transform their data to a usable format for the pipeline
################################################################################

#######################user setup options#######################################
#setup organism
if (pipestep!="Generate Special Input File"&pipestep!="Overview Figures") {
phylo=dlgList(c("Homo sapiens","Mus musculus","Rattus norvegicus",
"Saccharomyces cerevisiae","Caenorhabditis elegans","Danio rerio",
"Drosophila melanogaster","Escherichia coli","Arabidopsis thaliana",
"Other"),preselect="Homo sapiens",title="Select data organism source")$res
if(phylo=="Other") {phylo=dlgInput("Please enter organism scientific name: (Capitalize first 
letter and nothing else, space between words)")$res}}

#setup imputation and normalization
if (pipestep!="Generate Special Input File"&pipestep!="Annotation") {
impnorm=dlgList(c("Neither","Only normalize","Both normalize and impute"),
title="Imputation and normalization")$res}
################################################################################

#######################group setup##############################################
#perform group setup from groupSetup.R
if (pipestep!="Generate Special Input File"&pipestep!="Annotation") {
source(paste0(codepth,"groupSetup.R"))
}
################################################################################

#######################normalization and imputation#############################
if (pipestep!="Generate Special Input File"&pipestep!="Annotation") {
msgBox("Now proceeding to normalization and imputation.")
if (impnorm=="Both normalize and impute") {
#dlgMessage("Caution with imputation, if there's time, another run without 
#imputation is recommended as control")

rintfile[rintfile==0]=NA
#quantile normalize dataset
quantile_normalized_matrix = normalize_data_dm(as.matrix(rintfile),normalize_func = 'quantile',log_base = 2, offset = 1)

#if there are rows with less than 4 non-NA values, imputation cannot be done, user is forced to choose between neither or only normalize 
if (min(unlist(apply(quantile_normalized_matrix,1,FUN=function(x) sum(!is.na(x)))))<4) {
dlgMessage("Not all rows has 4 or more non-NA values,Imputation cannot be done. Please try row filtering in the next run if imputation is insisted.")
impnorm=character()
while (length(impnorm)==0) {
impnorm=dlgList(c("Neither","Only normalize"),
title="Imputation and normalization")$res}
} else {

impgrpyes=dlgMessage("Should group label be considered for imputation?

(i.e. genotype 1 missing values is part of genotype 1 value 
distribution, and not part of genotype 2 value distribution,
instead of genotype 1 missing values is part of all genotype
value distribution.)",type="yesno")$res
if (impgrpyes=="yes") {
impgrp=character()
while (length(impgrp)==0) {
impgrp=dlgList(names(group),title="Which group label?")$res}
#get group order info for selected grouping method 
impgrpo=as.factor(group[[which(names(group)==impgrp)]])
#perform imputation 
rtransifile=msImpute(quantile_normalized_matrix,method="v2-mnar",group=impgrpo)
} else {
#if (min(unlist(apply(quantile_normalized_matrix,1,FUN=function(x) sum(!is.na(x)))))<4) {
#dlgMessage("Not all rows has 4 or more non-NA values,group label is enforced.")
#impgrp=character()
#while (length(impgrp)==0) {
#impgrp=dlgList(names(group),title="Which group label?")$res}
#impgrpo=as.factor(group[[which(names(group)==impgrp)]])
#rtransifile=msImpute(quantile_normalized_matrix,method="v2-mnar",group=impgrpo)
#} else {

#if no grouping method chosen, use a different imputation method from the same package 
rtransifile=msImpute(quantile_normalized_matrix,method="v2")}
rawdatafile[,7:ncol(rawdatafile)]=2^rtransifile-1
}}
if (impnorm=="Only normalize") {
#median normalize dataset 
rtransifile=normalize_data_dm(rintfile,normalize_func = 'medianCentering',log_base = 2, offset = 1)
rawdatafile[,7:ncol(rawdatafile)]=2^rtransifile
} 
if (impnorm=="Neither") { 
#if missing value, or 0, exists, add 1 to all values 
if (sum(rawdatafile[,7:ncol(rawdatafile)]==0)>0|sum(is.na(rawdatafile[,7:ncol(rawdatafile)]))>0) {
rawdatafile[,7:ncol(rawdatafile)][is.na(rawdatafile[,7:ncol(rawdatafile)])]=0
rawdatafile[,7:ncol(rawdatafile)]=rawdatafile[,7:ncol(rawdatafile)]+1} 
 }
}
################################################################################




 
