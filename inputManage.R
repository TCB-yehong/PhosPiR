#######################download all libraries###################################
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#	BiocManager::install(version = BiocManager::version(),ask=F,update=F)
if(!require(proBatch)){BiocManager::install("proBatch",update=F,ask=F)}
library(proBatch)
if(!require(msImpute)){	BiocManager::install("msImpute",update=F,ask=F)}
library(msImpute)
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
################################################################################

#######################deposit location#########################################
setwd(dlg_dir(default = getwd(),title="Where should results be stored?")$res)
fdn=paste0("Phosphoproteomics_analysis_",Sys.Date())
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
} else {
fname=dlg_open(title = "Select input file with correct format",multiple=F)$res
}
################################################################################

#######################maxquant input formatting################################
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
rawdatafile=rawdatafile[which(rawdatafile$Reverse==""&rawdatafile$Potential.contaminant==""),-c(7,8)]
rawinfofile=rawinfofile[which(rawdatafile$Reverse==""&rawdatafile$Potential.contaminant==""),]
} else {
rawdatafile=read.csv(fname,stringsAsFactors=F)
colnames(rawdatafile)[1:6]=c("Protein","Protein.names","Gene.names","Amino.acid","Position","Sequence.window")
msg_box("Please make sure the file follows the described
format, and have reverse sequence and potential 
contaminants removed to ensure accurate analysis. 
Thank you!")
}
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
if (rcoff=="Remove rows with missing value >50%") {rcoffv=0.5} 
else if (rcoff=="Remove rows with missing value >30%") {rcoffv=0.3}
else {rcoffv=0.1}
rpass=apply(rintfile,1,FUN=function(x) (sum(x==0)/length(x))<rcoffv)
rintfile=rintfile[which(rpass),]
rawdatafile=rawdatafile[which(rpass),]
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
################################################################################

#######################universal format#########################################
#has more to do with non-maxquant input when it's not in the right format
#at the moment this is user resposibility
#will implement in newer versions some build-in data formatting to help user 
#transform their data to a usable format for the pipeline
################################################################################

#######################user setup options#######################################
#setup organism
phylo=dlgList(c("Homo sapiens","Mus musculus","Rattus norvegicus",
"Saccharomyces cerevisiae","Caenorhabditis elegans","Danio rerio",
"Drosophila melanogaster","Escherichia coli","Arabidopsis thaliana",
"Other"),preselect="Homo sapien",title="Select data organism source")$res
if(phylo=="Other") {phylo=dlgInput("Please enter organism scientific name: (Capitalize first 
letter and nothing else, space between words)")$res}

#setup imputation and normalization
impnorm=dlgList(c("Neither","Only normalize","Both normalize and impute"),
title="Imputation and normalization")$res

#setup which step to run
pipestep=dlgList(c("Annotation","Overview Figures","Differential Analysis","All steps"),
title="Select which analysis to run")$res
################################################################################

#######################group setup##############################################
source(paste0(codepth,"groupSetup.R"))
group=grouping()
grpcompare=grpcpare(group)
################################################################################

#######################normalization and imputation#############################
msgBox("Now proceeding to normalization and imputation.")
if (impnorm=="Both normalize and impute") {
#dlgMessage("Caution with imputation, if there's time, another run without 
#imputation is recommended as control")

rintfile[rintfile==0]=NA
quantile_normalized_matrix = normalize_data_dm(as.matrix(rintfile),normalize_func = 'quantile',log_base = 2, offset = 1)
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
impgrpo=as.factor(group[[which(names(group)==impgrp)]])
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
rtransifile=msImpute(quantile_normalized_matrix,method="v2")}
rawdatafile[,7:ncol(rawdatafile)]=2^rtransifile-1
}}
if (impnorm=="Only normalize") {
rtransifile=normalize_data_dm(rintfile,normalize_func = 'medianCentering',log_base = 2, offset = 1)
rawdatafile[,7:ncol(rawdatafile)]=2^rtransifile
} 
if (impnorm=="Neither") { rawdatafile[,7:ncol(rawdatafile)]=rawdatafile[,7:ncol(rawdatafile)]+1 }

################################################################################




 