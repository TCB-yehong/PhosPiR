#######################download all libraries###################################
if(!require(KinSwingR)){BiocManager::install("KinSwingR",update=F,ask=F)}
library(KinSwingR)
if(!require(BiocParallel)){BiocManager::install("BiocParallel",update=F,ask=F)}
library(BiocParallel)
################################################################################
#kinase analysis is available for only human, mouse, and rat. If target organism is not 1 of the 3, skip this step
kindbinfo=data.frame(organism=c("Homo sapiens","Mus musculus","Rattus norvegicus"),
dborganism=c("human","mouse","rat"))
kinphylock=match(phylo,kindbinfo$organism)
if (is.na(kinphylock)) {
msgBox("Organism not supported for kinase analysis. Moving on!")
} else {
dir.create("Kinase Analysis",showWarnings=F)
msgBox("Now proceeding to kinase analysis. Kinase activity changes will be predicted for all comparisons with fold change 
and p-value statistics. 

#choose which p-value to use for calculating directional changes for kinase 
Please choose a preferred p-value statistic, it will be used whenever available.")
pvalchoice=dlgList(c("T.test_Pvalue","T.test_FDR","Wilcox.test_Pvalue","Wilcox.test_FDR",
"ROTS_Pvalue","ROTS_FDR","RankProd_Pvalue","RankProd_PFP"),title="Preferred p-value statistic")$res

#extract kinase info for target organism from library originated from PhosphoSitePlus 
kinasedb=read.csv(paste0(codepth,"PhosphoSitePlus_Kinase_Dataset.csv"),stringsAsFactors=F)
kinorg=kindbinfo[kinphylock,2]
kinorgchoice=dlgList(c("All organism",kinorg),title="Preferred kinase organism range")$res
if (kinorgchoice=="All organism") {
kinaseinfo=kinasedb[which(kinasedb$SUB_ORGANISM==kinorg),]
} else {
kinaseinfo=kinasedb[which(kinasedb$KIN_ORGANISM==kinorg&kinasedb$SUB_ORGANISM==kinorg),]}
kinasefile=as.matrix(data.frame(kinase=kinaseinfo$GENE,substrate=kinaseinfo$SITE_...7_AA))
winlen=min(nchar(kinasefile[,2]))
pwms=buildPWM(kinasefile,substrate_length=winlen)

for (i in 1:length(grpcompare)) {
missinfo=c()
#for each comparison, if stats result of the comparison include user selected p-val and a fold change,
#put ID annotation, sequence window, fold change and selected pval for the dataset in a dataframe 
if (paste0("Comparison",i,"_FoldChange") %in% colnames(statsRes)&
paste0("Comparison",i,"_",pvalchoice) %in% colnames(statsRes)) {
kinstat=data.frame(
annotation=paste(statsRes$Protein,statsRes$Gene.names,statsRes$Position,seq.window2,sep="|"),
peptide=seq.window2,
fc=statsRes[,match(paste0("Comparison",i,"_FoldChange"),colnames(statsRes))],
pval=statsRes[,match(paste0("Comparison",i,"_",pvalchoice),colnames(statsRes))])
#if user selected p-value is not part of the stats result for this comparison, check if any p-vals are calculated
#if yes, get the first p-val for p-val stats in the dataframe 
} else if (sum(grepl(paste0("Comparison",i),colnames(statsRes))&grepl("Pvalue|FDR|PFP",colnames(statsRes)))>0&
paste0("Comparison",i,"_FoldChange") %in% colnames(statsRes)) {
kinstat=data.frame(
annotation=paste(statsRes$Protein,statsRes$Gene.names,statsRes$Position,seq.window2,sep="|"),
peptide=seq.window2,
fc=statsRes[,match(paste0("Comparison",i,"_FoldChange"),colnames(statsRes))],
pval=statsRes[,which(grepl(paste0("Comparison",i),colnames(statsRes))&grepl("Pvalue|FDR|PFP",colnames(statsRes)))[1]])
#if the comparison doesn't include enough stats parameters, return empty dataframe variable
} else {kinstat=NULL
missinfo=rbind(missinfo,paste("Comparison ",i," doesnt have enough info: missing FC or pval/FDR/PFP."))
}

if (!is.null(kinstat)) {
register(SnowParam(workers = detectCores()))
set.seed(123)
#set number of permutation performed 
if (nrow(kinstat)<1000) {
npermut=nrow(kinstat)} else {npermut=1000}
#use package function to calculate swing score for each kinase 
scores=scoreSequences(input_data = kinstat, pwm_in = pwms,n=npermut,force_trim=TRUE)
swing_out <- swing(input_data = kinstat, pwm_in = pwms, pwm_scores = scores,return_network=T)
#record result in csv file 
write.csv(swing_out$scores,paste0("Kinase Analysis/Comparison",i,"_swingScore.csv"),row.names=F)
if (nrow(swing_out$network)>0) {
#write.csv(swing_out$network,paste0("Kinase Analysis/Comparison",i,"_kinaseNetwork.csv"),row.names=F)

#extract significant kinase from output 
sigkinase=unique(swing_out$scores$kinase[which((swing_out$scores$p_greater<0.05|swing_out$scores$p_less<0.05)&swing_out$scores$swing_raw!=0)])
#for each sig kinase, make a motif plot using data from PhosphoSitePlus kinase library 
if (length(sigkinase)>0) {
for (j in sigkinase) {
tempkinsub=kinasefile[kinasefile[,1]==j,]
if (nrow(tempkinsub)>0) {
tempkinpwm=buildPWM(tempkinsub)
tiff(paste0("Kinase Analysis/Comparison",i,"_significantly_changing_kinase_",j,".tiff"), units="in", width=9, height=5, res=150)
tempkinfig=viewPWM(tempkinpwm,which_pwm=j,view_pwm=T,color_scheme="shapely")
dev.off()}}
}

#get score table for substrate in each predicted kinase 
pwm_pval=scores$peptide_p
#get input file 
input_data=kinstat
p_cut_pwm = 0.05
p_cut_fc = 0.05
#change score table to binary, where sig (by above cutoff) changes for a kinase-substrate pair is given 1, rest are 0
pwm_pval[, 3:ncol(pwm_pval)] <- ifelse(pwm_pval[, 3:ncol(pwm_pval)] > p_cut_pwm, 0, 1)
#setup annotation for input entries and substrate entries (from library)
input_data$annotation <- paste(input_data$annotation, input_data$peptide, sep = "::")
pwm_pval$annotation <- paste(pwm_pval$annotation, pwm_pval$peptide, sep = "::")
#get input data that's also a substrate for any predicted kinase, extract their scores in all predicted kinase 
data_merge <- unique(merge(input_data, pwm_pval, by = "annotation"))

#for every kinase-substrate pair that's significantly changing, get comparison stats result and put them together
datakinnet=c()
for (j in 1:nrow(data_merge)) {
for (k in 6:ncol(data_merge)) {
if (!is.na(data_merge[j,k])&data_merge[j,k]==1) {
datakinnet=rbind(datakinnet,cbind(colnames(data_merge)[k],data_merge[j,c("annotation","fc","pval")]))
}
}}

#polish dataframe format
datakinnet=as.data.frame(datakinnet)
if (nrow(datakinnet)>0) {
colnames(datakinnet)=c("Kinase","Substrate","FoldChange","Pvalue")
#extract a version with only significant comparison changes 
sigkinnet=datakinnet[datakinnet$Pvalue<0.05,]
#record both lists in csv file if they are not empty 
write.csv(datakinnet,paste0("Kinase Analysis/Comparison",i,"_kinaseNetwork.csv"),row.names=F)
if (nrow(sigkinnet)>0) {
write.csv(sigkinnet,paste0("Kinase Analysis/Comparison",i,"_significant_kinaseNetwork.csv"),row.names=F)}
}
}
}} 
missinfo=as.data.frame(missinfo)
write.csv(missinfo,"failed_due_to_information_missing.csv")
}


