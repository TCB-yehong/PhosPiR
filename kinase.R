#######################download all libraries###################################
if(!require(KinSwingR)){BiocManager::install("KinSwingR",update=F,ask=F)}
library(KinSwingR)
if(!require(BiocParallel)){BiocManager::install("BiocParallel",update=F,ask=F)}
library(BiocParallel)
################################################################################

kindbinfo=data.frame(organism=c("Homo sapiens","Mus musculus","Rattus norvegicus"),
dborganism=c("human","mouse","rat"))
kinphylock=match(phylo,kindbinfo$organism)
if (is.na(kinphylock)) {
msgBox("Organism not supported for kinase analysis. Moving on!")
} else {
dir.create("Kinase Analysis",showWarnings=F)
msgBox("Now proceeding to kinase analysis. Kinase activity changes will be predicted for all comparisons with fold change 
and p-value statistics. 

Please choose a preferred p-value statistic, it will be used whenever available.")
pvalchoice=dlgList(c("T.test_Pvalue","T.test_FDR","Wilcox.test_Pvalue","Wilcox.test_FDR",
"ROTS_Pvalue","ROTS_FDR","RankProd_Pvalue","RankProd_PFP"),title="Preferred p-value statistic")$res

kinasedb=read.csv(paste0(codepth,"PhosphoSitePlus_Kinase_Dataset.csv"),stringsAsFactors=F)
kinorg=kindbinfo[kinphylock,2]
kinaseinfo=kinasedb[which(kinasedb$KIN_ORGANISM==kinorg&kinasedb$SUB_ORGANISM==kinorg),]
kinasefile=as.matrix(data.frame(kinase=kinaseinfo$GENE,substrate=kinaseinfo$SITE_...7_AA))
pwms=buildPWM(kinasefile)


for (i in 1:length(grpcompare)) {
if (paste0("Comparison",i,"_FoldChange") %in% colnames(statsRes)&
paste0("Comparison",i,"_",pvalchoice) %in% colnames(statsRes)) {
kinstat=data.frame(
annotation=paste(statsRes$Protein,statsRes$Gene.names,statsRes$Position,seq.window2,sep="|"),
peptide=seq.window2,
fc=statsRes[,match(paste0("Comparison",i,"_FoldChange"),colnames(statsRes))],
pval=statsRes[,match(paste0("Comparison",i,"_",pvalchoice),colnames(statsRes))])
} else if (sum(grepl(paste0("Comparison",i),colnames(statsRes))&grepl("Pvalue|FDR|PFP",colnames(statsRes)))>0&
paste0("Comparison",i,"_FoldChange") %in% colnames(statsRes)) {
kinstat=data.frame(
annotation=paste(statsRes$Protein,statsRes$Gene.names,statsRes$Position,seq.window2,sep="|"),
peptide=seq.window2,
fc=statsRes[,match(paste0("Comparison",i,"_FoldChange"),colnames(statsRes))],
pval=statsRes[,which(grepl(paste0("Comparison",i),colnames(statsRes))&grepl("Pvalue|FDR|PFP",colnames(statsRes)))[1]])
} else {kinstat=NULL}

if (!is.null(kinstat)) {
register(SnowParam(workers = detectCores()))
set.seed(123)
if (nrow(kinstat)<1000) {
npermut=nrow(kinstat)} else {npermut=1000}
scores=scoreSequences(input_data = kinstat, pwm_in = pwms,n=npermut,force_trim=TRUE)
swing_out <- swing(input_data = kinstat, pwm_in = pwms, pwm_scores = scores,return_network=T)
write.csv(swing_out$scores,paste0("Kinase Analysis/Comparison",i,"_swingScore.csv"),row.names=F)
if (nrow(swing_out$network)>0) {
#write.csv(swing_out$network,paste0("Kinase Analysis/Comparison",i,"_kinaseNetwork.csv"),row.names=F)
sigkinase=unique(swing_out$scores$kinase[which((swing_out$scores$p_greater<0.05|swing_out$scores$p_less<0.05)&swing_out$scores$swing_raw!=0)])
if (length(sigkinase)>0) {
for (j in sigkinase) {
tempkinsub=kinasefile[kinasefile[,1]==j,]
tempkinpwm=buildPWM(tempkinsub)
tiff(paste0("Kinase Analysis/Comparison",i,"_significantly_changing_kinase_",j,".tiff"), units="in", width=9, height=5, res=150)
tempkinfig=viewPWM(tempkinpwm,which_pwm=j,view_pwm=T,color_scheme="shapely")
dev.off()}
}

pwm_pval=scores$peptide_p
input_data=kinstat
p_cut_pwm = 0.05
p_cut_fc = 0.05
pwm_pval[, 3:ncol(pwm_pval)] <- ifelse(pwm_pval[, 3:ncol(pwm_pval)] > p_cut_pwm, 0, 1)
input_data$annotation <- paste(input_data$annotation, input_data$peptide, sep = "::")
pwm_pval$annotation <- paste(pwm_pval$annotation, pwm_pval$peptide, sep = "::")
data_merge <- unique(merge(input_data, pwm_pval, by = "annotation"))

datakinnet=c()
for (j in 1:nrow(data_merge)) {
for (k in 6:ncol(data_merge)) {
if (!is.na(data_merge[j,k])&data_merge[j,k]==1) {
datakinnet=rbind(datakinnet,cbind(colnames(data_merge)[k],data_merge[j,c("annotation","fc","pval")]))
}
}}
datakinnet=as.data.frame(datakinnet)
colnames(datakinnet)=c("Kinase","Substrate","FoldChange","Pvalue")
sigkinnet=datakinnet[datakinnet$Pvalue<0.05,]
write.csv(datakinnet,paste0("Kinase Analysis/Comparison",i,"_kinaseNetwork.csv"),row.names=F)
if (nrow(sigkinnet)>0) {
write.csv(sigkinnet,paste0("Kinase Analysis/Comparison",i,"_significant_kinaseNetwork.csv"),row.names=F)}
}
}}
}


