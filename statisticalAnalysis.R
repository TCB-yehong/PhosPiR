#######################download all libraries###################################
if(!require(ggrepel)){install.packages("ggrepel")}
library(ggrepel)
if(!require(R.utils)){install.packages("R.utils")}
library(R.utils)
if(!require(gtools)){install.packages("gtools")}
library(gtools)
if(!require(gridExtra)){install.packages("gridExtra")}
library(gridExtra)
################################################################################

#######################test structure setup#####################################
#structure function for statistical testing and making volcano plots 
statFunc=function(transf,group,grpcompare) {
cpnum=length(grpcompare)
cat(cpnum,"comparisons will be analyzed.")

#make a matrix with tests as column and comparisons as row 
testsl=matrix(0,ncol=7,nrow=cpnum)
colnames(testsl)=c("FC","SEM","T-Test","Wilcoxon-Test","ROTS","RankProd","ANOVA")
rownames(testsl)=paste0("comparison",c(1:cpnum))
pairv=c()
#setup each comparison with user options
for (i in 1:cpnum) {
dlg_message(paste0("Comparison ",i," compares group ",paste(unique(grpcompare[[i]][,2]),collapse=", ")," from ",paste0(names(grpcompare)[i]," grouping."),"
"),type="ok")
grp2l=c()
#ask for pair information 
pairl=dlg_message(paste0("Is data paired in comparison ",i,"?"), "yesno")$res
pairv=c(pairv,pairl)
#if comparison includes only 2 groups, FC and SEM auto selected, then user select any or all of the stats tests
if (length(unique(grpcompare[[i]][,2]))==2) {
grp2l=c(grp2l,1)
cat("Fold change and SEM will be calculated.")
wtests=dlg_list(c("1.T-Test","2.Mann-Whitney-Wilcoxon Test","3.Reproducibility-Optimized Test Statistic (ROTS)",
"4.Rank Product (always assume unpaired)"),multiple=T,title=paste0("Which statistical tests for comparison ",i,"?"))$res
#record the number label (+2) of the selected tests 
testyes=unlist(lapply(wtests,FUN=function(x) as.numeric(unlist(strsplit(x,"\\."))[1])))+2
#mark 1 in the earlier matrix at the matching comparison and tests' row and column in the matrix 
testsl[i,c(1,2,testyes)]=rep(1,length(testyes)+2)} else {

#if comparison includes more than 2 groups, ANOVA test auto selected 
grp2l=c(grp2l,0)
cat("ANOVA will be performed.")
dlgMessage("Please make sure 'normalization and imputation' option was selected for valid ANOVA result")
testsl[i,7]=1
}}
#pairing info change to numerical binary
pairv[pairv=="yes"]=1
pairv[pairv=="no"]=2
pairv=as.numeric(pairv)
grp2l=as.logical(grp2l)


cat("Thank you! Analysis in progress, please wait...")
#keep only numeric columns 
data=transf[,-(1:6)]
transfana=transf
transfanan=c()
for (i in 1:nrow(testsl)) {
#for each comparison, check if pairing info is reasonably correct, if not change info to none pairing 
if (pairv[i]==1) {
pairlenl=sum(grpcompare[[i]][,2]==unique(grpcompare[[i]][,2])[2])==sum(grpcompare[[i]][,2]==unique(grpcompare[[i]][,2])[1])
if (!pairlenl) {cat("In comparison",i,"not all samples are paired, proceeding to unpaired calculations...")
pairv[i]=2}
}
#check for each test, if they are selected, if yes (1), perform test with corresponding functions 
if (testsl[i,1]) {
resc=fc(pairv[i],grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
if (testsl[i,2]) {
resc=sem(pairv[i],grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
if (testsl[i,3]) {
resc=tt(pairv[i],grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
if (testsl[i,4]) {
resc=wc(pairv[i],grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
if (testsl[i,5]) {
resc=rotz(pairv[i],grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
if (testsl[i,6]) {
resc=rp(grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
if (testsl[i,7]) {
resc=anva(pairv[i],grpcompare[[i]],data,i)
transfana=cbind(transfana,resc)
descb=colnames(resc)
transfanan=c(transfanan,descb)}
}

#save result of analysis 
transfana=as.data.frame(transfana)
dir.create("Statistical Analysis")
write.csv(transfana,"Statistical Analysis/data_and_analysis_results.csv",row.names=F)
if ("rawinfofile" %in% ls()) {
write.csv(rawinfofile,"Statistical Analysis/Data_background_info.csv",row.names=F)}

#prompt sig entries cutoff selection for stats result 
cat("Analysis complete. Please choose cutoffs for significantly changing entries
(your 'Significant Lists').")

testsep=unlist(lapply(transfanan,FUN=function(x) unlist(strsplit(x,"_"))[1]))
dir.create("Statistical Analysis/Significant Lists",showWarnings=F)
coinfo=c()
for (i in 1:length(grpcompare)) {
#get a vector of stats result for each comparison 
testrow=grep(as.character(i),testsep)

#prompt the list of stats for user to choose for setting up cutoff 
cotimes=0
coselec=1
while(coselec==1) {
cochoice=c()
while(length(cochoice)==0) {
cochoice=dlgList(transfanan[testrow],multiple=T,title=paste0("Select Comparison ",i," statistics cutoff"))$res}
cocol=match(cochoice,colnames(transfana))
#prompt for cutoff numerical value after stats options are chosen 
conum=c()
srcandi=c()
for (j in 1:length(cochoice)) {
conumeach=c()
while(length(conumeach)==0) {
conumeach=as.numeric(dlgInput(paste0("Please input an unlogged numerical cutoff value for ",cochoice[j]))$res)}
conum=c(conum,conumeach)
#set fold change cutoff to be smaller than selected value, and set pval cutoff to be larger than selected value 
if (grepl("FoldChange",cochoice[j])) {
srcandi=c(srcandi,which(abs(transfana[,cocol[j]])>conumeach))} else {
srcandi=c(srcandi,which(transfana[,cocol[j]]<conumeach))}
}
#display number of sig entries based on the earlier cutoff 
keepsigl=dlgMessage(paste0("Based on the selected cutoff, there are ",sum(table(srcandi)==length(cochoice)),
" significant entries for comparison ",i,". Keep significant list for further analysis?"),type="yesno")$res
#if user decided to keep this list, generate list info and append to a new row. 
#And extract data for this list with cutoff stats' columns and save to csv file 
if(keepsigl=="yes") {
cotimes=cotimes+1
if(length(cochoice)==1){
coinfo=rbind(coinfo,c(paste0("Comparison",i,"cutoff",cotimes),paste0(rownames(grpcompare[[i]]),collapse=","),paste0(grpcompare[[i]][,1]+6,collapse=","),names(grpcompare)[i],paste0(sort(unique(grpcompare[[i]][,2]),decreasing=T),collapse=" vs "),cochoice,conum,cocol))
} else {
coinfo=rbind(coinfo,c(paste0("Comparison",i,"cutoff",cotimes),paste0(rownames(grpcompare[[i]]),collapse=","),paste0(grpcompare[[i]][,1]+6,collapse=","),names(grpcompare)[i],paste0(sort(unique(grpcompare[[i]][,2]),decreasing=T),collapse=" vs "),paste0(cochoice,collapse=","),paste0(conum,sep="",collapse=","),paste0(cocol,sep="",collapse=",")))}
sigrec=transfana[as.numeric(names(table(srcandi))[which(table(srcandi)==length(cochoice))]),c(1:6,cocol)] #sig entries cbind with name rows in transf
write.csv(sigrec,paste0("Statistical Analysis/Significant Lists/Comparison",i,"_cutoff",cotimes,"_sigEntries.csv"),row.names=F)
}

#repeat if user would like to select another cutoff 
comore=dlgMessage(paste0("Would you like to select another set of cutoff for comparison ",i,"?"),type="yesno")$res
if (comore=="no"){
if (i==length(grpcompare)) {coselec=0
msgBox("Done! Thank you!")} else {
coselec=0
msgBox(paste0("Moving on to comparison ",i+1,"."))}}
}
}

#save all sig lists' info to csv file 
coinfo=as.data.frame(coinfo,stringsAsFactors=F)
colnames(coinfo)=c("Siglist","Samples_involved","Data_columns_involved","Group_category","Group_labels_involved","Cutoff_stats","Cutoff_values","Cutoff_value_columns")
write.csv(coinfo,"Statistical Analysis/Siglists_cutoff_info.csv",row.names=F)

#volcano plot
#choose max. 4 comparisons from 2 group comparisons (check if it's less than 4)
#choose 1 stats test result to plot (from stats that all comparisons have)
#plot volcano plot + save in tiff and eps
#ask if another one should be made

#ask if user would like to make volcano plots 
volcselec=dlgMessage("Make volcano plots?",type="yesno")$res
if (volcselec=="yes") {
volcyes=T} else {volcyes=F}
volccount=1
volcrecord=c()

#if user select yes, get comparison number and grouping method and groups involved in the comparison
while (volcyes) {
volcop=c()
grpcpidx=c()
for (i in 1:length(grpcompare)) {
#if (length(unique(grpcompare[[i]][,2]))>=2) {
volc1op=paste(
"Group",paste0("'",names(grpcompare)[i],"',"),
paste(sort(unique(grpcompare[[i]][,2]),decreasing=T),collapse=" vs "),collapse=" ")
volcop=c(volcop,volc1op)
grpcpidx=c(grpcpidx,i)}
#}

#prompt comparison selection GUI, if user chose more than 4 comparisons, reprompt 
if (length(volcop)>0) {
msgBox("Please select a maxinum of 4 comparisons to include in the volcano plot")
volcgrpselec=dlgList(volcop,title="Select comparisons for volcano plot",multiple=T)$res
while (length(volcgrpselec)>4) {
volcgrpselec=dlgList(volcop,title="Please choose maximum 4 comparisons",multiple=T)$res}
cpridx=grpcpidx[match(volcgrpselec,volcop)]
#check if all includes at least 1 stats test, if not, omit comparison that doesn't, give warning,
#if all doesn't, give warning and move on
cprgdidx=c()
for (i in cpridx) {
tempgrpcol=grep(paste0("Comparison",as.character(i),"_"),transfanan)
if (length(grep("Pvalue",transfanan[tempgrpcol]))+length(grep("FDR",transfanan[tempgrpcol]))>0) {
cprgdidx=c(cprgdidx,i)}}
if (length(cprgdidx)==0) {
msgBox("Note: no comparisons have P-value or FDR calculated, this volcano plot will not be made.")
} else if (length(cpridx)>length(cprgdidx)) {
msgBox("Note: some comparisons are missing P-value and FDR calcuation, they are omitted from plot.")} 

#if more than 1 stats test was performed, extract all columns that fits the stat test descriptions
if (length(cprgdidx)>0) {
tempgdgrpcol=transfanan[grep(paste(unlist(lapply(cprgdidx,FUN=function(x) paste0("Comparison",x,"_"))),collapse="|"),
#unlist(lapply(transfanan,FUN=function(x) strsplit(x,"_")[[1]][1])))]
transfanan)]

#processing tukey post-hoc specifically as so many sets of fold changes are included 
if (sum(grepl("Tukey",tempgdgrpcol))>0) {
multicol=grep(unlist(strsplit(tempgdgrpcol[grep("Tukey",tempgdgrpcol)[1]],"_"))[1],tempgdgrpcol)
multicolidx=match(c("Gene.names",tempgdgrpcol[grep("Tukey",tempgdgrpcol)]),colnames(transfana))
multivolcname=rep(transfana[,multicolidx[1]],(length(multicolidx)-1)/2)
multivolcfc=c()
for (i in seq(2,length(multicolidx),by=2)) {
multivolcfc=c(multivolcfc,transfana[,multicolidx[i]])}
multivolcpval=c()
for (i in seq(3,length(multicolidx),by=2)) {
multivolcpval=c(multivolcpval,transfana[,multicolidx[i]])}
multivolcgroup=rep(unique(substring(tempgdgrpcol[grep("Tukey",tempgdgrpcol)],1,
unlist(lapply(tempgdgrpcol[grep("Tukey",tempgdgrpcol)],FUN=function(x) gregexpr("\\_",x)[[1]][2]))-1)),
each=nrow(transfana))
multicpr=unique(unlist(lapply(tempgdgrpcol[multicol],FUN=function(x) strsplit(x,"_")[[1]][1])))
tempgdgrpcol=tempgdgrpcol[-multicol]
multicpridx=numeric()
for (i in cprgdidx) {
if (length(grep(i,multicpr))==0) {
multicpridx=c(multicpridx,i)
}}
cprgdidx=cprgdidx[multicpridx]
} 

#if there is at least 1 stats column for a comparison, extract info and generate parameters needed for volcano plot 
if (length(tempgdgrpcol)>0) {
tempgdsubname=substring(tempgdgrpcol,
unlist(lapply(tempgdgrpcol,FUN=function(x) gregexpr("\\_",x)[[1]][1]))+1,
nchar(tempgdgrpcol))
tempstatcolname=tempgdsubname[grep("Pvalue|FDR|PFP",tempgdsubname)]
statchoice=names(table(tempstatcolname))[which(table(tempstatcolname)==length(cprgdidx))]
volcstatselec=dlgList(statchoice,title="Select stat to plot")$res
volcstatcols=tempgdgrpcol[grep(paste(c(volcstatselec,"FoldChange"),collapse="|"),tempgdgrpcol)]
volcinfocolidx=match(c("Gene.names",volcstatcols),colnames(transfana))
#use these info to get volc info
#name column, FC column, Pval column, group column
volcname=rep(transfana[,volcinfocolidx[1]],(length(volcinfocolidx)-1)/2)
volcfc=c()
for (i in seq(2,length(volcinfocolidx),by=2)) {
volcfc=c(volcfc,transfana[,volcinfocolidx[i]])}
volcpval=c()
for (i in seq(3,length(volcinfocolidx),by=2)) {
volcpval=c(volcpval,transfana[,volcinfocolidx[i]])}
volcgroup=rep(paste0("Comparison",cprgdidx),each=nrow(transfana)) } else {
volcname=c()
volcfc=c()
volcpval=c()
volcgroup=c()
volcstatselec="Tukey post hoc pvalue"
}
if ("multicol" %in% ls()) {
volcname=c(volcname,multivolcname)
volcfc=c(volcfc,multivolcfc)
volcpval=c(volcpval,multivolcpval)
volcgroup=c(volcgroup,multivolcgroup) 
rm(multicol)
}


#calculate modified fc and modified pval
volcnlogpval=-log10(volcpval)
volclogFC=rep(0,length(volcfc))
volclogFC[volcfc>0]=log2(volcfc[volcfc>0])
volclogFC[volcfc<0]=-log2(abs(volcfc[volcfc<0]))

#setup fold change and p-value cutoff variable based on user input 
msgBox("Please choose fold change and p-value cutoff")
volcstatcf=dlgList(c("2 and 0.05, respectively", "1.5 and 0.05, respectively", 
"2 and 0.01, respectively", "1.5 and 0.01, respectively"),
title="cutoff choices")$res
volcstatcfc=match(volcstatcf,c("2 and 0.05, respectively", "1.5 and 0.05, respectively", 
"2 and 0.01, respectively", "1.5 and 0.01, respectively"))
if (volcstatcfc==1) {
volcpvalcf=0.05
volcfccf=2
} else if (volcstatcfc==2) {
volcpvalcf=0.05
volcfccf=1.5
} else if (volcstatcfc==3) {
volcpvalcf=0.01
volcfccf=2
} else {
volcpvalcf=0.01
volcfccf=1.5
}

#check if both side sig count < 30 (include gene label if yes, else no labels)
volclabon=sum(volcfc>volcfccf&(volcpval<volcpvalcf))<30&(sum(volcfc<(-volcfccf)&(volcpval<volcpvalcf))<30)
#organize data needed for volcano plots
dataf=data.frame(name=volcname,pval=volcpval,nlogpval=volcnlogpval,fc=volcfc,logFC=volclogFC,lab=volcgroup)
#check if the selected number of comparisons is less than 4, if not, ask user to re-select 
if (length(unique(dataf$lab))>4) {
volcgrpselec2=dlgList(unique(dataf$lab),title="Select 1 to 4 comparisons",multiple=T)$res
while (length(volcgrpselec2)>4) {
volcgrpselec2=dlgList(unique(dataf$lab),title="Please choose maxsimum 4 comparisons",multiple=T)$res}
dataf=dataf[dataf$lab %in% volcgrpselec2,]
}

#plot vocano plot with volcano plot function and save as tiff 
tiff(paste0("Statistical Analysis/Volcano_plot-",volccount,".tiff"), units="in", width=14, height=8, res=150)
volcalpha=0.6
volcplot(volcalpha,dataf, volcfccf,volcpvalcf,volclabon)
dev.off()

#plot vocano plot with volcano plot function and save as EPS 
setEPS()
postscript(paste0("Statistical Analysis/Volcano_plot-",volccount,".eps"), width=14, height=8)
volcalpha=1
volcplot(volcalpha,dataf, volcfccf,volcpvalcf,volclabon)
dev.off()

#record user choices for each volcano plot 
tempvolcrec=c(paste0(unique(dataf$lab),collapse=","),volcstatselec,volcstatcf)
names(tempvolcrec)=c("Included comparisons","Pvalue statistic","Cutoffs")
volcrecord=rbind(volcrecord,tempvolcrec)
volccount=volccount+1}

volcanoth=dlgMessage("Make another volcano plot?",type="yesno")$res
if (volcanoth=="yes") {
volcyes=T} else {volcyes=F}
#} else {msgBox("Sorry no available 2-group comparisons, volcano plot cannot be made.")
#volcyes=F}
}}

#save user choices for each volcano plot to csv file 
rownames(volcrecord)=1:nrow(volcrecord)
write.csv(as.data.frame(volcrecord),"Statistical Analysis/Volcano_plot_info.csv")
return(transfana)
}
################################################################################

#######################statistical tests functions##############################
#function for calculating fold change 
fc=function(pair,compareinfo,data,i) {
if (pair==2) {
noneg=apply(data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[2])]],1,mean)/apply(data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[1])]],1,mean)
noneg[noneg<1]=-1/noneg[noneg<1]
wneg=as.data.frame(round(noneg, 2))
colnames(wneg)=paste0("Comparison",i,"_FoldChange")
return(wneg)} else {
noneg=data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[2])]]/data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[1])]]
noneg[noneg<1]=-1/noneg[noneg<1]
wneg=as.data.frame(round(apply(noneg,1,mean),2))
colnames(wneg)=paste0("Comparison",i,"_FoldChange")
return(wneg)}
}

#function for calculating standard error of mean for fold change 
sem=function(pair,compareinfo,data,i) {
if (pair==2) {
grp1m=apply(data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[1])]],1,mean)
grp1sd=apply(data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[1])]],1,sd)
grp2m=apply(data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[2])]],1,mean)
grp2sd=apply(data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[2])]],1,sd)
semnp=as.data.frame(((grp1sd/grp1m)^2+(grp2sd/grp2m)^2)^0.5)
colnames(semnp)=paste0("StandardErrorofMean_","Comparison",i)
return(semnp)
} else {
res=data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[2])]]/data[,compareinfo[,1][which(compareinfo[,2]==unique(compareinfo[,2])[1])]]
semp=data.frame(sem=apply(res,1,FUN=function(x) sd(x)/(length(x)^0.5)))
colnames(semp)=paste0("StandardErrorofMean_","Comparison",i)
return(semp)
}
}

#function for T-test 
tt=function(pair,compareinfo,data,i) {
grp1=compareinfo[2]==names(table(compareinfo[2]))[1]
grp2=compareinfo[2]==names(table(compareinfo[2]))[2]
if (pair==2) {
ttres=apply(log2(data+1),1,FUN=function(x) tryCatch(t.test(x[compareinfo[1][grp1]],x[compareinfo[1][grp2]])$p.value,error=function(err) 1))
ttresad=p.adjust(as.numeric(ttres),method="BH")
ttdf=data.frame(ttpval=ttres,ttfdr=ttresad)
colnames(ttdf)=c(paste0("Comparison",i,"_T-test_Pvalue"),paste0("Comparison",i,"_T-test_FDR"))
return(ttdf)
} else{
ttres=apply(log2(data+1),1,FUN=function(x) tryCatch(t.test(x[compareinfo[1][grp1]],x[compareinfo[1][grp2]],paired=T)$p.value,error=function(err) 1))
ttresad=p.adjust(as.numeric(ttres),method="BH")
ttdf=data.frame(ttpval=ttres,ttfdr=ttresad)
colnames(ttdf)=c(paste0("Comparison",i,"_T-test_Pvalue"),paste0("Comparison",i,"_T-test_FDR"))
return(ttdf)
}
}

#function for wilcox test 
wc=function(pair,compareinfo,data,i) {
grp1=compareinfo[2]==names(table(compareinfo[2]))[1]
grp2=compareinfo[2]==names(table(compareinfo[2]))[2]
if (pair==2) {
wcres=apply(log2(data),1,FUN=function(x) tryCatch(wilcox.test(x[compareinfo[1][grp1]],x[compareinfo[1][grp2]])$p.value,error=function(err) 1))
wcresad=p.adjust(as.numeric(wcres),method="BH")
wcdf=data.frame(wcpval=wcres,wcfdr=wcresad)
colnames(wcdf)=c(paste0("Comparison",i,"_Wilcox-test_Pvalue"),paste0("Comparison",i,"_Wilcox-test_FDR"))
return(wcdf)
} else{
wcres=apply(log2(data),1,FUN=function(x) tryCatch(wilcox.test(x[compareinfo[1][grp1]],x[compareinfo[1][grp2]],paired=T)$p.value,error=function(err) 1))
wcresad=p.adjust(as.numeric(wcres),method="BH")
wcdf=data.frame(wcpval=wcres,wcfdr=wcresad)
colnames(wcdf)=c(paste0("Comparison",i,"_Wilcox-test_Pvalue"),paste0("Comparison",i,"_Wilcox-test_FDR"))
return(wcdf)
}
}

#function for ROTS test 
rotz=function(pair,compareinfo,data,i) {
if (!require(ROTS)) {BiocManager::install("ROTS",update=F,ask=F)}
library(ROTS)

if (pair==2) {
rotsres=ROTS(log2(data[,compareinfo[,1]]),groups=compareinfo[,2],B=3000,K=nrow(data)/2)
rotsdf=data.frame(rotspval=rotsres$pvalue,rotsfdr=rotsres$FDR)
colnames(rotsdf)=c(paste0("Comparison",i,"_ROTS_Pvalue"),paste0("Comparison",i,"_ROTS_FDR"))
return(rotsdf)
} else{
rotsres=ROTS(log2(data[,compareinfo[,1]]),groups=compareinfo[,2],B=3000,K=nrow(data)/2,paired=T)
rotsdf=data.frame(rotspval=rotsres$pvalue,rotsfdr=rotsres$FDR)
colnames(rotsdf)=c(paste0("Comparison",i,"_ROTS_Pvalue"),paste0("Comparison",i,"_ROTS_FDR"))
return(rotsdf)
}
}

#function for rank product test
rp=function(compareinfo,data,i) {
if (!require(RankProd)) {BiocManager::install("RankProd",update=F,ask=F)}
library(RankProd)

rpres=RPadvance(data=data[,compareinfo[,1]], cl=as.numeric(as.vector(factor(compareinfo[,2],labels=c(0,1)))), origin=rep(1,nrow(compareinfo)),gene.names=rownames(data), num.perm=3000, logged=F, na.rm=T, plot=F, rand=123)
rpdf=data.frame(rppval=apply(rpres$pval,1,min),rpfdr=apply(rpres$pfp,1,min))
colnames(rpdf)=c(paste0("Comparison",i,"_RankProd_Pvalue"),paste0("Comparison",i,"_RankProd_PFP"))
return(rpdf)
}

#function that performs ANOVA and Tukey post hoc test for non-paired comparison, and LME model for paired comparison 
anva=function(pair,compareinfo,data,i) {
if (pair==2) {
anvares=apply(data[,compareinfo[,1]],1,FUN=function(x) aov(x~as.factor(compareinfo[,2])))
anvastatsall=c()
for (j in 1:length(anvares)) {
tempanvares=anvares[[j]]
tempanvapval=summary(tempanvares)[[1]][["Pr(>F)"]][1]
names(tempanvapval)="ANOVA_Pvalue"
tempanvaaj=p.adjust(tempanvapval,method="BH",n=length(anvares))
names(tempanvaaj)="ANOVA_FDR"
temptukey=TukeyHSD(tempanvares)
temptukeystats=c()
for (k in 1:nrow(temptukey[[1]])) {
temptukeystats=c(temptukeystats,temptukey[[1]][k,c(1,4)]) #following: changed diff to FC
afac=as.numeric(unlist(strsplit(rownames(temptukey[[1]])[k],"-")))
temptukeystats[length(temptukeystats)-1]=foldchange(mean(as.numeric(data[j,][compareinfo[,1][compareinfo[,2]==afac[1]]])),
mean(as.numeric(data[j,][compareinfo[,1][compareinfo[,2]==afac[2]]])))
}
names(temptukeystats)=paste0("Tukey-",rep(rownames(temptukey[[1]]),each=2),c("_FoldChange","_Pvalue"))
tempanvastatsall=c(tempanvapval,tempanvaaj,temptukeystats)
anvastatsall=rbind(anvastatsall,tempanvastatsall) }
anvadf=as.data.frame(anvastatsall)
colnames(anvadf)=paste0("Comparison",i,"_",colnames(anvadf))
rownames(anvadf)=NULL
return(anvadf)
} else{
if(!require(nlme)){install.packages("nlme")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
library(multcompView)
library(lsmeans)
library(nlme)
#dt=log2(data[,compareinfo[,1]]+1)
dt=data[,compareinfo[,1]]
condi=compareinfo[,2]
samp=rep(0,length(condi))
for (j in unique(compareinfo[,2])){
samp[condi==j]=seq(from=1,by=1,length.out=sum(condi==j))}
anvapval=c()
posthocpval=c()
for (j in 1:nrow(dt)){
dtdf=data.frame(sample=samp,condition=condi,value=as.numeric(dt[j,]))
dtdf$condition=factor(dtdf$condition)
model.cor = lme(value ~ condition, 
            random = ~1|sample,
            data=dtdf)
if (ACF(model.cor)$ACF[2]<=(-1)){
ACFcor=-0.999} 
else if (ACF(model.cor)$ACF[2]>=1) {
ACFcor=0.999} else {ACFcor=ACF(model.cor)$ACF[2]}
model = lme(value ~ condition, 
            random = ~1 |sample,
            correlation = corAR1(form = ~ 1 | sample,value=ACFcor),
            data=dtdf,
            method="REML")
anvapval=c(anvapval,anova.lme(model)$'p-value'[2])
leastsquare = lsmeans::lsmeans(model,
                      pairwise ~ condition,
                      adjust="tukey") 
posthocpval=rbind(posthocpval,summary(leastsquare$contrasts)$p.value)
}
anvapvalaj=p.adjust(as.numeric(anvapval),method="BH")
anvadf=as.data.frame(cbind(anvapval,anvapvalaj,posthocpval))
colnames(anvadf)=c(paste0("Comparison",i,"_ANOVA_Pvalue"),paste0("Comparison",i,"_ANOVA_FDR"),paste0("Comparison",i,"_Post-hoc-Tukey-",summary(leastsquare$contrasts)$contrast,"_Pvalue"))
return(anvadf)
}
}

#function to plot volcano plot with x and y density plots 
volcplot=function(volcalpha,dataf,volcfccf,volcpvalcf,volclabon) {

#expand top of the figure
expandy = function(plot, ymin=0) {

  max.y = max(layer_data(plot)$y, na.rm=TRUE)
  min.log = floor(log(max.y,base=50))

  expand_limits(y=c(ymin, ceiling(max.y/50^min.log)*50^min.log*1.1))
}

#actual volcano plot 
volcano <- ggplot(dataf,aes(logFC, nlogpval, color=lab, solid=F)) + 
  geom_point(alpha=volcalpha) + 
#  theme(legend.position=c(0.85,0.3), legend.justification=c(0,1)) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16)) +
  xlab("Log2 Fold Change") + ylab("-Log10 P-value") +
  labs(colour="Comparison Group") +
  geom_vline(xintercept = log2(volcfccf),linetype="dashed") + geom_vline(xintercept = -log2(volcfccf),linetype="dashed") +
  geom_hline(yintercept = -log10(volcpvalcf),linetype="dashed") 
 volcano <- volcano + expandy(volcano)

#add labels if sig entries are less than 30 from either (increasing and decreasing) sides 
if (volclabon) {
volcano <- volcano +
			geom_label_repel(data=subset(dataf, abs(fc) > volcfccf & pval < volcpvalcf),
            aes(logFC,nlogpval,label=name),size=3, fontface = "bold",
            arrow = arrow(length = unit(0.03,"npc"), type = "closed", ends = "last", angle = 15),
            force = 5 ,xlim  = c(min(dataf$logFC),max(dataf$logFC)),max.overlaps=100)}

# Marginal density plot of x (top panel)
xdensity <- ggplot(dataf, aes(logFC, fill=lab)) + 
  scale_fill_discrete(name="Comparison Group") +
  geom_density(alpha=volcalpha) + 
  theme(legend.position=c(0.85,1), legend.justification=c(0,1)) +
#  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16)) +
  xlab("Log2 Fold Change") + ylab("Density") 
xdensity
# Marginal density plot of y (right panel)
ydensity <- ggplot(dataf, aes(nlogpval, fill=lab)) + 
  geom_density(alpha=volcalpha) + 
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=16)) +
  xlab("-Log10 P-value") + ylab("Density") 
ydensity

#top right corner's empty plot 
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(), 
   panel.border = element_blank(),
   panel.background = element_blank(),
   axis.title.x = element_blank(),
   axis.title.y = element_blank(),
   axis.text.x = element_blank(), 
   axis.text.y = element_blank(),
   axis.ticks = element_blank()
     )

return(
grid.arrange(xdensity, blankPlot, volcano, ydensity, 
        ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4)))}
################################################################################