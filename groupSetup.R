#read sample names
samplenames=substring(colnames(rintfile),
unlist(lapply(colnames(rintfile),FUN=function(x) gregexpr("\\.",x)[[1]][1]))+1,
nchar(colnames(rintfile)))

dir.create("Group Information")
write.table(samplenames,"Group Information/input_sample_names.csv",sep=",",col.names=F)

#############################group setup########################################
#choose setup method 
gsetmethod=dlgList(c("Read from file (.csv or .xlsx)","Go through group setup process"),
title="Select group setup method")$res

if ("Go through group setup process" %in% gsetmethod) {
grp=numeric()
#prompt number of group repeatedly if no number is entered
while(length(grp)==0){
grp=as.integer(dlg_input("How many ways can the samples be grouped (type a number): ","1")$res)}
#if value is not integer, default grouping method count =1
if(is.na(grp)) {grp=1}

grping=list()
#setup each grouping method's orders
for (i in 1:grp) {
#choose or type a grouping method description
grpby=dlgList(c("Genotype","Phenotype","Strain","Treatment 1","Treatment 2","Other"),title=paste("Grouping method",i,"is grouped by:",collapse=" "))$res
if (grpby=="Other") {
grpby=as.character(dlgInput(paste("Grouping method",i,"is grouped by:",collapse=" "))$res)}

#setup group number
newgroupyes="yes"
#current group number
n=1
#setup place holder
grpinp=rep(0,length(samplenames))
samplenameinput=samplenames
while (newgroupyes=="yes") {
#select sample names that belong to the current group 
sampleSelect=dlgList(samplenameinput,multiple=T,title=paste0("Select samples in group ",n))$res
#replace place holder by current group number 
grpinp[which(samplenameinput %in% sampleSelect)]=n
#mark selected samples with group number 
samplenameinput[which(samplenameinput %in% sampleSelect)]=paste(samplenameinput[which(samplenameinput %in% sampleSelect)]," --Group",n)
#whether there is another group in this grouping method, continue repeat the same steps until no more new groups
newgroupyes=dlgMessage("Is there another group in this grouping method?", type="yesno")$res
n=n+1
}
#save order in a list
grping[[i]]=as.numeric(grpinp)
names(grping)[[i]]=grpby
}
group=grping

#change to dataframe format and save this information for re-run's order input
grpingdf=as.data.frame(grping)
rownames(grpingdf)=samplenames
#get original dataset name
ofname=strsplit(fname,"/|\\.")[[1]]
ofnames=ofname[length(ofname)-1]
write.csv(grpingdf,paste0("Group Information/",ofnames,"_sample_grouping_order_input.csv"))
}

if ("Read from file (.csv or .xlsx)" %in% gsetmethod) {
#prompt to select input file 
grpfile=dlg_open(title = "Select group setup input file with correct format (.csv or .xlsx)",multiple=F)$res
#read input based on whether it is a csv or xlsx file 
grpfileparts=strsplit(grpfile,"\\.")[[1]]
if (grpfileparts[length(grpfileparts)]=="csv") {
grpfiledf=read.csv(grpfile,stringsAsFactors=F,row.names=1)
} else if (grpfileparts[length(grpfileparts)]=="xlsx") {
grpfiledf=read.xlsx(grpfile,sheet=1,check.names=T,rowNames=T)
} else {
dlg_message("Only .csv and .xslx formats are supported!")
}

#check if file included all necessary info, if not, reselect file 
while (!all(samplenames %in% rownames(grpfiledf))|sum(is.na(grpfiledf))!=0) {
dlg_message("Some samples are missing group information! Please select a file with all samples' information")
grpfile=dlg_open(title = "Select group setup input file with correct format (.csv or .xlsx)",multiple=F)$res
grpfileparts=strsplit(grpfile,"\\.")[[1]]
if (grpfileparts[length(grpfileparts)]=="csv") {
grpfiledf=read.csv(grpfile,stringsAsFactors=F,row.names=1)
} else if (grpfileparts[length(grpfileparts)]=="xlsx") {
grpfiledf=read.xlsx(grpfile,sheet=1,check.names=T,rowNames=T)
} else {
dlg_message("Only .csv and .xslx formats are supported!")
}
}
#in case sample order is different from dataset, match sample order 
grpfiledfo=grpfiledf[match(samplenames,rownames(grpfiledf)),,drop=F]
#record group order in a list variable with grouping method as list item name
grping=list()
for (i in 1:ncol(grpfiledfo)) {
grping[[i]]=as.numeric(grpfiledfo[,i])
names(grping)[[i]]=colnames(grpfiledfo)[i]
}
group=grping

#save this info as grouping order information 
grpingdf=as.data.frame(grping)
rownames(grpingdf)=samplenames
write.csv(grpingdf,"Group Information/sample_grouping_order.csv")
}
################################################################################


#############################comparison setup###################################
#choose setup method
gcprmethod=dlgList(c("Read from file (.csv or .xlsx)","Go through comparison setup process"),
title="Select group comparison setup method")$res

if ("Go through comparison setup process" %in% gcprmethod) {
grpcpare=function(group) {
#print in R each group's order 
for (i in names(group)) {
cat(i,"group groups the samples as follwing:",paste(group[[i]],collapse=","),"\n")}

grpcpr=list()
grpcprn=c()
idx=1
cnprnumset=c()
#select grouping method the compared groups are from, then enter group numbers
while (idx>0) {
grpnum=dlgList(names(group),title="Which group to compare?")$res
cnprnum=dlg_input(paste("Which 2 (or more)", grpnum, "types would you like to compare: ",collapse=" "))$res
#prompt again if input is not in the correct format (input can't match to group info or only 1 group in input)
while (any(!(as.integer(unlist(strsplit(cnprnum,","))) %in% unique(group[[grpnum]])),length(unique(unlist(strsplit(cnprnum,","))))==1)) {
cnprnum=dlg_input(paste("Please check input. Which 2 (or more)", grpnum, "types would you like to compare: ",collapse=" "))$res}
#create dataframe to record column number of the selected groups from just intensity columns, and record group numbers
grptemp=data.frame(intcol=which(group[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,",")))),group=group[[grpnum]][which(group[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,","))))])
colnames(grptemp)[2]=grpnum
rownames(grptemp)=samplenames[which(group[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,","))))]
#record this info to list, and append list for new comparisons 
grpcpr[[idx]]=grptemp[order(grptemp[,2]),]
#record grouping method and input group numbers 
grpcprn=c(grpcprn,grpnum)
cnprnumset=c(cnprnumset,cnprnum)
#repeat if user wants to make more comparisons 
lp=dlgMessage("Would you like to make another comparison?", type="yesno")$res
if (lp=="yes") {idx=idx+1
} else {cat(length(grpcpr),"comparisons have been set up.")
idx=0}}

names(grpcpr)=grpcprn
#making a information table for compared samples, and record in csv file 
cprinfo=data.frame(Comparisons=paste0("Comparison",1:length(grpcpr)))
cprinfo$Groups=grpcprn
cprinfo$Types=unlist(lapply(grpcpr, FUN=function(x) paste(sort(unique(x[,2]),decreasing=T),collapse=" vs ")))
cprinfo$Samples_Involved=unlist(lapply(grpcpr, FUN=function(x) paste(rownames(x),collapse=";")))
write.csv(cprinfo,"Group Information/Comparison_Information.csv",row.names=F)

#make dataframe and save comparison information for re-run purpose
cprinfoset=data.frame(Grouping_method=grpcprn,Comparison=cnprnumset)
#get original dataset name
ofname=strsplit(fname,"/|\\.")[[1]]
ofnames=ofname[length(ofname)-1]
write.csv(cprinfoset,paste0("Group Information/",ofnames,"_comparison_setup_input.csv"),row.names=F)

return(grpcpr)
}
grpcompare=grpcpare(grping)
}

if ("Read from file (.csv or .xlsx)" %in% gcprmethod) {
#prompt to select input file 
cprfile=dlg_open(title = "Select comparison input file with correct format (.csv or .xlsx)",multiple=F)$res
#read input based on whether it is a csv or xlsx file 
cprfileparts=strsplit(cprfile,"\\.")[[1]]
if (cprfileparts[length(cprfileparts)]=="csv") {
cprfiledf=read.csv(cprfile,stringsAsFactors=F)
} else if (cprfileparts[length(cprfileparts)]=="xlsx") {
cprfiledf=read.xlsx(cprfile,sheet=1,check.names=T)
} else {
dlg_message("Only .csv and .xslx formats are supported!")
}

grpcpr=list()
grpcprn=c()
cnprnumset=c()
for (i in 1:nrow(cprfiledf)) {
#get grouping method and comparison groups from each row of input 
grpnum=cprfiledf[i,1]
if (gsetmethod!="Go through group setup process") {
grpnum=gsub(" ",".",grpnum)
}
cnprnum=cprfiledf[i,2]
#create dataframe to record column number of the selected groups from just intensity columns, and record group numbers
grptemp=data.frame(intcol=which(grping[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,",")))),grping=grping[[grpnum]][which(grping[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,","))))])
colnames(grptemp)[2]=grpnum
rownames(grptemp)=samplenames[which(grping[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,","))))]
#record this info to list, and append list for new comparisons 
grpcpr[[i]]=grptemp[order(grptemp[,2]),]
#record grouping method and input group numbers 
grpcprn=c(grpcprn,grpnum)
cnprnumset=c(cnprnumset,cnprnum)
}

names(grpcpr)=grpcprn
#making a information table for compared samples, and record in csv file 
cprinfo=data.frame(Comparisons=paste0("Comparison",1:length(grpcpr)))
cprinfo$Groups=grpcprn
cprinfo$Types=unlist(lapply(grpcpr, FUN=function(x) paste(sort(unique(x[,2]),decreasing=T),collapse=" vs ")))
cprinfo$Samples_Involved=unlist(lapply(grpcpr, FUN=function(x) paste(rownames(x),collapse=";")))
write.csv(cprinfo,"Group Information/Comparison_Information.csv",row.names=F)

grpcompare=grpcpr
}
################################################################################

