#guide on/off switch for setting up group
guide=dlgMessage("Enable group setup guide? (warning: quite wordy)",type="yesno")$res
#recording sample name information
samplenames=substring(colnames(rintfile),
unlist(lapply(colnames(rintfile),FUN=function(x) gregexpr("\\.",x)[[1]][1]))+1,
nchar(colnames(rintfile)))
dir.create("Group Information")
write.table(samplenames,"Group Information/input_sample_names.csv",sep=",",col.names=F)

#group function, sets up sample groups
grouping=function() {
cat("Grouping input guidelines:
In your phosphoproteomics data file, samples are represented by intensity columns.
Each intensity column matches one sample in MS run. Some samples are repeats, some
could be time series sample grouped by experiment length. some could be grouped in 
different ways based on different experiment parameters. Please define your grouping
of samples here. simply enter how many ways your samples can be grouped, then for every
way of grouping, enter a series of numbers separated by comma (no space) to represent 
each sample in order. The total count of numbers should equal to sample number; and 
samples in the same group should be given the same number. Ex: 10 samples, 2 conditions 
with 5 repeats per condition could be grouped in 1 way, written as 1,1,1,1,1,2,2,2,2,2
")
if (guide=="yes") {
msg_box(message=c("Grouping input guidelines:

In your phosphoproteomics data file, samples are represented
by intensity columns. Each intensity column matches one 
sample in MS run. Some samples are repeats, some could be 
time series sample grouped by experiment length. some could 
be grouped in different ways based on different experiment 
parameters. 

Please define your grouping of samples here. simply enter 
how many ways your samples can be grouped, then for every 
way of grouping, enter a series of numbers separated by 
comma (no space) to represent each sample in order. The 
total count of numbers should equal to sample number; and
samples in the same group should be given the same number. 

Ex: 10 samples, 2 conditions with 5 repeats per condition 
could be grouped in 1 way, written as 1,1,1,1,1,2,2,2,2,2

IF you cannot remember your sample label and order from 
the data, it is extracted for you in file:", 
fdn,
"/Group Information/input_sample_names.csv"))}
grp=numeric()
while(length(grp)==0){
grp=as.integer(dlg_input("How many ways can the samples be grouped (type a number): ","1")$res)}
if(is.na(grp)) {grp=1}
grping=list()
for (i in 1:grp) {
grpby=dlgList(c("Genotype","Phenotype","Strain","Treatment 1","Treatment 2","Other"),title=paste("Group",i,"is grouped by:",collapse=" "))$res
if (grpby=="Other") {
grpby=as.character(dlgInput(paste("Group",i,"is grouped by:",collapse=" "))$res)}
grpinp=dlg_input(paste0("Please enter grouping order of group ",i, ": "))$res
while (length(unlist(strsplit(grpinp,",")))!=length(samplenames)) {
grpinp=dlg_input(paste0("Please double check sample number. Enter grouping order of group ",i, ": "))$res}
grping[[i]]=as.numeric(factor(unlist(strsplit(grpinp,","))))
names(grping)[[i]]=grpby
}
grpingdf=as.data.frame(grping)
rownames(grpingdf)=samplenames
#colnames(grpingdf)=paste0("grouping",1:length(grping))
write.csv(grpingdf,"Group Information/sample_grouping_order.csv")
return(grping)}

#comparison function, sets up comparisons within groups
grpcpare=function(group) {
if (guide=="yes") {
msg_box("Group comparison info input guidelines:

Comparison of different conditions will yield proteins or sites 
of interest (i.e. significantly differentially expressed proteins 
or phosphosites). 

Please indicate the intended comparisons (at the moment only
takes 1 grouping into consideration at a time) by stating 
which grouping would you want to do comparisons, and 
which conditions should be compared to each other.")

msg_box("Input format is similar to previous input format. Ex. 12
samples grouped in two ways: 1,1,1,1,2,2,2,2,3,3,3,3 and
1,1,1,1,1,1,2,2,2,2,2,2. In the frst 2 comparisons, 
condition 2 and condition 3 are compared to condition 1 
separately from grouping order 1, in the 3rd comparison, 
condition 1 is compared to condition 2 from grouping
order 2.

To set up the comparison, choose grouping 1, then write
'1,2' for the first comparison, then choose 'OK' to make
another comparison. Choose grouping 1 then write '1,3'
for the second comparison, then choose 'OK'. Choose
grouping 2 for the third comparison, then write
'1,2'. Choose 'Cancel' to finish comparison set up.

Note: If data is paired, i.e. condition 1 and condition 2
in the comparison used the same mice, please make sure the
pairs are in the same order.

")}

for (i in names(group)) {
cat(i,"group groups the samples as follwing:",paste(group[[i]],collapse=","),"\n")}
grpcpr=list()
grpcprn=c()
idx=1
while (idx>0) {
grpnum=dlgList(names(group),title="Which group to compare?")$res
cnprnum=dlg_input(paste("Which 2 (or more)", grpnum, "types would you like to compare: ",collapse=" "))$res
while (any(!(as.integer(unlist(strsplit(cnprnum,","))) %in% unique(group[[grpnum]])),length(unique(unlist(strsplit(cnprnum,","))))==1)) {
cnprnum=dlg_input(paste("Please check input. Which 2 (or more)", grpnum, "types would you like to compare: ",collapse=" "))$res}
grptemp=data.frame(intcol=which(group[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,",")))),group=group[[grpnum]][which(group[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,","))))])
colnames(grptemp)[2]=grpnum
rownames(grptemp)=samplenames[which(group[[grpnum]] %in% as.integer(unlist(strsplit(cnprnum,","))))]
grpcpr[[idx]]=grptemp[order(grptemp[,2]),]
grpcprn=c(grpcprn,grpnum)
lp=ok_cancel_box("Would you like to make another comparison?")
if (lp) {idx=idx+1
} else {cat(length(grpcpr),"comparisons have been set up.")
idx=0}}
names(grpcpr)=grpcprn
cprinfo=data.frame(Comparisons=paste0("Comparison",1:length(grpcpr)))
cprinfo$Groups=grpcprn
cprinfo$Types=unlist(lapply(grpcpr, FUN=function(x) paste(sort(unique(x[,2]),decreasing=T),collapse=" vs ")))
cprinfo$Samples_Involved=unlist(lapply(grpcpr, FUN=function(x) paste(rownames(x),collapse=";")))
write.csv(cprinfo,"Group Information/Comparison_Information.csv",row.names=F)
return(grpcpr)
}
