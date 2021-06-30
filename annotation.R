#######################download all libraries###################################
if(!require(biomaRt)){BiocManager::install("biomaRt",update=F,ask=F)}
library(biomaRt)
if(!require(Biostrings)){BiocManager::install("Biostrings",update=F,ask=F)}
library(Biostrings)
if(!require(GenomicAlignments)){BiocManager::install("GenomicAlignments",update=F,ask=F)}
library(GenomicAlignments)
if(!require(protr)){install.packages("protr")}
library(protr)
if(!require(UniprotR)){install.packages("UniprotR")}
library(UniprotR)
################################################################################

#######################getting reviewed accession###############################
ensembl=useMart("ensembl")
genedblist=listDatasets(ensembl)
orgsym=paste0(tolower(substr(phylo,1,1)),substr(phylo,gregexpr(" ",phylo)[[1]][1]+1,nchar(phylo)))
orgindb=grep(orgsym,genedblist$dataset)
if (length(orgindb) >0) {
dir.create("Annotation",showWarnings=F)
msgBox("Now proceeding to Annotation step (Can take very long to run if all steps are performed, can leave for overnight run).")
uniprotinfoyes=dlgMessage("Get Uniprot information (Can take very long to run)?",type="yesno")$res
if (phylo!="Homo sapiens") {
humaninfoyes=dlgMessage("Get human homolog information (Can take very long to run)?",type="yesno")$res}

getUniProt=function (id) 
{
    id <- as.character(id)
    n <- length(id)
    proteins <- vector("list", n)
    for (i in 1:n) {
        proteins[[i]] <- tryCatch(readFASTA(paste("https://www.uniprot.org/uniprot/", 
            id[i], ".fasta", sep = ""))[[1]], error=function(err) NA)
    }
    proteins
}

orgensdb=useDataset(genedblist$dataset[orgindb[1]],mart=ensembl)
#table(listAttributes(orgensdb)$page)
#listAttributes(orgensdb)[listAttributes(orgensdb)$page=="feature_page",]
#listAttributes(orgensdb)[listAttributes(orgensdb)$page=="homologs",]
if (sum(grepl("statsRes",ls()))==0) {statsRes=rawdatafile}
pepalllabcleaned=unlist(lapply(statsRes$Protein ,FUN=function(x) {
cleaned=c()
if (!is.na(x)) {
temp=strsplit(x,";")
for (k in 1:length(temp)) {
	if (sum(nchar(temp[[k]])<7)>0) {
cleaned=c(cleaned,temp[[k]][which(nchar(temp[[k]])<7)[1]])} 
	else if (grepl("-",temp[[k]][1])){
cleaned=c(cleaned,substr(temp[[k]][1],1,(gregexpr("-",temp[[k]][1])[[1]][1]-1)))}
	else {cleaned=c(cleaned,temp[[k]][1])}
#choose the first one that's uniprot 6-symbol acc., if none, return first, if first has sub category, remove sub
} } else {cleaned=c(cleaned,NA)}
return(cleaned)}))

statsRes$Protein=pepalllabcleaned
getgenename=getBM(attributes=c("uniprotsptrembl","external_gene_name","entrezgene_id","ensembl_gene_id","ensembl_peptide_id"),
filters="uniprotsptrembl",values=unique(pepalllabcleaned),mart=orgensdb)
if (nrow(getgenename)>0) {
getgenename$uniprotswissprot=NA
getswissprot=getBM(attributes=c("external_gene_name","uniprotswissprot"),
filters="external_gene_name",values=unique(getgenename$external_gene_name),mart=orgensdb)
getswissprot=getswissprot[which(getswissprot$uniprotswissprot!=""),]
getgenename$uniprotswissprot[match(getswissprot$external_gene_name,getgenename$external_gene_name)]=
getswissprot$uniprotswissprot
} else {getgenename$uniprotswissprot=character()}
getgenename2=getBM(attributes=c("uniprotswissprot","external_gene_name","entrezgene_id","ensembl_gene_id","ensembl_peptide_id"),
filters="uniprotswissprot",values=setdiff(unique(pepalllabcleaned),getgenename$uniprotsptrembl),mart=orgensdb)

#original,reviewed,entrezid,originalSitePosition,reviewedSitePosition<-need, now get first 3
idset=data.frame(original=c(getgenename$uniprotsptrembl,getgenename2$uniprotswissprot),
reviewed=c(getgenename$uniprotswissprot,getgenename2$uniprotswissprot),
entrezid=c(getgenename$entrezgene_id,getgenename2$entrezgene_id),
genename=c(getgenename$external_gene_name,getgenename2$external_gene_name),
geneid=c(getgenename$ensembl_gene_id,getgenename2$ensembl_gene_id),
peptideid=c(getgenename$ensembl_peptide_id,getgenename2$ensembl_peptide_id))
unireviewed=unique(idset$reviewed[!is.na(idset$reviewed)])
reviewedSeq=getUniProt(unireviewed)
idset$reviewSeq=reviewedSeq[match(idset$reviewed,unireviewed)]

#get all needed info in row order
accession_review_status=c()
reviewed_accession=c()
original_site=c()
reviewed_site=c()
entrez_id=c()
genename=c()
geneid=c()
peptideid=c()
for (i in 1:nrow(statsRes)) {
if (!(statsRes$Protein[i] %in% idset$original)) {
accession_review_status=c(accession_review_status,"No")
reviewed_accession=c(reviewed_accession,NA)
original_site=c(original_site,paste0(statsRes$Amino.acid[i],statsRes$Position[i]))
reviewed_site=c(reviewed_site,NA)
entrez_id=c(entrez_id,NA)
genename=c(genename,NA)
geneid=c(geneid,NA)
peptideid=c(peptideid,NA)} else {
tempinfo=idset[match(statsRes$Protein[i],idset$original),]

if (is.na(tempinfo$reviewed)) {
accession_review_status=c(accession_review_status,"No")
reviewed_accession=c(reviewed_accession,NA)
original_site=c(original_site,paste0(statsRes$Amino.acid[i],statsRes$Position[i]))
reviewed_site=c(reviewed_site,NA)} else if (tempinfo$original==tempinfo$reviewed) {
accession_review_status=c(accession_review_status,"Yes")
reviewed_accession=c(reviewed_accession,tempinfo$reviewed)
original_site=c(original_site,paste0(statsRes$Amino.acid[i],statsRes$Position[i]))
reviewed_site=c(reviewed_site,paste0(statsRes$Amino.acid[i],statsRes$Position[i])) } else {
accession_review_status=c(accession_review_status,"No")
reviewed_accession=c(reviewed_accession,tempinfo$reviewed)
original_site=c(original_site,paste0(statsRes$Amino.acid[i],statsRes$Position[i]))
if (sum(is.na(statsRes$Position))<nrow(statsRes)) {
temporiseq=statsRes$Sequence.window[i]
tempadd=0
if (grepl(";",temporiseq)) {
temporiseq=substr(temporiseq,1,gregexpr(";",temporiseq)[[1]][1]-1)}
if (grepl("_",temporiseq)) {
if (gregexpr("\\_",temporiseq)[[1]][1]==1) {
tempadd=length(gregexpr("\\_",temporiseq)[[1]])}
temporiseq=gsub("_","",temporiseq)}
sitealign=pairwiseAlignment(AAString(temporiseq),AAString(as.character(tempinfo$reviewSeq)),
type="overlap",substitutionMatrix="BLOSUM100")
if (sitealign@score>(nchar(temporiseq)-.25*nchar(temporiseq))*10) {
temppos=sitealign@subject@range@start+floor(nchar(temporiseq)/2)-tempadd } else {temppos=NULL}
if (is.null(temppos)) {reviewed_site=c(reviewed_site,NA)} else if (
substr(tempinfo$reviewSeq,temppos,temppos)!=statsRes$Amino.acid[i]) {
reviewed_site=c(reviewed_site,NA)
} else {reviewed_site=c(reviewed_site,paste0(statsRes$Amino.acid[i],temppos))} }}
if (is.na(tempinfo$peptideid)) {
peptideid=c(peptideid,NA)} else {peptideid=c(peptideid,tempinfo$peptideid)}
if (is.na(tempinfo$entrezid)) {
entrez_id=c(entrez_id,NA)} else {entrez_id=c(entrez_id,tempinfo$entrezid)}
if (is.na(tempinfo$genename)) {
genename=c(genename,NA)} else {genename=c(genename,tempinfo$genename)}
if (is.na(tempinfo$geneid)) {
geneid=c(geneid,NA)} else {geneid=c(geneid,tempinfo$geneid)}

}}

statsRes$accession_review_status=accession_review_status
statsRes$reviewed_accession=reviewed_accession
if (sum(is.na(statsRes$Position))<nrow(statsRes)) {
statsRes$original_site=original_site
statsRes$reviewed_site=reviewed_site
}
statsRes$entrez_id=entrez_id
write.csv(statsRes,"Annotation/data_results_and_reviewed_accession.csv",row.names=F)
################################################################################

#######################getting uniprot info#####################################
#Uniprot info file
#uniprotinfoyes=dlgMessage("Get Uniprot information (Can take very long to run)?",type="yesno")$res
if (uniprotinfoyes=="yes") {

dir.create("Annotation/Uniprot Information",showWarnings=F)
uniacc=unique(pepalllabcleaned)
GetExpression(uniacc, directorypath = "Annotation/Uniprot Information")
GetGeneral_Information(uniacc, directorypath = "Annotation/Uniprot Information")
GetMiscellaneous(uniacc, directorypath = "Annotation/Uniprot Information")
GetPathology_Biotech(uniacc , directorypath = "Annotation/Uniprot Information")
GetProteinFunction(uniacc, directorypath = "Annotation/Uniprot Information")
GetProteinGOInfo(uniacc, directorypath = "Annotation/Uniprot Information")
GetProteinInteractions(uniacc, directorypath = "Annotation/Uniprot Information")
GetPTM_Processing(uniacc, directorypath = "Annotation/Uniprot Information")
GetPublication(uniacc, directorypath = "Annotation/Uniprot Information")
GetSubcellular_location(uniacc, directorypath = "Annotation/Uniprot Information")
GetNamesTaxa(uniacc, directorypath = "Annotation/Uniprot Information")

}
################################################################################

#######################getting human homolog info###############################
if (phylo!="Homo sapiens") {
#humaninfoyes=dlgMessage("Get human homolog information (Can take very long to run)?",type="yesno")$res
if (humaninfoyes=="yes") {
homologattr=listAttributes(orgensdb)[listAttributes(orgensdb)$page=="homologs",]
hsaattrlist=homologattr[grep("Human",homologattr$description),]$name
accinput=list(unique(pepalllabcleaned[statsRes$accession_review_status=="No"]),
unique(pepalllabcleaned[statsRes$accession_review_status=="Yes"]))
gethsainfo=getBM(attributes=c("external_gene_name","chromosome_name","start_position","end_position","strand",hsaattrlist),filters="uniprotsptrembl",
values=accinput[[1]],mart=orgensdb)
gethsainfo2=getBM(attributes=c("external_gene_name","chromosome_name","start_position","end_position","strand",hsaattrlist),filters="uniprotswissprot",
values=accinput[[2]],mart=orgensdb)
resannot=statsRes[,c(1:6,match("accession_review_status",colnames(statsRes)):ncol(statsRes))]
resannot$geneid=geneid
resannot$peptideid=peptideid
hsainfo1=gethsainfo[match(genename[statsRes$accession_review_status=="No"],gethsainfo$external_gene_name),]
hsainfo2=gethsainfo2[match(genename[statsRes$accession_review_status=="Yes"],gethsainfo2$external_gene_name),]
hsainfo=as.data.frame(rbind(hsainfo1,hsainfo2))
hsarod=order(c(which(statsRes$accession_review_status=="No"),which(statsRes$accession_review_status=="Yes")))
hsainfoo=hsainfo[hsarod,]
hsaannot=as.data.frame(cbind(resannot,hsainfoo))

genereg=hsaannot[,match(c("geneid","peptideid","hsapiens_homolog_ensembl_gene","hsapiens_homolog_ensembl_peptide","accession_review_status"),colnames(hsaannot))]
genereg$accession=pepalllabcleaned
hsadb=useDataset(genedblist$dataset[grep("hsapiens",genedblist$dataset)],mart=ensembl)
goodrow1=which(genereg$accession_review_status=="No"&!is.na(genereg$accession))
goodrow2=which(genereg$accession_review_status=="Yes"&!is.na(genereg$accession))
orgseq=getSequence(id=genereg$accession[goodrow1],type="uniprotsptrembl",seqType="coding",mart=orgensdb,verbose=FALSE)
colnames(orgseq)[2]="accession"
orgseq2=getSequence(id=genereg$accession[goodrow2],type="uniprotswissprot",seqType="coding",mart=orgensdb,verbose=FALSE)
colnames(orgseq2)[2]="accession"
orgseqall=as.data.frame(rbind(orgseq,orgseq2))
hsaannot$coding_seq=orgseqall[match(genereg$accession,orgseqall$accession),"coding"]
orgseq3=getSequence(id=genereg$accession[goodrow1],type="uniprotsptrembl",seqType="gene_exon_intron",mart=orgensdb,verbose=FALSE)
colnames(orgseq3)[2]="accession"
orgseq4=getSequence(id=genereg$accession[goodrow2],type="uniprotswissprot",seqType="gene_exon_intron",mart=orgensdb,verbose=FALSE)
colnames(orgseq4)[2]="accession"
orgseqall2=as.data.frame(rbind(orgseq3,orgseq4))
#hsaannot$full_gene_seq=orgseqall2[match(genereg$accession,orgseqall2$accession),"gene_exon_intron"]
hsaaccinpclean=genereg$hsapiens_homolog_ensembl_peptide[!is.na(genereg$hsapiens_homolog_ensembl_peptide)&genereg$hsapiens_homolog_ensembl_peptide!=""]
hsaacc=getBM(attributes=c("uniprotsptrembl","uniprotswissprot","ensembl_peptide_id"),
filters="ensembl_peptide_id",values=hsaaccinpclean,mart=hsadb)

hsaseq=getSequence(id=hsaacc$uniprotsptrembl[hsaacc$uniprotsptrembl!=""],type="uniprotsptrembl",seqType="coding",mart=hsadb,verbose=FALSE)
hsaseq2=getSequence(id=hsaacc$uniprotswissprot[hsaacc$uniprotswissprot!=""],type="uniprotswissprot",seqType="coding",mart=hsadb,verbose=FALSE)
hsaseq3=getSequence(id=hsaacc$uniprotsptrembl[hsaacc$uniprotsptrembl!=""],type="uniprotsptrembl",seqType="gene_exon_intron",mart=hsadb,verbose=FALSE)
hsaseq4=getSequence(id=hsaacc$uniprotswissprot[hsaacc$uniprotswissprot!=""],type="uniprotswissprot",seqType="gene_exon_intron",mart=hsadb,verbose=FALSE)
colnames(hsaseq)[2]="accession"
colnames(hsaseq2)[2]="accession"
colnames(hsaseq3)[2]="accession"
colnames(hsaseq4)[2]="accession"
hsaseqall=as.data.frame(rbind(hsaseq,hsaseq2))
hsaseqall2=as.data.frame(rbind(hsaseq3,hsaseq4))
hsaannot$human_accession=hsaacc[match(hsaannot$hsapiens_homolog_ensembl_peptide,hsaacc$ensembl_peptide_id),"uniprotswissprot"]
haccblank=which(hsaannot$human_accession==""|is.na(hsaannot$human_accession))
hsaannot$human_accession[haccblank]=hsaacc[match(hsaannot$hsapiens_homolog_ensembl_peptide[haccblank],hsaacc$ensembl_peptide_id),"uniprotsptrembl"]
hsaannot$human_coding_seq=hsaseqall[match(hsaannot$human_accession,hsaseqall$accession),"coding"]
#hsaannot$human_full_gene_seq=hsaseqall2[match(hsaannot$human_accession,hsaseqall2$accession),"gene_exon_intron"]
hsaannot[hsaannot==""]=NA
#write.csv(hsaannot[,-match(c("coding_seq","human_coding_seq","full_gene_seq","human_full_gene_seq"),colnames(hsaannot))],"Annotation/human_homolog_information.csv",row.names=F)
write.csv(hsaannot[,-match(c("coding_seq","human_coding_seq"),colnames(hsaannot))],"Annotation/human_homolog_information.csv",row.names=F)
#} 
#human_homolog_alignment: 3 types, nucleotide coding_seq, peptide_seq, full_gene_seq
#protein sequence from: genereg$accession, hsaannot$human_accession,
#coding sequence: hsaannot$coding_seq, hsaannot$human_coding_seq
#full sequence: hsaannot$full_gene_seq, hsaannot$human_full_gene_seq
#file names: hsaannot$external_gene_name
#alignmentyes=dlgMessage("Generate alignment files? (For checking interested nucleotide or peptide coordinates in human. Output can be large.)",type="yesno")$res
#if (alignmentyes=="yes") {
dir.create("Annotation/Human Homolog Alignment",showWarnings=F)
#alignchoice=dlgList(c("Protein sequence alignment","Coding sequence alignment","Full gene sequence alignment"),
#title="Which alignments to generate?",multiple=T)$res

alignchoice=c("Protein sequence alignment","Coding sequence alignment")

if ("Protein sequence alignment" %in% alignchoice) {
dir.create("Annotation/Human Homolog Alignment/Protein Sequence",showWarnings=F)
alignprotidx=which(!is.na(genereg$accession)&!is.na(hsaannot$human_accession))
orgprotseq=getUniProt(genereg$accession[alignprotidx])
hsaprotseq=getUniProt(hsaannot$human_accession[alignprotidx])
nonaidx=which(!is.na(orgprotseq)&!is.na(hsaprotseq))
alignprotidx=alignprotidx[nonaidx]
orgprotseq=orgprotseq[nonaidx]
hsaprotseq=hsaprotseq[nonaidx]
for (i in 1:length(alignprotidx)) {
idperc=hsaannot$hsapiens_homolog_perc_id[alignprotidx][i]
if (!is.na(idperc)) {
if (idperc>90) {smatrix="BLOSUM100"} else if (idperc>80) {smatrix="BLOSUM80"} 
else if (idperc<45) {smatrix="BLOSUM45"} else {smatrix="BLOSUM62"}} else {smatrix="BLOSUM62"}
#temphomologalign=pairwiseAlignment(AAString(orgprotseq[[i]]),
#AAString(hsaprotseq[[i]]),type="overlap",substitutionMatrix=smatrix)
temphomologalign=pairwiseAlignment(AAString(gsub("U","C",orgprotseq[[i]])),
AAString(gsub("U","C",hsaprotseq[[i]])),type="overlap",substitutionMatrix=smatrix)
if (temphomologalign@score!=0) {
tempalignfname=hsaannot$external_gene_name[alignprotidx[i]]
if (is.na(tempalignfname)) {
tempalignfname=genereg$accession[alignprotidx[i]]}
writePairwiseAlignments(temphomologalign, file=paste0("Annotation/Human Homolog Alignment/Protein Sequence/",tempalignfname,".txt"))}
}}

if ("Coding sequence alignment" %in% alignchoice) {
dir.create("Annotation/Human Homolog Alignment/Coding Sequence",showWarnings=F)
aligncodingidx=which(!is.na(hsaannot$coding_seq)&!is.na(hsaannot$human_coding_seq))
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
for (i in aligncodingidx) {
tempcodingseq=gsub("[^ACGT]","",hsaannot$coding_seq[i])
temphsacodingseq=gsub("[^ACGT]","",hsaannot$human_coding_seq[i])
temphomologalign=pairwiseAlignment(AAString(tempcodingseq),AAString(temphsacodingseq),
type="overlap",substitutionMatrix = mat,gapOpening = 5, gapExtension = 2)
if (temphomologalign@score!=0) {
tempalignfname=hsaannot$external_gene_name[i]
if (is.na(tempalignfname)) {
tempalignfname=genereg$accession[i]}
writePairwiseAlignments(temphomologalign, file=paste0("Annotation/Human Homolog Alignment/Coding Sequence/",tempalignfname,".txt"))}
}}

#if ("Full gene sequence alignment" %in% alignchoice) {
#dir.create("Annotation/Human Homolog Alignment/Full Gene Sequence",showWarnings=F)
#aligngeneidx=which(!is.na(hsaannot$full_gene_seq)&!is.na(hsaannot$human_full_gene_seq))
#mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
#for (i in aligngeneidx) {
#temphomologalign=pairwiseAlignment(AAString(hsaannot$full_gene_seq[i]),
#AAString(hsaannot$human_full_gene_seq[i]),type="global",substitutionMatrix = mat,gapOpening = 5, gapExtension = 2)
#tempalignfname=hsaannot$external_gene_name[i]
#if (is.na(tempalignfname)) {
#tempalignfname=genereg$accession[i]}
#writePairwiseAlignments(temphomologalign, file=paste0("Annotation/Human Homolog Alignment/Full Gene Sequence/",tempalignfname,".txt"))
#}}
}
} 
################################################################################
} else {msgBox("Organism not supported for annotation analysis. Moving on!")} 
