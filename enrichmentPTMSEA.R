################################################################################
#code comes from broadinstitute/ssGSEA2.0
#content has been modified to fit this pipeline
#all source code necessary to run ssGSEA2.0 has been copied in the pipeline folder 
#they are files 'io.R', 'utils.R' and 'ssGSEA2.0.R'
#For the original code and more details on ssGSEA2.0,
#please visit https://github.com/broadinstitute/ssGSEA2.0

#License

#ssGSEA2.0 is distributed under the following BSD-style license:

#Copyright © 2018 – 2020 The Broad Institute, Inc.  All rights reserved.

#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#3. Neither the names of the copyright holders nor the names of their contributors may be used to endorse or promote products derived from this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED “AS IS.”  THE COPYRIGHT HOLDERS MAKE NO EXPRESS OR IMPLIED REPRESENTATIONS OR WARRANTIES OF ANY KIND REGARDING THE SOFTWARE AND COPYRIGHT, INCLUDING, BUT NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, CONFORMITY WITH ANY DOCUMENTATION, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF, HAVE REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF SUCH DAMAGE.

#If, by operation of law or otherwise, any of the aforementioned warranty disclaimers are determined inapplicable, your sole remedy, regardless of the form of action, including, but not limited to, negligence and strict liability, shall be replacement of the software with an updated version if one exists.

#Development of ssGSEA2.0 was supported by grant U24CA210979 from the National Cancer Institute Clinical Proteomic Tumor Analysis Consortium (CPTAC).  In addition, ssGSEA2.0 is distributed, in part, under and subject to the provisions of license(s) for: cmap / cmapR, © 2017 Connectivity Map (CMap) at the Broad Institute, Inc. (https://github.com/cmap/cmapR/blob/master/LICENSE, BSD 3-Clause License). All rights reserved.
################################################################################

if (!phylo %in% c("Homo sapiens","Mus musculus","Rattus norvegicus")) {
msgBox("Organism not supported for PTM-SEA analysis. Moving on!")
} else {
gmtdb=list.files(path=codepth,pattern="*.gmt")
if (phylo=="Homo sapiens") {
gene.set.databases=paste(codepth,gmtdb[grep("human",gmtdb)],sep="")
} else if (phylo=="Mus musculus") {
gene.set.databases=paste(codepth,gmtdb[grep("mouse",gmtdb)],sep="")
} else {
gene.set.databases=paste(codepth,gmtdb[grep("rat",gmtdb)],sep="")
}

dir.create("Enrichment/PhosphoSite enrichment",showWarning=F)
if(!require(cmapR)){BiocManager::install("cmapR",update=F,ask=F)}
library(cmapR)
if (nchar(statsRes$Sequence.window[1])>15) {
seqlenm=(nchar(statsRes$Sequence.window[1])-15)/2
seq.window2=substr(statsRes$Sequence.window,start=seqlenm+1,stop=nchar(statsRes$Sequence.window[1])-seqlenm)} else {
seq.window2=statsRes$Sequence.window}
ptmannot=data.frame(
geneSymbol=statsRes$Gene.names,
accession_number=statsRes$Protein,
id=paste0(seq.window2,"-p"),
id.uniprot=paste0(statsRes$Protein,";",statsRes$Amino.acid,statsRes$Position,"-p"),
id.flanking=seq.window2)
#ptmvalue=as.matrix(statsRes[,7:ncol(rawdatafile)])
#rownames(ptmvalue)=ptmannot$id
#ptmvalue=ptmvalue[!duplicated(ptmannot$id),]
ptmannotnodup=ptmannot[!duplicated(ptmannot$id),]
#menrich=new("GCT",mat=ptmvalue,rdesc=ptmannotnodup)
#write_gct(menrich, "Enrichment/PhosphoSite enrichment/gct_input")

for (i in 1:length(grpcompare)) {
tempptmvalue=as.matrix(statsRes[,grep(paste0("Comparison",i,"_"),colnames(statsRes))])
tempptmvalue=tempptmvalue[!duplicated(ptmannot$id),]
rownames(tempptmvalue)=ptmannotnodup$id
tempmenrich=new("GCT",mat=tempptmvalue,rdesc=ptmannotnodup)
write_gct(tempmenrich, paste0("Enrichment/PhosphoSite enrichment/gct_input_comparison",i))}

#rm(list=ls())
script.dir <- codepth ## get folder the script is located in
os <- Sys.info()['sysname'] ## determine operating system
if (!require("pacman")) install.packages ("pacman")
library(pacman)

## ##########################################################
##  define parameters below:
## ##########################################################

## ssGSEA / PTM-SEA parameters
sample.norm.type    = "rank"              ## "rank", "log", "log.rank", "none" 
weight              = 0.75                ## value between 0 (no weighting) and 1 (actual data counts)
statistic           = "area.under.RES"    ## "Kolmogorov-Smirnov"
output.score.type   = "NES"               ## 'ES' or 'NES'
nperm               = 1e3                 ## No. of permutations
min.overlap         = 3                  ## minimal overlap between gene set and data
correl.type         = "z.score"           ## 'rank', 'z.score', 'symm.rank'
par                 = T                   ## use 'doParallel' package?
spare.cores         = 1                   ## No. of cores to leave idle
export.signat.gct   = T                   ## if TRUE gene set GCT files will be exported 
extended.output     = T                   ## if TRUE the GCT files will contain stats on gene set overlaps etc.   

## #####################################################################
##   end paramaters
## - in a perfect world users don't have to worry about the stuff below...
## #####################################################################

## #################################
## directory with gct files
gct.dir.ok=T
gct.dir=paste0(getwd(),"/Enrichment/PhosphoSite enrichment")

## directory to write output
out.dir <- gct.dir

## MSigDB
db.ok=T

## ######################################################################
##                          START
## ######################################################################
#source(paste(script.dir, 'src/ssGSEA2.0.R', sep='/'))
source(file.path(script.dir, 'ssGSEA2.0.R'))

## #############################################
## prepare output folder
setwd(out.dir)

date.str <- paste(sub(' .*', '', Sys.time()), sep='_')
dir.create(date.str)
setwd(date.str)


## #############################################
## import signature database
signat.all <- unlist(lapply(gene.set.databases, readLines))
signat.all <- strsplit(signat.all, '\t')
names(signat.all) <- sapply(signat.all, function(x)x[1])
signat.all <- lapply(signat.all, function(x) x[-c(1,2)])

## save parameters used for ssGSEA
param.str = c(
    paste('##', Sys.time()),
    paste('gct.directory:', gct.dir, sep='\t'),
    paste('output.directory:', out.dir, sep='\t'),
    paste('gene.set.database:',gene.set.databases, sep='\t'),
    paste('sample.norm.type:', sample.norm.type, sep='\t'),
    paste('weight:', weight, sep='\t'),
    paste('statistic:', statistic, sep='\t'),
    paste('output.score.type', output.score.type, sep='\t'),
    paste('nperm:', nperm, sep='\t'),
    paste('min.overlap:', min.overlap, sep='\t'),
    paste('correl.type:', correl.type, sep='\t'),
    paste('run.parallel:', par, sep='\t')
   )
writeLines(param.str, con='parameters.txt')


## identify all gct files
gct.files <- dir(gct.dir, pattern='\\.gct$', full.names=T)
names(gct.files) <- paste(  sub('\\.gct$', '', sub('.*/','', gct.files)), 'ssGSEA', sep='_' )

#debug(ssGSEA2)

## #####################################
## loop over gct files and run ssGSEA
for(i in names(gct.files)){


    ## create sub folders if more than one gct file was found
    if(length(gct.files) > 1){
        subdir=sub(',|\\.|:|;|/', '_', i)
        dir.create(subdir)
        setwd(subdir)
    }

    ## ########################################
    ## ssGSEA

    ## input data set
    input.ds <- gct.files[i]

    cat('Running ssSGEA on:', sub('.*/', '', input.ds), '\n\n')

    ## run ssGSEA
    gsea.res <- tryCatch(ssGSEA2(input.ds, gene.set.databases=gene.set.databases, sample.norm.type=sample.norm.type, weight=weight,statistic=statistic, output.score.type = output.score.type, nperm  = nperm, min.overlap  = min.overlap, correl.type = correl.type, output.prefix = paste(i), par=par, 
                        spare.cores=spare.cores, param.file=F, export.signat.gct = export.signat.gct, extended.output = extended.output ),error=function(e) "An error has occurred in the calculation, the result came out empty.")

    ## save object
    save(gsea.res, file=paste(i, '.RData', sep=''))

    ## #########################################################
    ##              rank plots
    ## #########################################################
    ## flag to indicate presence of duplicated ids in GCT file
    ## e.g. in case of gene-centric-redundant signature analysis
    dups=F
    if(file.exists( sub('\\.gct', '_unique.gct', input.ds)))
      dups <- T
      
    ## input dataset
    if(dups)
       input.ds <-sub('\\.gct', '_unique.gct', input.ds)
    input.gct <- parse.gctx(input.ds) 
    
    ## gene/site ids
    gn.input <- input.gct@rid
    if(dups)
      gn.input <-  sub('_[0-9]{1,4}$', '', gn.input)
    
    ## sample names
    all.samp <- input.gct@cid
    
    ## expression data only
    input <- input.gct@mat
    
    if (file.exists(paste( i, '-scores(_[0-9]*x[0-9*]|)', '.gct', sep=''))) {	   
    ## import enrichment scores and p-values
    gsea.score.gct <- parse.gctx(dir('.', pattern=paste( i, '-scores(_[0-9]*x[0-9*]|)', '.gct', sep=''))) 
    gsea.score <- gsea.score.gct@mat
    gsea.pval.gct <- parse.gctx(dir('.', pattern=paste( i,  '-fdr-pvalues(_[0-9]*x[0-9*]|)', sep=''))) 
    gsea.pval <- gsea.pval.gct@mat
    
    ## gene set names
    all.gs <- rownames(gsea.score)

    ## keep only scored signatures
    signat <- signat.all[all.gs]

    ## create sub-folder
    dir.create('rank-plots')
    
    ## loop over gene sets
    for(gs in 1:length(all.gs)){

        gs.name <- all.gs[gs]

        pdf(paste('rank-plots/', make.names( chopString( gsub('\\:|\\/\\\t', ' ', gs.name), nChar=20, add.dots=F)) ,'_2.pdf', sep=''), 9.5, 9.5)
        par(mfrow=c(3, 3))
        for(samp in 1:length(all.samp)){

            ## extract results
            samp.name <- all.samp[samp]

            ## gsea results
            score <- gsea.score[gs.name, samp.name]
            pval <- gsea.pval[gs.name, samp.name]

            ## extract data
            data.expr <- input[, samp.name ]

            valid.idx <- which( !(is.na( data.expr ) | is.infinite(data.expr)) )

            data.expr <- data.expr[ valid.idx ]
            gn <- gn.input[ valid.idx ]

            ## order
            ord.idx <- order(data.expr, decreasing=T)

            ##gn <- row.names(input)[ord.idx]
            gn <- gn[ ord.idx ]
            data.expr <- data.expr[ ord.idx ]


            plot( data.expr, pch=20, col='darkgrey', lwd=4, type='l', xlab='Rank', ylab='Expression', main=paste(gs.name, samp.name, sep='\n'), ylim=range(data.expr), yaxs='i')
						abline(h=0, lty='dashed', lwd=2, col='grey70')

            ## #########################################################
            ##  ptm signatures?
            if(length(grep(';u$|;d$', signat[[gs.name]], value=T)) > 0){

                ## locations
                gsea.tmp.u <- sub(';u$','',grep(';u$', signat[[gs.name]], value=T))
                loc.u <- na.omit(match(gsea.tmp.u, gn))

                gsea.tmp.d <- sub(';d$','',grep(';d$',  signat[[gs.name]], value=T))
                loc.d <- na.omit(match(gsea.tmp.d, gn))

                if(!is.null(loc.u)){

                    rug(loc.u, col='darkred', side=3, lwd=3, ticksize=0.02)
                    points(loc.u, data.expr[loc.u], col=my.col2rgb('darkred',  150), pch=16, cex=2)
                }
                if(!is.null(loc.d)){
                    rug(loc.d, col='darkblue', side=1, lwd=3, ticksize=0.02)
                    points(loc.d, data.expr[loc.d], col=my.col2rgb('darkblue',  150), pch=16, cex=2)

                }
                ## some info
                legend('bottom', legend=c(paste('No. down-regulated in signature:', length(grep(';d$', signat[[gs.name]]))),
                                          paste('No. found in data set:', length(loc.d))
                                          ), inset=.05, bty='n', text.col='darkblue')

                legend('top', legend=c(paste('No. up-regulated in signature:', length(grep(';u$', signat[[gs.name]]))),
                                       paste('No. found in data set:', length(loc.u))
                                       ), inset=.05, bty='n', text.col='darkred')
            } else {## end if signature

                ## ####################################################
                ## regular gene set
                loc <- which(gn %in% signat[[gs.name]])
                rug(loc, col=my.col2rgb('darkred',  50), side=3, lwd=2, ticksize=0.02)
								points(loc, data.expr[loc], col=my.col2rgb('darkred',  150), pch=16, cex=2)

                ## box plot
                loc.quart <- quantile(loc)
                rug(loc.quart, col='darkblue', side=3, lwd=2, ticksize=0.03)
                rect( loc.quart[2], max(data.expr)-0.04*max(data.expr-min(data.expr)), loc.quart[4], max(data.expr), border='darkblue', lwd=2, col=NA )
                rect( loc.quart[1], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[2], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )
                rect( loc.quart[4], max(data.expr)-0.02*max(data.expr-min(data.expr)), loc.quart[5], max(data.expr)-0.02*max(data.expr-min(data.expr)), border='darkblue', lwd=2 )

                ## some info
                legend('bottom', legend=c(paste('No. in signature:', length( signat[[gs.name]])),
                                          paste('No. found in data set (non-redund.):', sum(signat[[gs.name]] %in% gn)),
                                          paste('No. found in data set (redundant):', length(loc))
                                       ), inset=.05, bty='n', text.col='darkred')

            }
            legend('right', legend=paste('NES=', round(score, 3), ' (p.adj=', round(pval, 5), ')', sep=''), bty='n', inset=.2, cex=1.5)
	}
	par(mfrow=c(1, 1))
	dev.off()
    }} ## end loop over gene sets

    if(length(gct.files) > 1)
        setwd('..')
  
    if(dups)
      file.remove(input.ds)
}

setwd(wddflt)}