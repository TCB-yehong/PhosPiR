#######################Set up source location###################################
if(!require(svDialogs)){install.packages("svDialogs")}
library(svDialogs)
dlgMessage("Please select run.R file again")
codedr=file.choose()
codepth=substr(codedr,1,gregexpr("\\\\",codedr)[[1]][length(gregexpr("\\\\",codedr)[[1]])])
################################################################################

#######################Running analysis#########################################
source(paste0(codepth,"inputManage.R"))

if (pipestep=="Generate Input File") { 
source(paste0(codepth,"inputGeneration.R"))
dlgMessage("All done, thank you!")
} else if (pipestep=="Annotation") {
source(paste0(codepth,"annotation.R"))
dlgMessage("All done, thank you!")

} else if (pipestep=="Overview Figures") {
source(paste0(codepth,"overviewFigure.R"))
dlgMessage("All done, thank you!")

} else if (pipestep=="Differential Analysis") {
source(paste0(codepth,"statisticalAnalysis.R"))
statsRes=statFunc(rawdatafile,group,grpcompare)
dlgMessage("All done, thank you!")

} else {
source(paste0(codepth,"overviewFigure.R"))
source(paste0(codepth,"statisticalAnalysis.R"))
statsRes=statFunc(rawdatafile,group,grpcompare)
siglistFiles=rownames(file.info(list.files("Statistical Analysis/Significant Lists")))
source(paste0(codepth,"annotation.R"))
source(paste0(codepth,"enrichment2.R"))
if (sum(is.na(statsRes$Position))<nrow(statsRes)) {
source(paste0(codepth,"enrichmentPTMSEA.R"))
source(paste0(codepth,"kinase.R"))
}
source(paste0(codepth,"Network.R"))
dlgMessage("All done, thank you!")
}
################################################################################
