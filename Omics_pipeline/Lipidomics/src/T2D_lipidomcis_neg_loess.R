options(bitmapType="cairo")
library(metaX)
para <- new("metaXpara")
pfile <- "T2D_lipidomics_neg_LipidBlast_and_HMDB_measurement.csv"
sfile <- "DM_lipid_neg_sample.list"
#sfile <- "zaochan_globle_pos(batch1-6-QC29)_remove3.list"
idres <- "T2D_lipidomics_neg_LipidBlast_and_HMDB_identification.csv"

para@outdir <- "20151118T2D_neg_loess_nor"
para@prefix <- "neg"

para@sampleListFile <- sfile
para@ratioPairs <- "T2D:NGT;Pre-DM:NGT;T2D:Pre-DM"
para <- importDataFromQI(para,file=pfile)
#para <- removeSample(para,rsamples=c("PDB13AY02035","PDB13CP00599","13B0149944"))

plsdaPara <- new("plsDAPara")
plsdaPara@scale = "pareto"
plsdaPara@kfold = 5
plsdaPara@cpu = 1
p <- metaXpipe(para,plsdaPara=plsdaPara,cvFilter=0.3,idres=idres,qcsc=TRUE,scale="pareto",remveOutlier=TRUE,nor.method="none",pclean = FALSE,doROC=FALSE)
#save(p,file="neg_p1.rda")

#sessionInfo()
