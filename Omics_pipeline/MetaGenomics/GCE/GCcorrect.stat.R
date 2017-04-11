args <- commandArgs(T)
cohort <- args[1]
sample <- args[2]
bac <- args[3]
ref <- args[4]
dir <- args[5]
if(length(args)<5){
	stop("Rscript $0 cohort sample bac ref dir\n")
}
###
#setwd('E:/SFTP/HCC')
#######

########
library("GCcorrect")
#bug fix
source("/ifs1/ST_MD/USER/fangchao/pipeline/Omics_pipeline/MetaGenomics/GCE/fix.function.R")
source("/ifs1/ST_MD/USER/fangchao/pipeline/Omics_pipeline/MetaGenomics/GCE/estimate.R")
dyn.load("/ifs1/ST_MD/USER/fangchao/pipeline/Omics_pipeline/MetaGenomics/GCE/makeGC.so")

if(!file.exists(dir)){dir.create(dir)}
if(!file.exists(paste0(dir,"/RData"))){dir.create(paste0(dir,"/RData"))}
if(!file.exists(paste0(dir,"/gcchr"))){dir.create(paste0(dir,"/gcchr"))}
if(!file.exists(paste0(dir,"/Results/"))){dir.create(paste0(dir,"/Results/"))}
sdir <- paste0(dir,"/Results/",sample)
if(!file.exists(sdir)){dir.create(sdir)}

chr1 = prepareChrom(paste0(dir,"/RData/",bac,".RData"),reference=paste0(ref,"/sp.",bac,".seq"),
repeats = paste0(ref,"/sp.",bac,".rep"),details = list(bac,1,50))

filename <- paste0(cohort,"/",sample,"/",bac,".pos")


reads <- rep(0,20)
reads[21:25] <- c(0,0,cohort,bac,sample)
rdat <- data.frame(sample=reads)
colnames(rdat) <- sample

if(!file.info(filename)$size){
	dat1_for = matrix(1,nc=2,nr=2)
}else{
	dat1_for = as.matrix(read.table(filename)[,2:3])
	if (nrow(dat1_for)<2){
#		dat1_for = matrix(1,nc=2,nr=2)
		print(paste0(filename,": Task finished.1"),quote=F)
		stop(paste0(filename," less than 2 points, abort."))
	}
}

dat1_rev = matrix(1,nc=2,nr=2);

dat1 = prepareReads(dat1_for, dat1_rev,chr_len = length(chr1$chrline), 1,40,'F')
if(length(dat1[["for10K"]])<2){
	print(paste0(filename,": Task finished.2"),quote=F)
	stop(paste0(filename,"reads number less than 2, abort."))
}

pdf(paste0(sdir,"/",bac,".",sample,".plot1.pdf"))
#print(
	basicForPlots(dat1,types = c(1,2,3),cohort = cohort,bac = bac,sample = sample)
#	)
dev.off()

forbin = sprintf("for%s", "10K")
dfor <- density(dat1[[forbin]],from = 0)
depth <- dfor$x[which(dfor$y==max(dfor$y))[1]]

set.seed(1000)

# removing zero stretches and pileups
	useSamp = !logical(length(chr1$isgc)) 
#zeros_10K = findLongZeros(dat1$for10K,2)
#	if(!is.null(ncol(zeros_10K)) && ncol(zeros_10K)>0){
#		for (i in 1:ncol(zeros_10K)) {
#			useSamp[((10000*(zeros_10K[1,i]-1)) + 1) : (10000*(zeros_10K[2,i]))] = F
#		}
#		useSamp = useSamp[1:length(chr1$isgc)]
#
#			high_10K = which(dat1$for10K>600)
#			for (i in high_10K) {
#				useSamp[10000*(i-1) + 1:10000] = F
#			}
#		useSamp = useSamp[1:length(chr1$isgc)]
#	} 
#####
	if(!length(which(useSamp==T))){
		reads <- rep(0,20)
			reads[21:25] <- c(0,0,cohort,bac,sample)
			rdat <- data.frame(sample=reads)
			colnames(rdat) <- sample
			return(rdat)
	}

sampCh1 = sampleChrom(chr1,dat1,n=length(chr1$chrline),margin = 0,len_range = c(), memoryopt = TRUE,useSamp = useSamp) 

	margin = 2

	forw <- dat1$forw[order(dat1$forw[,1]),]
#begdata = makeGCLens(chr1$isgc,forw,sampline = sampCh1$singleLocSamp,minlen = 0,maxlens=100,margin=5,max_frag_for_loc=10)
begdata = makeGCL2(chr1$isgc,forw,sampline = sampCh1$singleLocSamp,minlen = 1,maxlens=300,margin=2,max_frag_for_loc=10)

write.table(begdata$locs,paste0(sdir,"/",bac,".",sample,".locs.xls"),col.names=NA,sep="\t",quote=F)
write.table(begdata$frag,paste0(sdir,"/",bac,".",sample,".frag.xls"),col.names=NA,sep="\t",quote=F)
####
tvs = scoreGCLens(begdata, maxlen=300, minlen =  0,scale= T)
best_size = which.max(tvs)-1

pdf(paste0(sdir,"/",bac,".",sample,".plot2.pdf"))
#print(
	plotGCLens(tvs,lw=2,lt=1,ylim=c(0,0.15))
	title(main=sprintf("%s | %s | %s",cohort,bac,sample))
	abline(v=best_size)
	text(x=best_size,y=0,best_size,adj=-0.5)
#)
dev.off()

####
gcsize = best_size # I got 555 for my sample of the full chromosome
gcline = prepareGCChrom(chr1,gcsize,filename = paste0(dir,"/gcchr/",bac,".",sample,".gcline"))


##

cMeans= getCondMean(gcline[sampCh1$singleLocSamp+margin],sampCh1$forsamped,cutoff = 4,jump = 1,norm = FALSE)
pdf(paste0(sdir,"/",bac,".",sample,".plot3.pdf"))
#print(
	c.plotCondMean(cMean = cMeans,ci = TRUE,normRange = gcsize,meanLine=TRUE,lt = 1,col=4,cohort = cohort,bac = bac,sample = sample)
#	title(main=sprintf("%s | %s | %s \nTV = %g",cohort,bac,sample,cMean[[3]]))
#	)
dev.off()

library(plyr)
reads <- rep(0,20)
index <- round(as.numeric(names(cMeans$nreads))/gcsize/5,2)*5
tmp <- ddply(data.frame(index=index,nreads=cMeans$nreads,cts=cMeans$cts),
			"index",summarize,nreads=sum(nreads),cts=sum(cts),rate=sum(nreads)/sum(cts))

out <- data.frame(index=seq(0.05,1,0.05),nreads=reads,cts=reads,rate=reads)
if(tmp[1,1]==0){
		out[tmp$index*20,] <- tmp[-1,]
	}else{
		out[tmp$index*20,] <- tmp
	}
out[21,] <- c("gcsize",gcsize,gcsize,1)

write.table(out,paste0(sdir,"/",bac,".",sample,".stat.xls"),col.names=NA,sep="\t",quote=F)

print(paste0(filename,": Task all finished."),quote=F)
