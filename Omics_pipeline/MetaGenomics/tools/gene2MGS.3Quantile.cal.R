args <- commandArgs(T)
if(length(args)<5){
	stop("Rscript $0 [MGS reference] [gene abundance] [core num] [quantile] [MGS output]
			MGS reference is a file using NGS id as index, each of which following a list of gene ids belong to
			abundance file contained 2 columns (gene index\tabundance) with a title
			set core num for parallel calculation
			quantile should be 2(median) or 3, which depends on you)")
}

ref <- read.table(args[1],header=F,row.names=1,fill=NA)
input  <- read.table(args[2],header=T,row.names=1,sep="\t")
Cn <- args[3]
Qt  <- as.numeric(args[4])
out <- as.numeric(args[5])

valFun <- function(x,y,q){as.numeric(quantile(y[x[!is.na(x)],1],q*0.25))}

library(parallel)
cl <- makeCluster(spec=Cn,type="FORK")
res <- as.matrix(parApply(cl,ref,1,FUN=valFun,input,Qt))
#res <- as.matrix(apply(ref,1,FUN=valFun,input,Qt))
stopCluster(cl)

write.table(res,out,col.names =F,quote=F,sep="\t")
