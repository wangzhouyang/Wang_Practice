args <- commandArgs(T)
if (length(args) != 2) {stop("Rscript wilcox.test.KW.R [input] [output]
	input: dat_profile
	output: output file")
}

dat.dir <- args[1]
out.dir <- args[2]

dat <- as.data.frame(read.table(dat.dir,header=T,row.names=1))

ncol <- length(colnames(dat))
nrow <- length(rownames(dat))

res <- dat

for(i in 1:nrow){
	line <- as.numeric(dat[i,])
	for(j in 1:ncol){head 
		res[i,j] = (dat[i,j] - min(line))/(max(line) - min(line))
	}
}

write.table(res,out.dir,quote=F,sep="\t",col.names=NA,row.names=T)
