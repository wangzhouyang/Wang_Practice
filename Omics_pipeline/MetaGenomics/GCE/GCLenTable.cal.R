args <- commandArgs(T)
if(length(args)<2){
	stop("Rscript $0 Target(sample)Directory outputPrefix\n")
}

dir <- args[1]
outpfx <- args[2]

f <- paste0("files.",dir,".stat")
dat <- read.table(f)
locs <- matrix(nrow = 700,ncol = 700)
colnames(locs) <- seq(0,699)
rownames(locs) <- seq(0,699)
locs[is.na(locs)] <- 0
frag <- locs
for (i in 1:nrow(dat)) {
	if (dat[i,1] == 6) {
		tlocs <- as.matrix(read.table(paste0(dir,"/",dat[i,2],".locs.xls")))
			tfrag <-
			as.matrix(read.table(paste0(dir,"/",dat[i,2],".frag.xls")))
			locs = locs + tlocs
			frag = frag + tfrag
	}
}
write.table(locs,paste0(outpfx,".locs.xls"),sep="\t",quote=F,col.names=NA)
write.table(frag,paste0(outpfx,".frag.xls"),sep="\t",quote=F,col.names=NA)
