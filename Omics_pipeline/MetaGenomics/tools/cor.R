args <- commandArgs(T)
if (length(args) < 1) {
	stop("Rscript $0 [function] [arguments]
		Rscript $0 cal <file1> <file2> <suffix> <output> # calculate cor for 2 individual files
		Rscript $0 combine <cor list> <output> # transform cor list into matrix
		Rscript $0 samllcal <file> <core> <output> # for calcualte cor for small profile")
}
open <- function(x){
	if(length(grep(".gz",x,fixed=T))){
		dat <- read.table(gzfile(x),header=T,row.names=1,sep="\t");
	}else{
		dat <- read.table(x,header=T,row.names=1,sep="\t");
	}
	return(dat)
}
cal <- function(file1,file2,suffix,out){
	base1 <- sub(suffix,"",basename(file1));
	base2 <- sub(suffix,"",basename(file2));

	dat1 <- read.table(gzfile(file1),header=T,row.names=1,sep="\t");
	dat2 <- read.table(gzfile(file2),header=T,row.names=1,sep="\t");
	res <- cor.test(dat1[,3],dat2[,3]);
	write.table(list(base1,base2,res$p.value,as.numeric(res$estimate)),out,sep="\t",quote=F,row.names=F,col.names=F);
}

combine <- function(file,out){
	dat <- read.table(file,header=F,sep="\t");
	name <- levels(dat$V1)
	num <- length(name)
	outP <- matrix(ncol=num,nrow=num,NA)
	colnames(outP) <- name
	rownames(outP) <- name
	outE <- outP
	for(k in 1:nrow(dat)){
		i <- which(name==dat[k,1])
		j <- which(name==dat[k,2])
		p <- as.numeric(dat[k,3])
		e <- as.numeric(dat[k,4])
		outP[i,j] <- p
		outP[j,i] <- p
		outE[i,j] <- e
		outE[j,i] <- e
	}
	write.table(outP,paste0(out,".pval"),sep="\t",quote=F,row.names=T,col.names=NA);
	write.table(outE,paste0(out,".esti"),sep="\t",quote=F,row.names=T,col.names=NA);
}

smallcal <- function(file,core,out){
	dat <- open(file);
	num <- ncol(dat);
	list <- NULL;
	for (i in 1:num){
		for(j in i:num){
			list <- rbind(list,c(i,j))
		}
	}
	pcal <- function(x,list,dat){
		res <- cor.test(dat[,as.numeric(list[x,1])],dat[,as.numeric(list[x,2])])
		return(c(list[x,1],list[x,2],res$p.value,as.numeric(res$estimate)))
	}
	library(parallel)
	cl <- makeCluster(core,"FORK")
	res <- parLapply(cl,1:nrow(list),pcal,list,dat)
	res.df <- do.call('rbind',res)
	stopCluster(cl)
	outP <- matrix(ncol=num,nrow=num,NA);
	colnames(outP) <- colnames(dat)
	rownames(outP) <- colnames(dat)
	outE <- outP;
	for(k in 1:nrow(res.df)){
		i <- as.numeric(res.df[k,1])
		j <- as.numeric(res.df[k,2])
		p <- as.numeric(res.df[k,3])
		e <- as.numeric(res.df[k,4])
		outP[i,j] <- p
		outP[j,i] <- p
		outE[i,j] <- e
		outE[j,i] <- e
		res.df[k,1] <- colnames(dat)[i]
		res.df[k,2] <- colnames(dat)[j]
	}
	write.table(outP,paste0(out,".pval"),sep="\t",quote=F,row.names=T,col.names=NA);
	write.table(outE,paste0(out,".esti"),sep="\t",quote=F,row.names=T,col.names=NA);
	write.table(res.df,paste0(out,".list"),sep="\t",quote=F,row.names=F,col.names=F);
}

if(args[1] == "smallcal"){
	smallcal(args[2], args[3], args[4])
}else if (args[1] == "cal"){
	cal(args[2], args[3], args[4], args[5])
}else if (args[1] == "combine"){
	combine(args[2], args[3])
}
