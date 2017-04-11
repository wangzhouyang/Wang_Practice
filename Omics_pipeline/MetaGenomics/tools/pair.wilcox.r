argv<-commandArgs(T)
dat <- as.matrix(read.table(argv[1],head=T))
dat <- dat[rowSums(dat)!=0,]
num <- nrow(dat)


fun=function(data)
{
x=data[seq(1,length(data),by=2)]
y=data[seq(2,length(data),by=2)]
wilcox.test(x,y,paired=T,exact=F)$p.value
}
table<-apply(dat, 1, fun)
write.table(table,file=argv[2])
