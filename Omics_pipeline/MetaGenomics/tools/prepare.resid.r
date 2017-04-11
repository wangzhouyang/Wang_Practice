argv<-commandArgs(T)
if(length(argv<2){
	stop("Rscript $0 <input> <output>")
}
fIN <- argv[1]
fOUT <- argv[2]

pr <- read.table(fIN,header=T)

#library(car)
#source("/home/xiehailiang/bin/heritability/heritability_accessory.R")
pr <- as.matrix(pr)
pr <- pr[rowSums(pr)!=0, ]
pr<-na.omit(pr)
asin.dat<-t(apply(pr, 1, function(x) {x<-as.numeric(x);x<-sqrt(x);x<-asin(x)}))
d<-c()
fun=function(x){    
	x<-as.numeric(x)
	a1<-0
	a2<-0
	x1<-x[seq(1,length(x),by=2)]
	if(min(x1)!=max(x1)){
		a1<-shapiro.test(x1)$p.value
	}
	x2<-x[seq(2,length(x),by=2)]
	if(min(x2)!=max(x2)) {
		a2<-shapiro.test(x2)$p.value
	}
	d<-c(a1,a2)     
} 
asin.dat.p <-t(apply(asin.dat, 1, fun))
	bat1.p<-apply(asin.dat, 1, function(x) {x<-as.numeric(x);bartlett.test(x,g=gl(2,1,length(x)))$p.value})
	bat2.p<-apply(asin.dat, 1, function(x) {x<-as.numeric(x);bartlett.test(x,g=gl(length(x)/2,2,length(x)))$p.value})
	ap<-cbind(asin.dat.p,bat1.p,bat2.p)
p<-apply(ap, 1, function(x) {min(x)})
	pr <-pr[p>0.05,]
write.table(pr,file=argv[2],quote=FALSE)
