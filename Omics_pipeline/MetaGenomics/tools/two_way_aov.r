argv <- commandArgs(T)
test <- read.table(argv[1],header=T)
test<-as.matrix(test)
test<-test[rowSums(test)!=0,]
fun=function(y1)
 {yy<-as.numeric(y1)
 dat <- data.frame(yy, a = gl(2, 1,length(yy)), b= gl(length(yy)/2,2,length(yy)))
 myaov<-aov(yy~a+b, data=dat)
s<-summary(myaov)
sa<-s[[1]][[2]][1]#Sa
sb<-s[[1]][[2]][2]#Sb

msa<-s[[1]][[3]][1]#mSa
msb<-s[[1]][[3]][2]#mSb

fa<-s[[1]][[4]][1]#fa

p<-s[[1]][[5]][1]#p of a factor
p2<-s[[1]][[5]][2]#p of b factor
c(sa,sb,msa,msb,p,p2)
 }
matrix1 <-t(apply(test,1,fun))
colnames(matrix1)=c("SSA","SSB","MSA","MSB","PA","PB")
write.table(matrix1,file=argv[2],quote=F)




