#############adjust-age-bmi-drug##########
args=commandArgs(T)
	if(length(args)!=3){stop("Rscript  Logistic.R [input1]  [input2] [prefix]")
		input1:  metabolite_profile
			input2:  phenotype[Diagnosis_Gender_Age_BMI_drug]
			prefix:  output
	}
metabo <- args[1]
Phe    <- args[2]
prefix <- args[3]
metabo_dat   <- t(read.table(metabo,header=T,row.names=1,sep="\t",check.names = F))
phe_dat      <- read.table(Phe,header=T,row.names=1,sep="\t",check.names=F)

#########function####################
logis <-function(x,phe){
	dat  <- data.frame(x,phe)
		fit  <- glm(Diagnosis~.,family=binomial(link='logit'),data=dat)
		coef <- summary(fit)$coefficients
		out  <- gather(data.frame(t(coef[2:5,1:4])))$value
		return(out)               
}

result <- t(apply(metabo_dat,2,logis,phe=phe_dat))
colnames(result) <- as.vector(sapply(rownames(phe)[2:5],function(set) {out<-c(set[1],paste0("pvalue_",set[1]))}))
out_pre<- paste0(prefix,".logist.sta")
write.table(result,out_pre,sep="\t",quote=F,col.name=NA)
