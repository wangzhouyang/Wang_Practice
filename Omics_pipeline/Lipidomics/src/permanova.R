args <- commandArgs(T)
if (length(args) != 4) {stop("Rscript permanova.R [input1] [input2] [input3] [prefix]
input1: Distance Matrix
input2: Phenotype
input3: no. of levels no more than input3 change as factor, such as: 0, 1, ...
prefix: eg: Bray Euclid
make sure sample of phe.pr and dist.pr have the same rownames")
}

dist.dir <- args[1]
phe.dir <- args[2]
cutoff <- as.numeric(args[3])
prefix <- args[4]
out.dir <- paste(prefix, "_adonis.txt", sep = "") 

library(vegan)

dist.pr <- read.table(dist.dir)
phe.pr <- read.table(phe.dir)
phe.cn <- colnames(phe.pr)
phe.nc <- ncol(phe.pr)

dist.rn <- rownames(dist.pr)
phe.rn <- rownames(phe.pr)

rn <- intersect(phe.rn, dist.rn)

flag <- pmatch(rn, dist.rn)
dist.pr <- dist.pr[flag, flag, drop = F]
flag <- pmatch(rn, phe.rn)
phe.pr <- phe.pr[flag, , drop = F]

out <- matrix(nrow=phe.nc,ncol=9)
colnames(out) <- list("phenotype", "SampleNum", "Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2", "Pr(>F)", "FDR")

for (i in 1:phe.nc) {
  phe <- phe.pr[, i]
  flag <- which(!is.na(phe))
  len <- length(flag)
  phe <- phe[flag]

  if (len == 0 | length(unique(phe)) == 1) {
     out <- list(phe.cn[i], len, NA, NA, NA, NA, NA, NA)
     write.table(out, out.dir, quote = F, sep = "\t", col.names = F, row.names = F, append = T)
     next
   }

  set.seed(0)
  dist <- as.dist(dist.pr[flag, flag, drop = F])
  if (length(unique(phe)) <= cutoff) {
    phe <- as.factor(phe)
  }
  res <- adonis(dist ~ phe, permutations = 9999)
  tmp <- res$aov.tab[1, ]
  out[i,1:8] <- c(phe.cn[i], len, as.character(tmp))
}

p <- out[,8]
fdr <- p.adjust(p, method = "BH")
out[,9] <- fdr

write.table(out, out.dir, quote = F, sep = "\t", col.names = T, row.names = F, append = F)
