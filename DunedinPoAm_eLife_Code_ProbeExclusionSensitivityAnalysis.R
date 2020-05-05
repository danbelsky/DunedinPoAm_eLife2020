###-------------------------------------------------------------------------------------------------------------------------------------------------------###
###             testing DunedinPoAm scores calculated by omitting certain probes
###             This analysis is to explore whether two of the probes selected by the Elastic net algorithm are unduly biasing the score
###             These probes are:
###             cg11897887: this probe appears to index a SNP (Beta values in the Dunedin study are either ~0, ~0.5 or ~0.9)
###             cg05575921: this probe is a well-known smoking-related locus
###
###             To calculate the scores minus each of the probes, the DunedinPoAm projector was run twice on a matrix of Beta values first omitting data for
###             cg11897887, then again omiiting data for cg05575921 (N.B. the projector did not include correction for missing data).
###-------------------------------------------------------------------------------------------------------------------------------------------------------###

#### 1. Create data frame of DunedimPoAm scores and phenotypic data:


### Import mPoA scores

### scores required for:

# 1. "Original" p38-trained
# 2. p38-trained minus SNP probe cg11897887
# 3. p38-trained minus AHRR probe cg05575921


library(dplyr)

## Scores in 38 (i.e. "original") DunedinPoAm

Dunedinp38_mpoas <- read.csv("Dunedinp38_mpoa.csv")
Dunedinp38_mpoas <- rename(Dunedinp38_mpoas, DunedinPoAm_Orig = Joeympoa38trained)

### import ammended scores (remove AHRR probe cg05575921 and SNP probe cg11897887 beta values during calcuation)

Dunedin38_mpoasproberemove <- read.csv("Dunedin38_mpoasproberemove.csv")
Dunedin38_mpoasproberemove <- rename(Dunedin38_mpoasproberemove, DunedinPoAm_noSNP = Joeyp38_noSNP, DunedinPoAm_noAHRR = Joeyp38_noAHRR)
Dunedinp38_mpoas <- inner_join(Dunedinp38_mpoas, Dunedin38_mpoasproberemove, by = c("snum", "Age"))
Dunedinp38_mpoas <- filter(Dunedinp38_mpoas, Age == "38")

library(psych)

describe(select(Dunedinp38_mpoas, -snum))

#                     vars   n  mean   sd median trimmed  mad   min   max range skew kurtosis se
# DunedinPoAm_Orig      2 816  0.98 0.09   0.97    0.98 0.08  0.66  1.40  0.74 0.63     0.90  0
# DunedinPoAm_noSNP     3 816  0.98 0.09   0.97    0.98 0.08  0.65  1.39  0.74 0.63     0.90  0
# DunedinPoAm_noAHRR    4 816  1.19 0.07   1.19    1.19 0.06  0.88  1.54  0.66 0.24     1.46  0


###### Import pace of aging data

library(haven)
PGS_check_8aug2019 <- read_sav("PGS check_8aug2019.sav")
p45PoA <- select(PGS_check_8aug2019, snum, ZFacialAge45, POA45s, PackYrLifTm45)

AgingProject_11Nov2015 <- read.csv("AgingProject_11Nov2015.csv")
p38PoA <- select(AgingProject_11Nov2015, snum, BioAgeKD38, PaceOfAging)
p38PoA$snum <- as.numeric(p38PoA$snum)

DunedinSmokingStatusCategoriesRecode20180507 <- read.delim("DunedinSmokingStatusCategoriesRecode20180507.txt")
smoking38 <- select(DunedinSmokingStatusCategoriesRecode20180507, snum, packyears38)
smoking38$snum <- as.numeric(smoking38$snum)

POA <- full_join(p45PoA, p38PoA, by="snum")
POA <- full_join(POA, smoking38, by="snum")

###### Join all together

Dunedinp38_mpoas$snum <- as.numeric(Dunedinp38_mpoas$snum)

Dunedin_POA <- full_join(Dunedinp38_mpoas, POA, by="snum")

###### Run correlations

### A. all datapoints (minus p26 and alternative aging measures)
library(corrplot)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(select(Dunedin_POA, -snum, -Age))
head(p.mat[, 1:5])

m <- cor(select(Dunedin_POA, -snum, -Age), use="pairwise.complete.obs")
corrplot(m, method="circle", type="lower",  
         order="hclust",  addCoef.col = "white", tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.01, insig = "blank", tl.cex = .7, number.cex=0.4, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

PoAcor <- as.data.frame(m)

write.table(PoAcor, file="DunedinPoAm_minusProbes_correlations.txt", sep="/t", col.names=T, row.names=T)

