start <- Sys.time()

library(simstudy)
library(MCTM)
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(misty)
library(MASS)

## MCTM library is removed from R. 
## Formerly available versions can be obtained from the archive.

set.seed(123456)

# Functions 

## Calculation of transition frequencies

transi<- function(x) {
  TransMatrix(as.numeric(x),order= 1,probs= F)
}



## Calculation of residuals

resid_ln<- function(x){
  RS <- rowSums(x)
  CS <- colSums(x) 
  GT <- sum(x)
  res<- loglm(~1+2,data= x,fitted=T) 
  num<- (x - res$fitted)
  den <-  sqrt(res$fitted*((1 - RS / GT) %*% t(1 - CS / GT))) 
  ASR<- num/den
  ASR
}

# Significance of the residuals


zsign<-function(x){
  ifelse(abs(x) > 1.96,1,0) 
}


zsignc<-function(x){
  ifelse(x > 1.96, "Act",ifelse(x < -1.96, "Inh", "Ind"))
}



# Definition of matrices with 8 states


# Matrices with 8 states


## Matrix A

m8A <- matrix(c(1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
              8, 8, byrow = T)

# Matrix B with p = 6/8, 1 significant cell [B,C]


m8B <- matrix(c(1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 6/8,1/28,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
              8, 8, byrow = T)


# Matrix B with p = 6/8, 2 significant cellS [B,C] y [D,C]


m8C <- matrix(c(1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 6/8,1/28,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 6/8,1/28,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
              8, 8, byrow = T)




# Generation of the sequences

## Sequences with 8 states


seq8A50 <- genMarkov(n = 15000, transMat = m8A, 
                     chainLen = 50, wide = TRUE)
seq8B50 <- genMarkov(n = 20000, transMat = m8B, 
                     chainLen = 50, wide = TRUE)
seq8C50 <- genMarkov(n = 20000, transMat = m8C, 
                     chainLen = 50, wide = TRUE)


seq8A100 <- genMarkov(n = 10000, transMat = m8A, 
                      chainLen = 100, wide = TRUE)
seq8B100 <- genMarkov(n = 15000, transMat = m8B, 
                      chainLen = 100, wide = TRUE)
seq8C100 <- genMarkov(n = 15000, transMat = m8C, 
                      chainLen = 100, wide = TRUE)



seq8A150 <- genMarkov(n = 10000, transMat = m8A, 
                      chainLen = 150, wide = TRUE)
seq8B150 <- genMarkov(n = 15000, transMat = m8B, 
                      chainLen = 150, wide = TRUE)
seq8C150 <- genMarkov(n = 15000, transMat = m8C, 
                      chainLen = 150, wide = TRUE)



seq8A50<- seq8A50[,-1]
seq8B50<- seq8B50[,-1]
seq8C50<- seq8C50[,-1]

seq8A100<- seq8A100[,-1]
seq8B100<- seq8B100[,-1]
seq8C100<- seq8C100[,-1]


seq8A150<- seq8A150[,-1]
seq8B150<- seq8B150[,-1]
seq8C150<- seq8C150[,-1]


## names8 

nombres8<- c("AA","AB","AC","AD","AE","AF","AG","AH",
             "BA","BB","BC","BD","BE","BF","BG","BH",
             "CA","CB","CC","CD","CE","CF","CG","CH",
             "DA","DB","DC","DD","DE","DF","DG","DH",
             "EA","EB","EC","ED","EE","EF","EG","EH",
             "FA","FB","FC","FD","FE","FF","FG","FH",
             "GA","GB","GC","GD","GE","GF","GG","GH",
             "HA","HB","HC","HD","HE","HF","HG","HH")


# Generation of the sequences


## Sequences with 8 states

seq8A50 <- genMarkov(n = 15000, transMat = m8A, 
                     chainLen = 50, wide = TRUE)
seq8B50 <- genMarkov(n = 20000, transMat = m8B, 
                     chainLen = 50, wide = TRUE)
seq8C50 <- genMarkov(n = 20000, transMat = m8C, 
                     chainLen = 50, wide = TRUE)


seq8A100 <- genMarkov(n = 10000, transMat = m8A, 
                      chainLen = 100, wide = TRUE)
seq8B100 <- genMarkov(n = 15000, transMat = m8B, 
                      chainLen = 100, wide = TRUE)
seq8C100 <- genMarkov(n = 15000, transMat = m8C, 
                      chainLen = 100, wide = TRUE)



seq8A150 <- genMarkov(n = 10000, transMat = m8A, 
                      chainLen = 150, wide = TRUE)
seq8B150 <- genMarkov(n = 15000, transMat = m8B, 
                      chainLen = 150, wide = TRUE)
seq8C150 <- genMarkov(n = 15000, transMat = m8C, 
                      chainLen = 150, wide = TRUE)



seq8A50<- seq8A50[,-1]
seq8B50<- seq8B50[,-1]
seq8C50<- seq8C50[,-1]

seq8A100<- seq8A100[,-1]
seq8B100<- seq8B100[,-1]
seq8C100<- seq8C100[,-1]


seq8A150<- seq8A150[,-1]
seq8B150<- seq8B150[,-1]
seq8C150<- seq8C150[,-1]



# Residuals of the 8-category condition

## N = 50, matrix A

transiz8A50<-alply(seq8A50,1,transi)
dim8<- ldply(transiz8A50,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

transiz8A50 <- transiz8A50[-falsos]
chisq8A50<- llply(transiz8A50,function(x) {loglm(~1+2, data=x)$pearson})
chi8A50<- llply(chisq8A50, as.vector)
chi28A50<- unlist(chi8A50)
buenos<- which(chi8A50 <= qchisq(.95,49))
length(buenos)

transiz8A50<- transiz8A50[buenos]

resid8A50<- llply(transiz8A50, function(x) resid_ln(x))

resid8A50v<- ldply(resid8A50,as.vector)
buenos<- complete.cases(resid8A50v)
resid8A50v<- resid8A50v[buenos,]
resid8A50v<- resid8A50v[1:5000,]


resid8A50v$c11 <- sapply(X = resid8A50v$V1,zsignc)
resid8A50v$c12 <- sapply(X = resid8A50v$V9,zsignc)
resid8A50v$c13 <- sapply(X = resid8A50v$V17,zsignc)
resid8A50v$c14 <- sapply(X = resid8A50v$V25,zsignc)
resid8A50v$c15 <- sapply(X = resid8A50v$V33,zsignc)
resid8A50v$c16 <- sapply(X = resid8A50v$V41,zsignc)
resid8A50v$c17 <- sapply(X = resid8A50v$V49,zsignc)
resid8A50v$c18 <- sapply(X = resid8A50v$V57,zsignc)


resid8A50v$c21 <- sapply(X = resid8A50v$V2,zsignc)
resid8A50v$c22 <- sapply(X = resid8A50v$V10,zsignc)
resid8A50v$c23 <- sapply(X = resid8A50v$V18,zsignc)
resid8A50v$c24 <- sapply(X = resid8A50v$V26,zsignc)
resid8A50v$c25 <- sapply(X = resid8A50v$V34,zsignc)
resid8A50v$c26 <- sapply(X = resid8A50v$V42,zsignc)
resid8A50v$c27 <- sapply(X = resid8A50v$V50,zsignc)
resid8A50v$c28 <- sapply(X = resid8A50v$V58,zsignc)


resid8A50v$c31 <- sapply(X = resid8A50v$V3,zsignc)
resid8A50v$c32 <- sapply(X = resid8A50v$V11,zsignc)
resid8A50v$c33 <- sapply(X = resid8A50v$V19,zsignc)
resid8A50v$c34 <- sapply(X = resid8A50v$V27,zsignc)
resid8A50v$c35 <- sapply(X = resid8A50v$V35,zsignc)
resid8A50v$c36 <- sapply(X = resid8A50v$V43,zsignc)
resid8A50v$c37 <- sapply(X = resid8A50v$V51,zsignc)
resid8A50v$c38 <- sapply(X = resid8A50v$V59,zsignc)


resid8A50v$c41 <- sapply(X = resid8A50v$V4,zsignc)
resid8A50v$c42 <- sapply(X = resid8A50v$V12,zsignc)
resid8A50v$c43 <- sapply(X = resid8A50v$V20,zsignc)
resid8A50v$c44 <- sapply(X = resid8A50v$V28,zsignc)
resid8A50v$c45 <- sapply(X = resid8A50v$V36,zsignc)
resid8A50v$c46 <- sapply(X = resid8A50v$V44,zsignc)
resid8A50v$c47 <- sapply(X = resid8A50v$V52,zsignc)
resid8A50v$c48 <- sapply(X = resid8A50v$V60,zsignc)


resid8A50v$c51 <- sapply(X = resid8A50v$V5,zsignc)
resid8A50v$c52 <- sapply(X = resid8A50v$V13,zsignc)
resid8A50v$c53 <- sapply(X = resid8A50v$V21,zsignc)
resid8A50v$c54 <- sapply(X = resid8A50v$V29,zsignc)
resid8A50v$c55 <- sapply(X = resid8A50v$V37,zsignc)
resid8A50v$c56 <- sapply(X = resid8A50v$V45,zsignc)
resid8A50v$c57 <- sapply(X = resid8A50v$V53,zsignc)
resid8A50v$c58 <- sapply(X = resid8A50v$V61,zsignc)


resid8A50v$c61 <- sapply(X = resid8A50v$V6,zsignc)
resid8A50v$c62 <- sapply(X = resid8A50v$V14,zsignc)
resid8A50v$c63 <- sapply(X = resid8A50v$V22,zsignc)
resid8A50v$c64 <- sapply(X = resid8A50v$V30,zsignc)
resid8A50v$c65 <- sapply(X = resid8A50v$V38,zsignc)
resid8A50v$c66 <- sapply(X = resid8A50v$V46,zsignc)
resid8A50v$c67 <- sapply(X = resid8A50v$V54,zsignc)
resid8A50v$c68 <- sapply(X = resid8A50v$V62,zsignc)


resid8A50v$c71 <- sapply(X = resid8A50v$V7,zsignc)
resid8A50v$c72 <- sapply(X = resid8A50v$V15,zsignc)
resid8A50v$c73 <- sapply(X = resid8A50v$V23,zsignc)
resid8A50v$c74 <- sapply(X = resid8A50v$V31,zsignc)
resid8A50v$c75 <- sapply(X = resid8A50v$V39,zsignc)
resid8A50v$c76 <- sapply(X = resid8A50v$V47,zsignc)
resid8A50v$c77 <- sapply(X = resid8A50v$V55,zsignc)
resid8A50v$c78 <- sapply(X = resid8A50v$V63,zsignc)


resid8A50v$c81 <- sapply(X = resid8A50v$V8,zsignc)
resid8A50v$c82 <- sapply(X = resid8A50v$V16,zsignc)
resid8A50v$c83 <- sapply(X = resid8A50v$V24,zsignc)
resid8A50v$c84 <- sapply(X = resid8A50v$V32,zsignc)
resid8A50v$c85 <- sapply(X = resid8A50v$V40,zsignc)
resid8A50v$c86 <- sapply(X = resid8A50v$V48,zsignc)
resid8A50v$c87 <- sapply(X = resid8A50v$V56,zsignc)
resid8A50v$c88 <- sapply(X = resid8A50v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8A50<- resid8A50v[,66:129]

datosLong8A50<-pivot_longer(zsig8A50, cols = 1:64,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong8A50$N <- gl(1,320000, labels = 150)
datosLong8A50$CATEGORY<- rep(nombres8,5000)

tabla8A50<-table(datosLong8A50$CATEGORY,datosLong8A50$MEASURE)
tabla8A50p<- as.data.frame(prop.table(tabla8A50,1))
colnames(tabla8A50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8A50p$CATEGORY<- rep(8,192)
tabla8A50p$MATRIX<- rep("A",192)
tabla8A50p$N<- rep(50,192)


## N = 50, matrix B

transiz8B50<-alply(seq8B50,1,transi)
dim8<- ldply(transiz8B50,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

transiz8B50 <- transiz8B50[-falsos]
chisq8B50<- llply(transiz8B50,function(x) {loglm(~1+2, data=x)$pearson})
chi8B50<- llply(chisq8B50, as.vector)
chi28B50<- unlist(chi8B50)
buenos<- which(chi8B50 > qchisq(.95,49))
length(buenos)

transiz8B50<- transiz8B50[buenos]

resid8B50<- llply(transiz8B50, function(x) resid_ln(x))

resid8B50v<- ldply(resid8B50,as.vector)
buenos<- complete.cases(resid8B50v)
resid8B50v<- resid8B50v[buenos,]
resid8B50v<- resid8B50v[1:5000,]

resid8B50v$c11 <- sapply(X = resid8B50v$V1,zsignc)
resid8B50v$c12 <- sapply(X = resid8B50v$V9,zsignc)
resid8B50v$c13 <- sapply(X = resid8B50v$V17,zsignc)
resid8B50v$c14 <- sapply(X = resid8B50v$V25,zsignc)
resid8B50v$c15 <- sapply(X = resid8B50v$V33,zsignc)
resid8B50v$c16 <- sapply(X = resid8B50v$V41,zsignc)
resid8B50v$c17 <- sapply(X = resid8B50v$V49,zsignc)
resid8B50v$c18 <- sapply(X = resid8B50v$V57,zsignc)


resid8B50v$c21 <- sapply(X = resid8B50v$V2,zsignc)
resid8B50v$c22 <- sapply(X = resid8B50v$V10,zsignc)
resid8B50v$c23 <- sapply(X = resid8B50v$V18,zsignc)
resid8B50v$c24 <- sapply(X = resid8B50v$V26,zsignc)
resid8B50v$c25 <- sapply(X = resid8B50v$V34,zsignc)
resid8B50v$c26 <- sapply(X = resid8B50v$V42,zsignc)
resid8B50v$c27 <- sapply(X = resid8B50v$V50,zsignc)
resid8B50v$c28 <- sapply(X = resid8B50v$V58,zsignc)


resid8B50v$c31 <- sapply(X = resid8B50v$V3,zsignc)
resid8B50v$c32 <- sapply(X = resid8B50v$V11,zsignc)
resid8B50v$c33 <- sapply(X = resid8B50v$V19,zsignc)
resid8B50v$c34 <- sapply(X = resid8B50v$V27,zsignc)
resid8B50v$c35 <- sapply(X = resid8B50v$V35,zsignc)
resid8B50v$c36 <- sapply(X = resid8B50v$V43,zsignc)
resid8B50v$c37 <- sapply(X = resid8B50v$V51,zsignc)
resid8B50v$c38 <- sapply(X = resid8B50v$V59,zsignc)

resid8B50v$c41 <- sapply(X = resid8B50v$V4,zsignc)
resid8B50v$c42 <- sapply(X = resid8B50v$V12,zsignc)
resid8B50v$c43 <- sapply(X = resid8B50v$V20,zsignc)
resid8B50v$c44 <- sapply(X = resid8B50v$V28,zsignc)
resid8B50v$c45 <- sapply(X = resid8B50v$V36,zsignc)
resid8B50v$c46 <- sapply(X = resid8B50v$V44,zsignc)
resid8B50v$c47 <- sapply(X = resid8B50v$V52,zsignc)
resid8B50v$c48 <- sapply(X = resid8B50v$V60,zsignc)


resid8B50v$c51 <- sapply(X = resid8B50v$V5,zsignc)
resid8B50v$c52 <- sapply(X = resid8B50v$V13,zsignc)
resid8B50v$c53 <- sapply(X = resid8B50v$V21,zsignc)
resid8B50v$c54 <- sapply(X = resid8B50v$V29,zsignc)
resid8B50v$c55 <- sapply(X = resid8B50v$V37,zsignc)
resid8B50v$c56 <- sapply(X = resid8B50v$V45,zsignc)
resid8B50v$c57 <- sapply(X = resid8B50v$V53,zsignc)
resid8B50v$c58 <- sapply(X = resid8B50v$V61,zsignc)


resid8B50v$c61 <- sapply(X = resid8B50v$V6,zsignc)
resid8B50v$c62 <- sapply(X = resid8B50v$V14,zsignc)
resid8B50v$c63 <- sapply(X = resid8B50v$V22,zsignc)
resid8B50v$c64 <- sapply(X = resid8B50v$V30,zsignc)
resid8B50v$c65 <- sapply(X = resid8B50v$V38,zsignc)
resid8B50v$c66 <- sapply(X = resid8B50v$V46,zsignc)
resid8B50v$c67 <- sapply(X = resid8B50v$V54,zsignc)
resid8B50v$c68 <- sapply(X = resid8B50v$V62,zsignc)


resid8B50v$c71 <- sapply(X = resid8B50v$V7,zsignc)
resid8B50v$c72 <- sapply(X = resid8B50v$V15,zsignc)
resid8B50v$c73 <- sapply(X = resid8B50v$V23,zsignc)
resid8B50v$c74 <- sapply(X = resid8B50v$V31,zsignc)
resid8B50v$c75 <- sapply(X = resid8B50v$V39,zsignc)
resid8B50v$c76 <- sapply(X = resid8B50v$V47,zsignc)
resid8B50v$c77 <- sapply(X = resid8B50v$V55,zsignc)
resid8B50v$c78 <- sapply(X = resid8B50v$V63,zsignc)


resid8B50v$c81 <- sapply(X = resid8B50v$V8,zsignc)
resid8B50v$c82 <- sapply(X = resid8B50v$V16,zsignc)
resid8B50v$c83 <- sapply(X = resid8B50v$V24,zsignc)
resid8B50v$c84 <- sapply(X = resid8B50v$V32,zsignc)
resid8B50v$c85 <- sapply(X = resid8B50v$V40,zsignc)
resid8B50v$c86 <- sapply(X = resid8B50v$V48,zsignc)
resid8B50v$c87 <- sapply(X = resid8B50v$V56,zsignc)
resid8B50v$c88 <- sapply(X = resid8B50v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8B50<- resid8B50v[,66:129]

datosLong8B50<-pivot_longer(zsig8B50, cols = 1:64,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong8B50$N <- gl(1,320000, labels = 50)

datosLong8B50$CATEGORY<- rep(nombres8,5000)

tabla8B50<-table(datosLong8B50$CATEGORY,datosLong8B50$MEASURE)
tabla8B50p<- as.data.frame(prop.table(tabla8B50,1))
colnames(tabla8B50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8B50p$CATEGORY<- rep(8,192)
tabla8B50p$MATRIX<- rep("B",192)
tabla8B50p$N<- rep(50,192)


## N = 50, matrix C

transiz8C50<-alply(seq8C50,1,transi)
dim8<- ldply(transiz8C50,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

transiz8C50 <- transiz8C50[-falsos]
chisq8C50<- llply(transiz8C50,function(x) {loglm(~1+2, data=x)$pearson})
chi8C50<- llply(chisq8C50, as.vector)
chi28C50<- unlist(chi8C50)
buenos<- which(chi8C50 > qchisq(.95,49))
length(buenos)

transiz8C50<- transiz8C50[buenos]
resid8C50<- llply(transiz8C50, function(x) resid_ln(x))


resid8C50v<- ldply(resid8C50,as.vector)
buenos<- complete.cases(resid8C50v)
resid8C50v<- resid8C50v[buenos,]
resid8C50v<- resid8C50v[1:5000,]


resid8C50v$c11 <- sapply(X = resid8C50v$V1,zsignc)
resid8C50v$c12 <- sapply(X = resid8C50v$V9,zsignc)
resid8C50v$c13 <- sapply(X = resid8C50v$V17,zsignc)
resid8C50v$c14 <- sapply(X = resid8C50v$V25,zsignc)
resid8C50v$c15 <- sapply(X = resid8C50v$V33,zsignc)
resid8C50v$c16 <- sapply(X = resid8C50v$V41,zsignc)
resid8C50v$c17 <- sapply(X = resid8C50v$V49,zsignc)
resid8C50v$c18 <- sapply(X = resid8C50v$V57,zsignc)


resid8C50v$c21 <- sapply(X = resid8C50v$V2,zsignc)
resid8C50v$c22 <- sapply(X = resid8C50v$V10,zsignc)
resid8C50v$c23 <- sapply(X = resid8C50v$V18,zsignc)
resid8C50v$c24 <- sapply(X = resid8C50v$V26,zsignc)
resid8C50v$c25 <- sapply(X = resid8C50v$V34,zsignc)
resid8C50v$c26 <- sapply(X = resid8C50v$V42,zsignc)
resid8C50v$c27 <- sapply(X = resid8C50v$V50,zsignc)
resid8C50v$c28 <- sapply(X = resid8C50v$V58,zsignc)


resid8C50v$c31 <- sapply(X = resid8C50v$V3,zsignc)
resid8C50v$c32 <- sapply(X = resid8C50v$V11,zsignc)
resid8C50v$c33 <- sapply(X = resid8C50v$V19,zsignc)
resid8C50v$c34 <- sapply(X = resid8C50v$V27,zsignc)
resid8C50v$c35 <- sapply(X = resid8C50v$V35,zsignc)
resid8C50v$c36 <- sapply(X = resid8C50v$V43,zsignc)
resid8C50v$c37 <- sapply(X = resid8C50v$V51,zsignc)
resid8C50v$c38 <- sapply(X = resid8C50v$V59,zsignc)


resid8C50v$c41 <- sapply(X = resid8C50v$V4,zsignc)
resid8C50v$c42 <- sapply(X = resid8C50v$V12,zsignc)
resid8C50v$c43 <- sapply(X = resid8C50v$V20,zsignc)
resid8C50v$c44 <- sapply(X = resid8C50v$V28,zsignc)
resid8C50v$c45 <- sapply(X = resid8C50v$V36,zsignc)
resid8C50v$c46 <- sapply(X = resid8C50v$V44,zsignc)
resid8C50v$c47 <- sapply(X = resid8C50v$V52,zsignc)
resid8C50v$c48 <- sapply(X = resid8C50v$V60,zsignc)


resid8C50v$c51 <- sapply(X = resid8C50v$V5,zsignc)
resid8C50v$c52 <- sapply(X = resid8C50v$V13,zsignc)
resid8C50v$c53 <- sapply(X = resid8C50v$V21,zsignc)
resid8C50v$c54 <- sapply(X = resid8C50v$V29,zsignc)
resid8C50v$c55 <- sapply(X = resid8C50v$V37,zsignc)
resid8C50v$c56 <- sapply(X = resid8C50v$V45,zsignc)
resid8C50v$c57 <- sapply(X = resid8C50v$V53,zsignc)
resid8C50v$c58 <- sapply(X = resid8C50v$V61,zsignc)


resid8C50v$c61 <- sapply(X = resid8C50v$V6,zsignc)
resid8C50v$c62 <- sapply(X = resid8C50v$V14,zsignc)
resid8C50v$c63 <- sapply(X = resid8C50v$V22,zsignc)
resid8C50v$c64 <- sapply(X = resid8C50v$V30,zsignc)
resid8C50v$c65 <- sapply(X = resid8C50v$V38,zsignc)
resid8C50v$c66 <- sapply(X = resid8C50v$V46,zsignc)
resid8C50v$c67 <- sapply(X = resid8C50v$V54,zsignc)
resid8C50v$c68 <- sapply(X = resid8C50v$V62,zsignc)


resid8C50v$c71 <- sapply(X = resid8C50v$V7,zsignc)
resid8C50v$c72 <- sapply(X = resid8C50v$V15,zsignc)
resid8C50v$c73 <- sapply(X = resid8C50v$V23,zsignc)
resid8C50v$c74 <- sapply(X = resid8C50v$V31,zsignc)
resid8C50v$c75 <- sapply(X = resid8C50v$V39,zsignc)
resid8C50v$c76 <- sapply(X = resid8C50v$V47,zsignc)
resid8C50v$c77 <- sapply(X = resid8C50v$V55,zsignc)
resid8C50v$c78 <- sapply(X = resid8C50v$V63,zsignc)


resid8C50v$c81 <- sapply(X = resid8C50v$V8,zsignc)
resid8C50v$c82 <- sapply(X = resid8C50v$V16,zsignc)
resid8C50v$c83 <- sapply(X = resid8C50v$V24,zsignc)
resid8C50v$c84 <- sapply(X = resid8C50v$V32,zsignc)
resid8C50v$c85 <- sapply(X = resid8C50v$V40,zsignc)
resid8C50v$c86 <- sapply(X = resid8C50v$V48,zsignc)
resid8C50v$c87 <- sapply(X = resid8C50v$V56,zsignc)
resid8C50v$c88 <- sapply(X = resid8C50v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8C50<- resid8C50v[,66:129]

datosLong8C50<-pivot_longer(zsig8C50, cols = 1:64,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong8C50$N <- gl(1,320000, labels = 50)

datosLong8C50$CATEGORY<- rep(nombres8,5000)

tabla8C50<-table(datosLong8C50$CATEGORY,datosLong8C50$MEASURE)
tabla8C50p<- as.data.frame(prop.table(tabla8C50,1))
colnames(tabla8C50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8C50p$CATEGORY<- rep(8,192)
tabla8C50p$MATRIX<- rep("C",192)
tabla8C50p$N<- rep(50,192)



## N = 100, matrix A

transiz8A100<-alply(seq8A100,1,transi)
dim8<- ldply(transiz8A100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

#transiz8A100 <- transiz8A100[-falsos]
chisq8A100<- llply(transiz8A100,function(x) {loglm(~1+2, data=x)$pearson})
chi8A100<- llply(chisq8A100, as.vector)
chi28A100<- unlist(chi8A100)
buenos<- which(chi8A100 <= qchisq(.95,49))
length(buenos)

transiz8A100<- transiz8A100[buenos]

resid8A100<- llply(transiz8A100, function(x) resid_ln(x))

resid8A100v<- ldply(resid8A100,as.vector)
buenos<- complete.cases(resid8A100v)
resid8A100v<- resid8A100v[buenos,]
resid8A100v<- resid8A100v[1:5000,]



resid8A100v$c11 <- sapply(X = resid8A100v$V1,zsignc)
resid8A100v$c12 <- sapply(X = resid8A100v$V9,zsignc)
resid8A100v$c13 <- sapply(X = resid8A100v$V17,zsignc)
resid8A100v$c14 <- sapply(X = resid8A100v$V25,zsignc)
resid8A100v$c15 <- sapply(X = resid8A100v$V33,zsignc)
resid8A100v$c16 <- sapply(X = resid8A100v$V41,zsignc)
resid8A100v$c17 <- sapply(X = resid8A100v$V49,zsignc)
resid8A100v$c18 <- sapply(X = resid8A100v$V57,zsignc)


resid8A100v$c21 <- sapply(X = resid8A100v$V2,zsignc)
resid8A100v$c22 <- sapply(X = resid8A100v$V10,zsignc)
resid8A100v$c23 <- sapply(X = resid8A100v$V18,zsignc)
resid8A100v$c24 <- sapply(X = resid8A100v$V26,zsignc)
resid8A100v$c25 <- sapply(X = resid8A100v$V34,zsignc)
resid8A100v$c26 <- sapply(X = resid8A100v$V42,zsignc)
resid8A100v$c27 <- sapply(X = resid8A100v$V50,zsignc)
resid8A100v$c28 <- sapply(X = resid8A100v$V58,zsignc)


resid8A100v$c31 <- sapply(X = resid8A100v$V3,zsignc)
resid8A100v$c32 <- sapply(X = resid8A100v$V11,zsignc)
resid8A100v$c33 <- sapply(X = resid8A100v$V19,zsignc)
resid8A100v$c34 <- sapply(X = resid8A100v$V27,zsignc)
resid8A100v$c35 <- sapply(X = resid8A100v$V35,zsignc)
resid8A100v$c36 <- sapply(X = resid8A100v$V43,zsignc)
resid8A100v$c37 <- sapply(X = resid8A100v$V51,zsignc)
resid8A100v$c38 <- sapply(X = resid8A100v$V59,zsignc)


resid8A100v$c41 <- sapply(X = resid8A100v$V4,zsignc)
resid8A100v$c42 <- sapply(X = resid8A100v$V12,zsignc)
resid8A100v$c43 <- sapply(X = resid8A100v$V20,zsignc)
resid8A100v$c44 <- sapply(X = resid8A100v$V28,zsignc)
resid8A100v$c45 <- sapply(X = resid8A100v$V36,zsignc)
resid8A100v$c46 <- sapply(X = resid8A100v$V44,zsignc)
resid8A100v$c47 <- sapply(X = resid8A100v$V52,zsignc)
resid8A100v$c48 <- sapply(X = resid8A100v$V60,zsignc)


resid8A100v$c51 <- sapply(X = resid8A100v$V5,zsignc)
resid8A100v$c52 <- sapply(X = resid8A100v$V13,zsignc)
resid8A100v$c53 <- sapply(X = resid8A100v$V21,zsignc)
resid8A100v$c54 <- sapply(X = resid8A100v$V29,zsignc)
resid8A100v$c55 <- sapply(X = resid8A100v$V37,zsignc)
resid8A100v$c56 <- sapply(X = resid8A100v$V45,zsignc)
resid8A100v$c57 <- sapply(X = resid8A100v$V53,zsignc)
resid8A100v$c58 <- sapply(X = resid8A100v$V61,zsignc)


resid8A100v$c61 <- sapply(X = resid8A100v$V6,zsignc)
resid8A100v$c62 <- sapply(X = resid8A100v$V14,zsignc)
resid8A100v$c63 <- sapply(X = resid8A100v$V22,zsignc)
resid8A100v$c64 <- sapply(X = resid8A100v$V30,zsignc)
resid8A100v$c65 <- sapply(X = resid8A100v$V38,zsignc)
resid8A100v$c66 <- sapply(X = resid8A100v$V46,zsignc)
resid8A100v$c67 <- sapply(X = resid8A100v$V54,zsignc)
resid8A100v$c68 <- sapply(X = resid8A100v$V62,zsignc)


resid8A100v$c71 <- sapply(X = resid8A100v$V7,zsignc)
resid8A100v$c72 <- sapply(X = resid8A100v$V15,zsignc)
resid8A100v$c73 <- sapply(X = resid8A100v$V23,zsignc)
resid8A100v$c74 <- sapply(X = resid8A100v$V31,zsignc)
resid8A100v$c75 <- sapply(X = resid8A100v$V39,zsignc)
resid8A100v$c76 <- sapply(X = resid8A100v$V47,zsignc)
resid8A100v$c77 <- sapply(X = resid8A100v$V55,zsignc)
resid8A100v$c78 <- sapply(X = resid8A100v$V63,zsignc)


resid8A100v$c81 <- sapply(X = resid8A100v$V8,zsignc)
resid8A100v$c82 <- sapply(X = resid8A100v$V16,zsignc)
resid8A100v$c83 <- sapply(X = resid8A100v$V24,zsignc)
resid8A100v$c84 <- sapply(X = resid8A100v$V32,zsignc)
resid8A100v$c85 <- sapply(X = resid8A100v$V40,zsignc)
resid8A100v$c86 <- sapply(X = resid8A100v$V48,zsignc)
resid8A100v$c87 <- sapply(X = resid8A100v$V56,zsignc)
resid8A100v$c88 <- sapply(X = resid8A100v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8A100<- resid8A100v[,66:129]

datosLong8A100<-pivot_longer(zsig8A100, cols = 1:64,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong8A100$N <- gl(1,320000, labels = 100)
datosLong8A100$CATEGORY<- rep(nombres8,5000)

tabla8A100<-table(datosLong8A100$CATEGORY,datosLong8A100$MEASURE)
tabla8A100p<- as.data.frame(prop.table(tabla8A100,1))
colnames(tabla8A100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8A100p$CATEGORY<- rep(8,192)
tabla8A100p$MATRIX<- rep("A",192)
tabla8A100p$N<- rep(100,192)

##


## N = 100, matrix B

transiz8B100<-alply(seq8B100,1,transi)

dim8<- ldply(transiz8B100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

transiz8B100 <- transiz8B100[-falsos]
chisq8B100<- llply(transiz8B100,function(x) {loglm(~1+2, data=x)$pearson})
chi8B100<- llply(chisq8B100, as.vector)
chi28B100<- unlist(chi8B100)
buenos<- which(chi8B100 > qchisq(.95,49))
length(buenos)

transiz8B100<- transiz8B100[buenos]
transiz8B100<- transiz8B100[1:5000]
resid8B100<- llply(transiz8B100, function(x) resid_ln(x))


resid8B100v<- ldply(resid8B100,as.vector)
buenos<- complete.cases(resid8B100v)
resid8B100v<- resid8B100v[buenos,]
resid8B100v<- resid8B100v[1:5000,]



resid8B100v$c11 <- sapply(X = resid8B100v$V1,zsignc)
resid8B100v$c12 <- sapply(X = resid8B100v$V9,zsignc)
resid8B100v$c13 <- sapply(X = resid8B100v$V17,zsignc)
resid8B100v$c14 <- sapply(X = resid8B100v$V25,zsignc)
resid8B100v$c15 <- sapply(X = resid8B100v$V33,zsignc)
resid8B100v$c16 <- sapply(X = resid8B100v$V41,zsignc)
resid8B100v$c17 <- sapply(X = resid8B100v$V49,zsignc)
resid8B100v$c18 <- sapply(X = resid8B100v$V57,zsignc)


resid8B100v$c21 <- sapply(X = resid8B100v$V2,zsignc)
resid8B100v$c22 <- sapply(X = resid8B100v$V10,zsignc)
resid8B100v$c23 <- sapply(X = resid8B100v$V18,zsignc)
resid8B100v$c24 <- sapply(X = resid8B100v$V26,zsignc)
resid8B100v$c25 <- sapply(X = resid8B100v$V34,zsignc)
resid8B100v$c26 <- sapply(X = resid8B100v$V42,zsignc)
resid8B100v$c27 <- sapply(X = resid8B100v$V50,zsignc)
resid8B100v$c28 <- sapply(X = resid8B100v$V58,zsignc)


resid8B100v$c31 <- sapply(X = resid8B100v$V3,zsignc)
resid8B100v$c32 <- sapply(X = resid8B100v$V11,zsignc)
resid8B100v$c33 <- sapply(X = resid8B100v$V19,zsignc)
resid8B100v$c34 <- sapply(X = resid8B100v$V27,zsignc)
resid8B100v$c35 <- sapply(X = resid8B100v$V35,zsignc)
resid8B100v$c36 <- sapply(X = resid8B100v$V43,zsignc)
resid8B100v$c37 <- sapply(X = resid8B100v$V51,zsignc)
resid8B100v$c38 <- sapply(X = resid8B100v$V59,zsignc)


resid8B100v$c41 <- sapply(X = resid8B100v$V4,zsignc)
resid8B100v$c42 <- sapply(X = resid8B100v$V12,zsignc)
resid8B100v$c43 <- sapply(X = resid8B100v$V20,zsignc)
resid8B100v$c44 <- sapply(X = resid8B100v$V28,zsignc)
resid8B100v$c45 <- sapply(X = resid8B100v$V36,zsignc)
resid8B100v$c46 <- sapply(X = resid8B100v$V44,zsignc)
resid8B100v$c47 <- sapply(X = resid8B100v$V52,zsignc)
resid8B100v$c48 <- sapply(X = resid8B100v$V60,zsignc)


resid8B100v$c51 <- sapply(X = resid8B100v$V5,zsignc)
resid8B100v$c52 <- sapply(X = resid8B100v$V13,zsignc)
resid8B100v$c53 <- sapply(X = resid8B100v$V21,zsignc)
resid8B100v$c54 <- sapply(X = resid8B100v$V29,zsignc)
resid8B100v$c55 <- sapply(X = resid8B100v$V37,zsignc)
resid8B100v$c56 <- sapply(X = resid8B100v$V45,zsignc)
resid8B100v$c57 <- sapply(X = resid8B100v$V53,zsignc)
resid8B100v$c58 <- sapply(X = resid8B100v$V61,zsignc)


resid8B100v$c61 <- sapply(X = resid8B100v$V6,zsignc)
resid8B100v$c62 <- sapply(X = resid8B100v$V14,zsignc)
resid8B100v$c63 <- sapply(X = resid8B100v$V22,zsignc)
resid8B100v$c64 <- sapply(X = resid8B100v$V30,zsignc)
resid8B100v$c65 <- sapply(X = resid8B100v$V38,zsignc)
resid8B100v$c66 <- sapply(X = resid8B100v$V46,zsignc)
resid8B100v$c67 <- sapply(X = resid8B100v$V54,zsignc)
resid8B100v$c68 <- sapply(X = resid8B100v$V62,zsignc)


resid8B100v$c71 <- sapply(X = resid8B100v$V7,zsignc)
resid8B100v$c72 <- sapply(X = resid8B100v$V15,zsignc)
resid8B100v$c73 <- sapply(X = resid8B100v$V23,zsignc)
resid8B100v$c74 <- sapply(X = resid8B100v$V31,zsignc)
resid8B100v$c75 <- sapply(X = resid8B100v$V39,zsignc)
resid8B100v$c76 <- sapply(X = resid8B100v$V47,zsignc)
resid8B100v$c77 <- sapply(X = resid8B100v$V55,zsignc)
resid8B100v$c78 <- sapply(X = resid8B100v$V63,zsignc)


resid8B100v$c81 <- sapply(X = resid8B100v$V8,zsignc)
resid8B100v$c82 <- sapply(X = resid8B100v$V16,zsignc)
resid8B100v$c83 <- sapply(X = resid8B100v$V24,zsignc)
resid8B100v$c84 <- sapply(X = resid8B100v$V32,zsignc)
resid8B100v$c85 <- sapply(X = resid8B100v$V40,zsignc)
resid8B100v$c86 <- sapply(X = resid8B100v$V48,zsignc)
resid8B100v$c87 <- sapply(X = resid8B100v$V56,zsignc)
resid8B100v$c88 <- sapply(X = resid8B100v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns 

zsig8B100<- resid8B100v[,66:129]

datosLong8B100<-pivot_longer(zsig8B100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong8B100$N <- gl(1,320000, labels = 100)

datosLong8B100$CATEGORY<- rep(nombres8,5000)

tabla8B100<-table(datosLong8B100$CATEGORY,datosLong8B100$Medida)
tabla8B100p<- as.data.frame(prop.table(tabla8B100,1))
colnames(tabla8B100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8B100p$CATEGORY<- rep(8,192)
tabla8B100p$MATRIX<- rep("B",192)
tabla8B100p$N<- rep(100,192)


####


## N = 100, matrix C

transiz8C100<-alply(seq8C100,1,transi)

dim8<- ldply(transiz8C100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

#transiz8C100 <- transiz8C100[-falsos]
chisq8C100<- llply(transiz8C100,function(x) {loglm(~1+2, data=x)$pearson})
chi8C100<- llply(chisq8C100, as.vector)
chi28C100<- unlist(chi8C100)
buenos<- which(chi8C100 > qchisq(.95,49))
length(buenos)

transiz8C100<- transiz8C100[buenos]
transiz8C100<- transiz8C100[1:5000]
resid8C100<- llply(transiz8C100, function(x) resid_ln(x))


resid8C100v<- ldply(resid8C100,as.vector)
buenos<- complete.cases(resid8C100v)
resid8C100v<- resid8C100v[buenos,]
resid8C100v<- resid8C100v[1:5000,]


resid8C100v$c11 <- sapply(X = resid8C100v$V1,zsignc)
resid8C100v$c12 <- sapply(X = resid8C100v$V9,zsignc)
resid8C100v$c13 <- sapply(X = resid8C100v$V17,zsignc)
resid8C100v$c14 <- sapply(X = resid8C100v$V25,zsignc)
resid8C100v$c15 <- sapply(X = resid8C100v$V33,zsignc)
resid8C100v$c16 <- sapply(X = resid8C100v$V41,zsignc)
resid8C100v$c17 <- sapply(X = resid8C100v$V49,zsignc)
resid8C100v$c18 <- sapply(X = resid8C100v$V57,zsignc)


resid8C100v$c21 <- sapply(X = resid8C100v$V2,zsignc)
resid8C100v$c22 <- sapply(X = resid8C100v$V10,zsignc)
resid8C100v$c23 <- sapply(X = resid8C100v$V18,zsignc)
resid8C100v$c24 <- sapply(X = resid8C100v$V26,zsignc)
resid8C100v$c25 <- sapply(X = resid8C100v$V34,zsignc)
resid8C100v$c26 <- sapply(X = resid8C100v$V42,zsignc)
resid8C100v$c27 <- sapply(X = resid8C100v$V50,zsignc)
resid8C100v$c28 <- sapply(X = resid8C100v$V58,zsignc)


resid8C100v$c31 <- sapply(X = resid8C100v$V3,zsignc)
resid8C100v$c32 <- sapply(X = resid8C100v$V11,zsignc)
resid8C100v$c33 <- sapply(X = resid8C100v$V19,zsignc)
resid8C100v$c34 <- sapply(X = resid8C100v$V27,zsignc)
resid8C100v$c35 <- sapply(X = resid8C100v$V35,zsignc)
resid8C100v$c36 <- sapply(X = resid8C100v$V43,zsignc)
resid8C100v$c37 <- sapply(X = resid8C100v$V51,zsignc)
resid8C100v$c38 <- sapply(X = resid8C100v$V59,zsignc)


resid8C100v$c41 <- sapply(X = resid8C100v$V4,zsignc)
resid8C100v$c42 <- sapply(X = resid8C100v$V12,zsignc)
resid8C100v$c43 <- sapply(X = resid8C100v$V20,zsignc)
resid8C100v$c44 <- sapply(X = resid8C100v$V28,zsignc)
resid8C100v$c45 <- sapply(X = resid8C100v$V36,zsignc)
resid8C100v$c46 <- sapply(X = resid8C100v$V44,zsignc)
resid8C100v$c47 <- sapply(X = resid8C100v$V52,zsignc)
resid8C100v$c48 <- sapply(X = resid8C100v$V60,zsignc)


resid8C100v$c51 <- sapply(X = resid8C100v$V5,zsignc)
resid8C100v$c52 <- sapply(X = resid8C100v$V13,zsignc)
resid8C100v$c53 <- sapply(X = resid8C100v$V21,zsignc)
resid8C100v$c54 <- sapply(X = resid8C100v$V29,zsignc)
resid8C100v$c55 <- sapply(X = resid8C100v$V37,zsignc)
resid8C100v$c56 <- sapply(X = resid8C100v$V45,zsignc)
resid8C100v$c57 <- sapply(X = resid8C100v$V53,zsignc)
resid8C100v$c58 <- sapply(X = resid8C100v$V61,zsignc)


resid8C100v$c61 <- sapply(X = resid8C100v$V6,zsignc)
resid8C100v$c62 <- sapply(X = resid8C100v$V14,zsignc)
resid8C100v$c63 <- sapply(X = resid8C100v$V22,zsignc)
resid8C100v$c64 <- sapply(X = resid8C100v$V30,zsignc)
resid8C100v$c65 <- sapply(X = resid8C100v$V38,zsignc)
resid8C100v$c66 <- sapply(X = resid8C100v$V46,zsignc)
resid8C100v$c67 <- sapply(X = resid8C100v$V54,zsignc)
resid8C100v$c68 <- sapply(X = resid8C100v$V62,zsignc)


resid8C100v$c71 <- sapply(X = resid8C100v$V7,zsignc)
resid8C100v$c72 <- sapply(X = resid8C100v$V15,zsignc)
resid8C100v$c73 <- sapply(X = resid8C100v$V23,zsignc)
resid8C100v$c74 <- sapply(X = resid8C100v$V31,zsignc)
resid8C100v$c75 <- sapply(X = resid8C100v$V39,zsignc)
resid8C100v$c76 <- sapply(X = resid8C100v$V47,zsignc)
resid8C100v$c77 <- sapply(X = resid8C100v$V55,zsignc)
resid8C100v$c78 <- sapply(X = resid8C100v$V63,zsignc)


resid8C100v$c81 <- sapply(X = resid8C100v$V8,zsignc)
resid8C100v$c82 <- sapply(X = resid8C100v$V16,zsignc)
resid8C100v$c83 <- sapply(X = resid8C100v$V24,zsignc)
resid8C100v$c84 <- sapply(X = resid8C100v$V32,zsignc)
resid8C100v$c85 <- sapply(X = resid8C100v$V40,zsignc)
resid8C100v$c86 <- sapply(X = resid8C100v$V48,zsignc)
resid8C100v$c87 <- sapply(X = resid8C100v$V56,zsignc)
resid8C100v$c88 <- sapply(X = resid8C100v$V64,zsignc)


# Calculation of activation, inhibition, and independence patterns

zsig8C100<- resid8C100v[,66:129]

datosLong8C100<-pivot_longer(zsig8C100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong8C100$N <- gl(1,320000, labels = 100)

datosLong8C100$CATEGORY<- rep(nombres8,5000)

tabla8C100<-table(datosLong8C100$CATEGORY,datosLong8C100$Medida)
tabla8C100p<- as.data.frame(prop.table(tabla8C100,1))
colnames(tabla8C100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8C100p$CATEGORY<- rep(8,192)
tabla8C100p$MATRIX<- rep("C",192)
tabla8C100p$N<- rep(100,192)


## N = 150, matrix A

transiz8A150<-alply(seq8A150,1,transi)
dim8<- ldply(transiz8A150,dim)
falsos<-which(dim8$V1 != 8| dim8$V2 != 8)
length(falsos)

#transiz8A150 <- transiz8A150[-falsos]
chisq8A150<- llply(transiz8A150,function(x) {loglm(~1+2, data=x)$pearson})
chi8A150<- llply(chisq8A150, as.vector)
chi28A150<- unlist(chi8A150)
buenos<- which(chi8A150 <= qchisq(.95,49))
length(buenos)

transiz8A150<- transiz8A150[buenos]

resid8A150<- llply(transiz8A150, function(x) resid_ln(x))

resid8A150v<- ldply(resid8A150,as.vector)
buenos<- complete.cases(resid8A150v)
resid8A150v<- resid8A150v[buenos,]
resid8A150v<- resid8A150v[1:5000,]


resid8A150v$c11 <- sapply(X = resid8A150v$V1,zsignc)
resid8A150v$c12 <- sapply(X = resid8A150v$V9,zsignc)
resid8A150v$c13 <- sapply(X = resid8A150v$V17,zsignc)
resid8A150v$c14 <- sapply(X = resid8A150v$V25,zsignc)
resid8A150v$c15 <- sapply(X = resid8A150v$V33,zsignc)
resid8A150v$c16 <- sapply(X = resid8A150v$V41,zsignc)
resid8A150v$c17 <- sapply(X = resid8A150v$V49,zsignc)
resid8A150v$c18 <- sapply(X = resid8A150v$V57,zsignc)


resid8A150v$c21 <- sapply(X = resid8A150v$V2,zsignc)
resid8A150v$c22 <- sapply(X = resid8A150v$V10,zsignc)
resid8A150v$c23 <- sapply(X = resid8A150v$V18,zsignc)
resid8A150v$c24 <- sapply(X = resid8A150v$V26,zsignc)
resid8A150v$c25 <- sapply(X = resid8A150v$V34,zsignc)
resid8A150v$c26 <- sapply(X = resid8A150v$V42,zsignc)
resid8A150v$c27 <- sapply(X = resid8A150v$V50,zsignc)
resid8A150v$c28 <- sapply(X = resid8A150v$V58,zsignc)


resid8A150v$c31 <- sapply(X = resid8A150v$V3,zsignc)
resid8A150v$c32 <- sapply(X = resid8A150v$V11,zsignc)
resid8A150v$c33 <- sapply(X = resid8A150v$V19,zsignc)
resid8A150v$c34 <- sapply(X = resid8A150v$V27,zsignc)
resid8A150v$c35 <- sapply(X = resid8A150v$V35,zsignc)
resid8A150v$c36 <- sapply(X = resid8A150v$V43,zsignc)
resid8A150v$c37 <- sapply(X = resid8A150v$V51,zsignc)
resid8A150v$c38 <- sapply(X = resid8A150v$V59,zsignc)


resid8A150v$c41 <- sapply(X = resid8A150v$V4,zsignc)
resid8A150v$c42 <- sapply(X = resid8A150v$V12,zsignc)
resid8A150v$c43 <- sapply(X = resid8A150v$V20,zsignc)
resid8A150v$c44 <- sapply(X = resid8A150v$V28,zsignc)
resid8A150v$c45 <- sapply(X = resid8A150v$V36,zsignc)
resid8A150v$c46 <- sapply(X = resid8A150v$V44,zsignc)
resid8A150v$c47 <- sapply(X = resid8A150v$V52,zsignc)
resid8A150v$c48 <- sapply(X = resid8A150v$V60,zsignc)


resid8A150v$c51 <- sapply(X = resid8A150v$V5,zsignc)
resid8A150v$c52 <- sapply(X = resid8A150v$V13,zsignc)
resid8A150v$c53 <- sapply(X = resid8A150v$V21,zsignc)
resid8A150v$c54 <- sapply(X = resid8A150v$V29,zsignc)
resid8A150v$c55 <- sapply(X = resid8A150v$V37,zsignc)
resid8A150v$c56 <- sapply(X = resid8A150v$V45,zsignc)
resid8A150v$c57 <- sapply(X = resid8A150v$V53,zsignc)
resid8A150v$c58 <- sapply(X = resid8A150v$V61,zsignc)


resid8A150v$c61 <- sapply(X = resid8A150v$V6,zsignc)
resid8A150v$c62 <- sapply(X = resid8A150v$V14,zsignc)
resid8A150v$c63 <- sapply(X = resid8A150v$V22,zsignc)
resid8A150v$c64 <- sapply(X = resid8A150v$V30,zsignc)
resid8A150v$c65 <- sapply(X = resid8A150v$V38,zsignc)
resid8A150v$c66 <- sapply(X = resid8A150v$V46,zsignc)
resid8A150v$c67 <- sapply(X = resid8A150v$V54,zsignc)
resid8A150v$c68 <- sapply(X = resid8A150v$V62,zsignc)


resid8A150v$c71 <- sapply(X = resid8A150v$V7,zsignc)
resid8A150v$c72 <- sapply(X = resid8A150v$V15,zsignc)
resid8A150v$c73 <- sapply(X = resid8A150v$V23,zsignc)
resid8A150v$c74 <- sapply(X = resid8A150v$V31,zsignc)
resid8A150v$c75 <- sapply(X = resid8A150v$V39,zsignc)
resid8A150v$c76 <- sapply(X = resid8A150v$V47,zsignc)
resid8A150v$c77 <- sapply(X = resid8A150v$V55,zsignc)
resid8A150v$c78 <- sapply(X = resid8A150v$V63,zsignc)


resid8A150v$c81 <- sapply(X = resid8A150v$V8,zsignc)
resid8A150v$c82 <- sapply(X = resid8A150v$V16,zsignc)
resid8A150v$c83 <- sapply(X = resid8A150v$V24,zsignc)
resid8A150v$c84 <- sapply(X = resid8A150v$V32,zsignc)
resid8A150v$c85 <- sapply(X = resid8A150v$V40,zsignc)
resid8A150v$c86 <- sapply(X = resid8A150v$V48,zsignc)
resid8A150v$c87 <- sapply(X = resid8A150v$V56,zsignc)
resid8A150v$c88 <- sapply(X = resid8A150v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8A150<- resid8A150v[,66:129]

datosLong8A150<-pivot_longer(zsig8A150, cols = 1:64,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong8A150$N <- gl(1,320000, labels = 150)
datosLong8A150$CATEGORY<- rep(nombres8,5000)

tabla8A150<-table(datosLong8A150$CATEGORY,datosLong8A150$MEASURE)
tabla8A150p<- as.data.frame(prop.table(tabla8A150,1))
colnames(tabla8A150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8A150p$CATEGORY<- rep(8,192)
tabla8A150p$MATRIX<- rep("A",192)
tabla8A150p$N<- rep(150,192)



## N = 150, matrix B

transiz8B150<-alply(seq8B150,1,transi)
dim8<- ldply(transiz8B150,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

#transiz8B150 <- transiz8B150[-falsos]
chisq8B150<- llply(transiz8B150,function(x) {loglm(~1+2, data=x)$pearson})
chi8B150<- llply(chisq8B150, as.vector)
chi28B150<- unlist(chi8B150)
buenos<- which(chi8B150 > qchisq(.95,49))
length(buenos)

transiz8B150<- transiz8B150[buenos]

resid8B150<- llply(transiz8B150, function(x) resid_ln(x))

resid8B150v<- ldply(resid8B150,as.vector)
buenos<- complete.cases(resid8B150v)
resid8B150v<- resid8B150v[buenos,]
resid8B150v<- resid8B150v[1:5000,]

resid8B150v$c11 <- sapply(X = resid8B150v$V1,zsignc)
resid8B150v$c12 <- sapply(X = resid8B150v$V9,zsignc)
resid8B150v$c13 <- sapply(X = resid8B150v$V17,zsignc)
resid8B150v$c14 <- sapply(X = resid8B150v$V25,zsignc)
resid8B150v$c15 <- sapply(X = resid8B150v$V33,zsignc)
resid8B150v$c16 <- sapply(X = resid8B150v$V41,zsignc)
resid8B150v$c17 <- sapply(X = resid8B150v$V49,zsignc)
resid8B150v$c18 <- sapply(X = resid8B150v$V57,zsignc)


resid8B150v$c21 <- sapply(X = resid8B150v$V2,zsignc)
resid8B150v$c22 <- sapply(X = resid8B150v$V10,zsignc)
resid8B150v$c23 <- sapply(X = resid8B150v$V18,zsignc)
resid8B150v$c24 <- sapply(X = resid8B150v$V26,zsignc)
resid8B150v$c25 <- sapply(X = resid8B150v$V34,zsignc)
resid8B150v$c26 <- sapply(X = resid8B150v$V42,zsignc)
resid8B150v$c27 <- sapply(X = resid8B150v$V50,zsignc)
resid8B150v$c28 <- sapply(X = resid8B150v$V58,zsignc)


resid8B150v$c31 <- sapply(X = resid8B150v$V3,zsignc)
resid8B150v$c32 <- sapply(X = resid8B150v$V11,zsignc)
resid8B150v$c33 <- sapply(X = resid8B150v$V19,zsignc)
resid8B150v$c34 <- sapply(X = resid8B150v$V27,zsignc)
resid8B150v$c35 <- sapply(X = resid8B150v$V35,zsignc)
resid8B150v$c36 <- sapply(X = resid8B150v$V43,zsignc)
resid8B150v$c37 <- sapply(X = resid8B150v$V51,zsignc)
resid8B150v$c38 <- sapply(X = resid8B150v$V59,zsignc)


resid8B150v$c41 <- sapply(X = resid8B150v$V4,zsignc)
resid8B150v$c42 <- sapply(X = resid8B150v$V12,zsignc)
resid8B150v$c43 <- sapply(X = resid8B150v$V20,zsignc)
resid8B150v$c44 <- sapply(X = resid8B150v$V28,zsignc)
resid8B150v$c45 <- sapply(X = resid8B150v$V36,zsignc)
resid8B150v$c46 <- sapply(X = resid8B150v$V44,zsignc)
resid8B150v$c47 <- sapply(X = resid8B150v$V52,zsignc)
resid8B150v$c48 <- sapply(X = resid8B150v$V60,zsignc)


resid8B150v$c51 <- sapply(X = resid8B150v$V5,zsignc)
resid8B150v$c52 <- sapply(X = resid8B150v$V13,zsignc)
resid8B150v$c53 <- sapply(X = resid8B150v$V21,zsignc)
resid8B150v$c54 <- sapply(X = resid8B150v$V29,zsignc)
resid8B150v$c55 <- sapply(X = resid8B150v$V37,zsignc)
resid8B150v$c56 <- sapply(X = resid8B150v$V45,zsignc)
resid8B150v$c57 <- sapply(X = resid8B150v$V53,zsignc)
resid8B150v$c58 <- sapply(X = resid8B150v$V61,zsignc)


resid8B150v$c61 <- sapply(X = resid8B150v$V6,zsignc)
resid8B150v$c62 <- sapply(X = resid8B150v$V14,zsignc)
resid8B150v$c63 <- sapply(X = resid8B150v$V22,zsignc)
resid8B150v$c64 <- sapply(X = resid8B150v$V30,zsignc)
resid8B150v$c65 <- sapply(X = resid8B150v$V38,zsignc)
resid8B150v$c66 <- sapply(X = resid8B150v$V46,zsignc)
resid8B150v$c67 <- sapply(X = resid8B150v$V54,zsignc)
resid8B150v$c68 <- sapply(X = resid8B150v$V62,zsignc)


resid8B150v$c71 <- sapply(X = resid8B150v$V7,zsignc)
resid8B150v$c72 <- sapply(X = resid8B150v$V15,zsignc)
resid8B150v$c73 <- sapply(X = resid8B150v$V23,zsignc)
resid8B150v$c74 <- sapply(X = resid8B150v$V31,zsignc)
resid8B150v$c75 <- sapply(X = resid8B150v$V39,zsignc)
resid8B150v$c76 <- sapply(X = resid8B150v$V47,zsignc)
resid8B150v$c77 <- sapply(X = resid8B150v$V55,zsignc)
resid8B150v$c78 <- sapply(X = resid8B150v$V63,zsignc)


resid8B150v$c81 <- sapply(X = resid8B150v$V8,zsignc)
resid8B150v$c82 <- sapply(X = resid8B150v$V16,zsignc)
resid8B150v$c83 <- sapply(X = resid8B150v$V24,zsignc)
resid8B150v$c84 <- sapply(X = resid8B150v$V32,zsignc)
resid8B150v$c85 <- sapply(X = resid8B150v$V40,zsignc)
resid8B150v$c86 <- sapply(X = resid8B150v$V48,zsignc)
resid8B150v$c87 <- sapply(X = resid8B150v$V56,zsignc)
resid8B150v$c88 <- sapply(X = resid8B150v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8B150<- resid8B150v[,66:129]

datosLong8B150<-pivot_longer(zsig8B150, cols = 1:64,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong8B150$N <- gl(1,320000, labels = 150)
datosLong8B150$CATEGORY<- rep(nombres8,5000)

tabla8B150<-table(datosLong8B150$CATEGORY,datosLong8B150$MEASURE)
tabla8B150p<- as.data.frame(prop.table(tabla8B150,1))
colnames(tabla8B150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8B150p$CATEGORY<- rep(8,192)
tabla8B150p$MATRIX<- rep("B",192)
tabla8B150p$N<- rep(150,192)


# N = 150 matrix C

transiz8C150 <-alply(seq8C150,1,transi)

dim8<- ldply(transiz8C150,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

#transiz8C150 <- transiz8C150[-falsos]
chisq8C150<- llply(transiz8C150,function(x) {loglm(~1+2, data=x)$pearson})
chi8C150<- llply(chisq8C150, as.vector)
chi28C150<- unlist(chi8C150)
buenos<- which(chi8C150 > qchisq(.95,25))
length(buenos)

transiz8C150<- transiz8C150[buenos]

resid8C150<- llply(transiz8C150, function(x) resid_ln(x))

resid8C150v<- ldply(resid8C150,as.vector)
buenos<- complete.cases(resid8C150v)
length(buenos)
resid8C150v<- resid8C150v[buenos,]
resid8C150v<- resid8C150v[1:5000,]


resid8C150v$c11 <- sapply(X = resid8C150v$V1,zsignc)
resid8C150v$c12 <- sapply(X = resid8C150v$V9,zsignc)
resid8C150v$c13 <- sapply(X = resid8C150v$V17,zsignc)
resid8C150v$c14 <- sapply(X = resid8C150v$V25,zsignc)
resid8C150v$c15 <- sapply(X = resid8C150v$V33,zsignc)
resid8C150v$c16 <- sapply(X = resid8C150v$V41,zsignc)
resid8C150v$c17 <- sapply(X = resid8C150v$V49,zsignc)
resid8C150v$c18 <- sapply(X = resid8C150v$V57,zsignc)


resid8C150v$c21 <- sapply(X = resid8C150v$V2,zsignc)
resid8C150v$c22 <- sapply(X = resid8C150v$V10,zsignc)
resid8C150v$c23 <- sapply(X = resid8C150v$V18,zsignc)
resid8C150v$c24 <- sapply(X = resid8C150v$V26,zsignc)
resid8C150v$c25 <- sapply(X = resid8C150v$V34,zsignc)
resid8C150v$c26 <- sapply(X = resid8C150v$V42,zsignc)
resid8C150v$c27 <- sapply(X = resid8C150v$V50,zsignc)
resid8C150v$c28 <- sapply(X = resid8C150v$V58,zsignc)


resid8C150v$c31 <- sapply(X = resid8C150v$V3,zsignc)
resid8C150v$c32 <- sapply(X = resid8C150v$V11,zsignc)
resid8C150v$c33 <- sapply(X = resid8C150v$V19,zsignc)
resid8C150v$c34 <- sapply(X = resid8C150v$V27,zsignc)
resid8C150v$c35 <- sapply(X = resid8C150v$V35,zsignc)
resid8C150v$c36 <- sapply(X = resid8C150v$V43,zsignc)
resid8C150v$c37 <- sapply(X = resid8C150v$V51,zsignc)
resid8C150v$c38 <- sapply(X = resid8C150v$V59,zsignc)


resid8C150v$c41 <- sapply(X = resid8C150v$V4,zsignc)
resid8C150v$c42 <- sapply(X = resid8C150v$V12,zsignc)
resid8C150v$c43 <- sapply(X = resid8C150v$V20,zsignc)
resid8C150v$c44 <- sapply(X = resid8C150v$V28,zsignc)
resid8C150v$c45 <- sapply(X = resid8C150v$V36,zsignc)
resid8C150v$c46 <- sapply(X = resid8C150v$V44,zsignc)
resid8C150v$c47 <- sapply(X = resid8C150v$V52,zsignc)
resid8C150v$c48 <- sapply(X = resid8C150v$V60,zsignc)


resid8C150v$c51 <- sapply(X = resid8C150v$V5,zsignc)
resid8C150v$c52 <- sapply(X = resid8C150v$V13,zsignc)
resid8C150v$c53 <- sapply(X = resid8C150v$V21,zsignc)
resid8C150v$c54 <- sapply(X = resid8C150v$V29,zsignc)
resid8C150v$c55 <- sapply(X = resid8C150v$V37,zsignc)
resid8C150v$c56 <- sapply(X = resid8C150v$V45,zsignc)
resid8C150v$c57 <- sapply(X = resid8C150v$V53,zsignc)
resid8C150v$c58 <- sapply(X = resid8C150v$V61,zsignc)


resid8C150v$c61 <- sapply(X = resid8C150v$V6,zsignc)
resid8C150v$c62 <- sapply(X = resid8C150v$V14,zsignc)
resid8C150v$c63 <- sapply(X = resid8C150v$V22,zsignc)
resid8C150v$c64 <- sapply(X = resid8C150v$V30,zsignc)
resid8C150v$c65 <- sapply(X = resid8C150v$V38,zsignc)
resid8C150v$c66 <- sapply(X = resid8C150v$V46,zsignc)
resid8C150v$c67 <- sapply(X = resid8C150v$V54,zsignc)
resid8C150v$c68 <- sapply(X = resid8C150v$V62,zsignc)


resid8C150v$c71 <- sapply(X = resid8C150v$V7,zsignc)
resid8C150v$c72 <- sapply(X = resid8C150v$V15,zsignc)
resid8C150v$c73 <- sapply(X = resid8C150v$V23,zsignc)
resid8C150v$c74 <- sapply(X = resid8C150v$V31,zsignc)
resid8C150v$c75 <- sapply(X = resid8C150v$V39,zsignc)
resid8C150v$c76 <- sapply(X = resid8C150v$V47,zsignc)
resid8C150v$c77 <- sapply(X = resid8C150v$V55,zsignc)
resid8C150v$c78 <- sapply(X = resid8C150v$V63,zsignc)


resid8C150v$c81 <- sapply(X = resid8C150v$V8,zsignc)
resid8C150v$c82 <- sapply(X = resid8C150v$V16,zsignc)
resid8C150v$c83 <- sapply(X = resid8C150v$V24,zsignc)
resid8C150v$c84 <- sapply(X = resid8C150v$V32,zsignc)
resid8C150v$c85 <- sapply(X = resid8C150v$V40,zsignc)
resid8C150v$c86 <- sapply(X = resid8C150v$V48,zsignc)
resid8C150v$c87 <- sapply(X = resid8C150v$V56,zsignc)
resid8C150v$c88 <- sapply(X = resid8C150v$V64,zsignc)



# Calculation of activation, inhibition, and independence patterns

zsig8C150<- resid8C150v[,66:129]

datosLong8C150<-pivot_longer(zsig8C150, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong8C150$N <- gl(1,320000, labels = 150)
datosLong8C150$CATEGORY<- rep(nombres8,5000)

tabla8C150<-table(datosLong8C150$CATEGORY,datosLong8C150$Medida)
tabla8C150p<- as.data.frame(prop.table(tabla8C150,1))
colnames(tabla8C150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla8C150p$CATEGORY<- rep(8,192)
tabla8C150p$MATRIX<- rep("C",192)
tabla8C150p$N<- rep(150,192)


# Data summary table


tabla8_lm<- rbind(tabla8A50p,tabla8A100p,tabla8A150p,
                  tabla8B50p,tabla8B100p,tabla8B150p, 
                  tabla8C50p,tabla8C100p,tabla8C150p)



# Activation patterns

tab8bc1<- subset(tabla8_lm,PATTERN == "BC" & DEPENDENCY == "Act" )
tab8dc1<- subset(tabla8_lm,PATTERN == "DC" & DEPENDENCY == "Act" )
tab8_lm<- rbind(tab8bc1,tab8dc1)
tab8_lm$MODEL<- gl(1,18, labels = "LOGL")
tab8_lm$N<- rep(gl(3,1, labels = c("50","100","150")),6)



# Inhibition patterns


tab8ba1<- subset(tabla8_lm,PATTERN == "BA" & DEPENDENCY == "Inh" )
tab8bb1<- subset(tabla8_lm,PATTERN == "BB" & DEPENDENCY == "Inh" )
tab8bd1<- subset(tabla8_lm,PATTERN == "BD" & DEPENDENCY == "Inh" )
tab8be1<- subset(tabla8_lm,PATTERN == "BE" & DEPENDENCY == "Inh" )
tab8bf1<- subset(tabla8_lm,PATTERN == "BF" & DEPENDENCY == "Inh" )
tab8bg1<- subset(tabla8_lm,PATTERN == "BG" & DEPENDENCY == "Inh" )
tab8bh1<- subset(tabla8_lm,PATTERN == "BH" & DEPENDENCY == "Inh" )


tab8da1<- subset(tabla8_lm,PATTERN == "DA" & DEPENDENCY == "Inh" )
tab8db1<- subset(tabla8_lm,PATTERN == "DB" & DEPENDENCY == "Inh" )
tab8dd1<- subset(tabla8_lm,PATTERN == "DD" & DEPENDENCY == "Inh" )
tab8de1<- subset(tabla8_lm,PATTERN == "DE" & DEPENDENCY == "Inh" )
tab8df1<- subset(tabla8_lm,PATTERN == "DF" & DEPENDENCY == "Inh" )
tab8dg1<- subset(tabla8_lm,PATTERN == "DG" & DEPENDENCY == "Inh" )
tab8dh1<- subset(tabla8_lm,PATTERN == "DH" & DEPENDENCY == "Inh" )

tab8_inh_lm<- rbind(tab8ba1, tab8bb1,tab8bd1,tab8be1,tab8bf1,
                     tab8bg1, tab8bh1,tab8da1,tab8db1,tab8dd1,
                     tab8de1,tab8df1,tab8dg1,tab8dh1)

tab8_inh_lm$MODEL<- gl(1,126, labels = "LOGL")
tab8_inh_lm$N<- rep(gl(3,1, labels = c("50","100","150")),42)

# Spureous pattern

tab8ac1<- subset(tabla8_lm,PATTERN == "AC" & DEPENDENCY == "Inh" )
tab8cc1<- subset(tabla8_lm,PATTERN == "CC" & DEPENDENCY == "Inh" )

tab8_sp_lm<- rbind(tab8ac1,tab8cc1)

tab8_sp_lm$MODEL<- gl(1,18, labels = "LOGL")
tab8_sp_lm$N<- rep(gl(3,1, labels = c("50","100","150")),6)



end<- Sys.time()
print(end-start)
