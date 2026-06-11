inicio <- Sys.time()

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
  ifelse(abs(x) > 3.197,1,0) 
}


zsignc<-function(x){
  ifelse(x > 3.197, "Act",ifelse(x < -3.197, "Inh", "Ind"))
}


# Matrices with 6 states


## Matrix   A

m6A <-matrix(c(1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6), 
             6, 6, 
             byrow = T)

## Matrix B p = .75 states = 4, 1 cell significant BC


m6B <-matrix(c(1/6,  1/6,  1/6,  1/6,   1/6,  1/6,
               .05,  .05,  .75,  .05,   .05, .05,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6), 
             6, 6, 
             byrow = T)


## Matrix B p = .75 states = 4, 2 cells significant BC y DC

m6C <-matrix(c(1/6,  1/6,  1/6,  1/6,   1/6,  1/6,
               .05,  .05,  .75,  .05,   .05, .05,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               .05,  .05,  .75,  .05,   .05, .05,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6), 
             6, 6, 
             byrow = T)

## names

nombres6<- c("AA","AB","AC","AD","AE","AF",
            "BA","BB","BC","BD","BE","BF",
            "CA","CB","CC","CD","CE","CF",
            "DA","DB","DC","DD","DE","DF",
            "EA","EB","EC","ED","EE","EF",
            "FA","FB","FC","FD","FE","FF")



# Generation of the sequences

## Sequences with 6 states

seq6A50 <- genMarkov(n = 15000, transMat = m6A, 
                     chainLen = 50, wide = TRUE)
seq6B50 <- genMarkov(n = 15000, transMat = m6B, 
                     chainLen = 50, wide = TRUE)
seq6C50 <- genMarkov(n = 15000, transMat = m6C, 
                     chainLen = 50, wide = TRUE)


seq6A100 <- genMarkov(n = 15000, transMat = m6A, 
                      chainLen = 100, wide = TRUE)
seq6B100 <- genMarkov(n = 15000, transMat = m6B, 
                      chainLen = 100, wide = TRUE)
seq6C100 <- genMarkov(n = 15000, transMat = m6C, 
                      chainLen = 100, wide = TRUE)

seq6A150 <- genMarkov(n = 10000, transMat = m6A,
                      chainLen = 150, wide = TRUE)
seq6B150 <- genMarkov(n = 10000, transMat = m6B,
                      chainLen = 150, wide = TRUE)
seq6C150 <- genMarkov(n = 10000, transMat = m6C,
                      chainLen = 150, wide = TRUE)



seq6A50<- seq6A50[,-1]
seq6B50<- seq6B50[,-1]
seq6C50<- seq6C50[,-1]

seq6A100<- seq6A100[,-1]
seq6B100<- seq6B100[,-1]
seq6C100<- seq6C100[,-1]


seq6A150<- seq6A150[,-1]
seq6B150<- seq6B150[,-1]
seq6C150<- seq6C150[,-1]



# Residuals of the 6-category condition



# Residuals of the 6-category condition

## N = 50, matrix A

transiz6A50<-alply(seq6A50,1,transi)
dim6<- ldply(transiz6A50,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

transiz6A50 <- transiz6A50[-falsos]
chisq6A50<- llply(transiz6A50,function(x) {loglm(~1+2, data=x)$pearson})
chi6A50<- llply(chisq6A50, as.vector)
chi26A50<- unlist(chi6A50)
buenos<- which(chi6A50 <= qchisq(.95,25))
length(buenos)

transiz6A50<- transiz6A50[buenos]

resid6A50<- llply(transiz6A50, function(x) resid_ln(x))

resid6A50v<- ldply(resid6A50,as.vector)
buenos<- complete.cases(resid6A50v)
resid6A50v<- resid6A50v[buenos,]
resid6A50v<- resid6A50v[1:5000,]


resid6A50v$c11 <- sapply(X = resid6A50v$V1,zsignc)
resid6A50v$c12 <- sapply(X = resid6A50v$V7,zsignc)
resid6A50v$c13 <- sapply(X = resid6A50v$V13,zsignc)
resid6A50v$c14 <- sapply(X = resid6A50v$V19,zsignc)
resid6A50v$c15 <- sapply(X = resid6A50v$V25,zsignc)
resid6A50v$c16 <- sapply(X = resid6A50v$V31,zsignc)


resid6A50v$c21 <- sapply(X = resid6A50v$V2,zsignc)
resid6A50v$c22 <- sapply(X = resid6A50v$V8,zsignc)
resid6A50v$c23 <- sapply(X = resid6A50v$V14,zsignc)
resid6A50v$c24 <- sapply(X = resid6A50v$V20,zsignc)
resid6A50v$c25 <- sapply(X = resid6A50v$V26,zsignc)
resid6A50v$c26 <- sapply(X = resid6A50v$V28,zsignc)


resid6A50v$c31 <- sapply(X = resid6A50v$V3,zsignc)
resid6A50v$c32 <- sapply(X = resid6A50v$V9,zsignc)
resid6A50v$c33 <- sapply(X = resid6A50v$V15,zsignc)
resid6A50v$c34 <- sapply(X = resid6A50v$V21,zsignc)
resid6A50v$c35 <- sapply(X = resid6A50v$V27,zsignc)
resid6A50v$c36 <- sapply(X = resid6A50v$V33,zsignc)


resid6A50v$c41 <- sapply(X = resid6A50v$V4,zsignc)
resid6A50v$c42 <- sapply(X = resid6A50v$V10,zsignc)
resid6A50v$c43 <- sapply(X = resid6A50v$V16,zsignc)
resid6A50v$c44 <- sapply(X = resid6A50v$V22,zsignc)
resid6A50v$c45 <- sapply(X = resid6A50v$V28,zsignc)
resid6A50v$c46 <- sapply(X = resid6A50v$V34,zsignc)


resid6A50v$c51 <- sapply(X = resid6A50v$V5,zsignc)
resid6A50v$c52 <- sapply(X = resid6A50v$V11,zsignc)
resid6A50v$c53 <- sapply(X = resid6A50v$V17,zsignc)
resid6A50v$c54 <- sapply(X = resid6A50v$V23,zsignc)
resid6A50v$c55 <- sapply(X = resid6A50v$V29,zsignc)
resid6A50v$c56 <- sapply(X = resid6A50v$V32,zsignc)



resid6A50v$c61 <- sapply(X = resid6A50v$V6,zsignc)
resid6A50v$c62 <- sapply(X = resid6A50v$V12,zsignc)
resid6A50v$c63 <- sapply(X = resid6A50v$V18,zsignc)
resid6A50v$c64 <- sapply(X = resid6A50v$V24,zsignc)
resid6A50v$c65 <- sapply(X = resid6A50v$V30,zsignc)
resid6A50v$c66 <- sapply(X = resid6A50v$V36,zsignc)



zsig6A50<- resid6A50v[,38:73]

datosLong6A50<-pivot_longer(zsig6A50, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6A50$N <- gl(1,180000, labels = "50")

datosLong6A50$CATEGORY<- rep(nombres6,5000)

tabla6A50<-table(datosLong6A50$CATEGORY,datosLong6A50$MEASURE)
tabla6A50p<- as.data.frame(prop.table(tabla6A50,1))
colnames(tabla6A50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6A50p$CATEGORY<- rep(6,72)
tabla6A50p$MATRIX<- rep("A",72)
tabla6A50p$N<- rep(50,72)


## N = 50, matrix B

transiz6B50<-alply(seq6B50,1,transi)
dim6<- ldply(transiz6B50,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

transiz6B50 <- transiz6B50[-falsos]
chisq6B50<- llply(transiz6B50,function(x) {loglm(~1+2, data=x)$pearson})
chi6B50<- llply(chisq6B50, as.vector)
chi26B50<- unlist(chi6B50)
buenos<- which(chi6B50 > qchisq(.95,25))
length(buenos)

transiz6B50<- transiz6B50[buenos]

resid6B50<- llply(transiz6B50, function(x) resid_ln(x))

resid6B50v<- ldply(resid6B50,as.vector)
buenos<- complete.cases(resid6B50v)
resid6B50v<- resid6B50v[buenos,]
resid6B50v<- resid6B50v[1:5000,]

resid6B50v$c11 <- sapply(X = resid6B50v$V1,zsignc)
resid6B50v$c12 <- sapply(X = resid6B50v$V7,zsignc)
resid6B50v$c13 <- sapply(X = resid6B50v$V13,zsignc)
resid6B50v$c14 <- sapply(X = resid6B50v$V19,zsignc)
resid6B50v$c15 <- sapply(X = resid6B50v$V25,zsignc)
resid6B50v$c16 <- sapply(X = resid6B50v$V31,zsignc)


resid6B50v$c21 <- sapply(X = resid6B50v$V2,zsignc)
resid6B50v$c22 <- sapply(X = resid6B50v$V8,zsignc)
resid6B50v$c23 <- sapply(X = resid6B50v$V14,zsignc)
resid6B50v$c24 <- sapply(X = resid6B50v$V20,zsignc)
resid6B50v$c25 <- sapply(X = resid6B50v$V26,zsignc)
resid6B50v$c26 <- sapply(X = resid6B50v$V28,zsignc)


resid6B50v$c31 <- sapply(X = resid6B50v$V3,zsignc)
resid6B50v$c32 <- sapply(X = resid6B50v$V9,zsignc)
resid6B50v$c33 <- sapply(X = resid6B50v$V15,zsignc)
resid6B50v$c34 <- sapply(X = resid6B50v$V21,zsignc)
resid6B50v$c35 <- sapply(X = resid6B50v$V27,zsignc)
resid6B50v$c36 <- sapply(X = resid6B50v$V33,zsignc)


resid6B50v$c41 <- sapply(X = resid6B50v$V4,zsignc)
resid6B50v$c42 <- sapply(X = resid6B50v$V10,zsignc)
resid6B50v$c43 <- sapply(X = resid6B50v$V16,zsignc)
resid6B50v$c44 <- sapply(X = resid6B50v$V22,zsignc)
resid6B50v$c45 <- sapply(X = resid6B50v$V28,zsignc)
resid6B50v$c46 <- sapply(X = resid6B50v$V34,zsignc)


resid6B50v$c51 <- sapply(X = resid6B50v$V5,zsignc)
resid6B50v$c52 <- sapply(X = resid6B50v$V11,zsignc)
resid6B50v$c53 <- sapply(X = resid6B50v$V17,zsignc)
resid6B50v$c54 <- sapply(X = resid6B50v$V23,zsignc)
resid6B50v$c55 <- sapply(X = resid6B50v$V29,zsignc)
resid6B50v$c56 <- sapply(X = resid6B50v$V32,zsignc)



resid6B50v$c61 <- sapply(X = resid6B50v$V6,zsignc)
resid6B50v$c62 <- sapply(X = resid6B50v$V12,zsignc)
resid6B50v$c63 <- sapply(X = resid6B50v$V18,zsignc)
resid6B50v$c64 <- sapply(X = resid6B50v$V24,zsignc)
resid6B50v$c65 <- sapply(X = resid6B50v$V30,zsignc)
resid6B50v$c66 <- sapply(X = resid6B50v$V36,zsignc)



zsig6B50<- resid6B50v[,38:73]

datosLong6B50<-pivot_longer(zsig6B50, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6B50$N <- gl(1,180000, labels = "50")

datosLong6B50$CATEGORY<- rep(nombres6,5000)

tabla6B50<-table(datosLong6B50$CATEGORY,datosLong6B50$MEASURE)
tabla6B50p<- as.data.frame(prop.table(tabla6B50,1))
colnames(tabla6B50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6B50p$CATEGORY<- rep(6,108)
tabla6B50p$MATRIX<- rep("B",108)
tabla6B50p$N<- rep(50,108)



## N = 50, matrix C

transiz6C50<-alply(seq6C50,1,transi)
dim6<- ldply(transiz6C50,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

transiz6C50 <- transiz6C50[-falsos]
chisq6C50<- llply(transiz6C50,function(x) {loglm(~1+2, data=x)$pearson})
chi6C50<- llply(chisq6C50, as.vector)
chi26C50<- unlist(chi6C50)
buenos<- which(chi6C50 > qchisq(.95,25))
length(buenos)

transiz6C50<- transiz6C50[buenos]
resid6C50<- llply(transiz6C50, function(x) resid_ln(x))


resid6C50v<- ldply(resid6C50,as.vector)
buenos<- complete.cases(resid6C50v)
resid6C50v<- resid6C50v[buenos,]
resid6C50v<- resid6C50v[1:5000,]


resid6C50v$c11 <- sapply(X = resid6C50v$V1,zsignc)
resid6C50v$c12 <- sapply(X = resid6C50v$V7,zsignc)
resid6C50v$c13 <- sapply(X = resid6C50v$V13,zsignc)
resid6C50v$c14 <- sapply(X = resid6C50v$V19,zsignc)
resid6C50v$c15 <- sapply(X = resid6C50v$V25,zsignc)
resid6C50v$c16 <- sapply(X = resid6C50v$V31,zsignc)


resid6C50v$c21 <- sapply(X = resid6C50v$V2,zsignc)
resid6C50v$c22 <- sapply(X = resid6C50v$V8,zsignc)
resid6C50v$c23 <- sapply(X = resid6C50v$V14,zsignc)
resid6C50v$c24 <- sapply(X = resid6C50v$V20,zsignc)
resid6C50v$c25 <- sapply(X = resid6C50v$V26,zsignc)
resid6C50v$c26 <- sapply(X = resid6C50v$V28,zsignc)


resid6C50v$c31 <- sapply(X = resid6C50v$V3,zsignc)
resid6C50v$c32 <- sapply(X = resid6C50v$V9,zsignc)
resid6C50v$c33 <- sapply(X = resid6C50v$V15,zsignc)
resid6C50v$c34 <- sapply(X = resid6C50v$V21,zsignc)
resid6C50v$c35 <- sapply(X = resid6C50v$V27,zsignc)
resid6C50v$c36 <- sapply(X = resid6C50v$V33,zsignc)


resid6C50v$c41 <- sapply(X = resid6C50v$V4,zsignc)
resid6C50v$c42 <- sapply(X = resid6C50v$V10,zsignc)
resid6C50v$c43 <- sapply(X = resid6C50v$V16,zsignc)
resid6C50v$c44 <- sapply(X = resid6C50v$V22,zsignc)
resid6C50v$c45 <- sapply(X = resid6C50v$V28,zsignc)
resid6C50v$c46 <- sapply(X = resid6C50v$V34,zsignc)


resid6C50v$c51 <- sapply(X = resid6C50v$V5,zsignc)
resid6C50v$c52 <- sapply(X = resid6C50v$V11,zsignc)
resid6C50v$c53 <- sapply(X = resid6C50v$V17,zsignc)
resid6C50v$c54 <- sapply(X = resid6C50v$V23,zsignc)
resid6C50v$c55 <- sapply(X = resid6C50v$V29,zsignc)
resid6C50v$c56 <- sapply(X = resid6C50v$V32,zsignc)



resid6C50v$c61 <- sapply(X = resid6C50v$V6,zsignc)
resid6C50v$c62 <- sapply(X = resid6C50v$V12,zsignc)
resid6C50v$c63 <- sapply(X = resid6C50v$V18,zsignc)
resid6C50v$c64 <- sapply(X = resid6C50v$V24,zsignc)
resid6C50v$c65 <- sapply(X = resid6C50v$V30,zsignc)
resid6C50v$c66 <- sapply(X = resid6C50v$V36,zsignc)



zsig6C50<- resid6C50v[,38:73]

datosLong6C50<-pivot_longer(zsig6C50, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6C50$N <- gl(1,180000, labels = "50")

datosLong6C50$CATEGORY<- rep(nombres6,5000)

tabla6C50<-table(datosLong6C50$CATEGORY,datosLong6C50$MEASURE)
tabla6C50p<- as.data.frame(prop.table(tabla6C50,1))
colnames(tabla6C50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6C50p$CATEGORY<- rep(6,108)
tabla6C50p$MATRIX<- rep("C",108)
tabla6C50p$N<- rep(50,108)



########

## N = 100, matrix A

transiz6A100<-alply(seq6A100,1,transi)
dim6<- ldply(transiz6A100,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

#transiz6A100 <- transiz6A100[-falsos]
chisq6A100<- llply(transiz6A100,function(x) {loglm(~1+2, data=x)$pearson})
chi6A100<- llply(chisq6A100, as.vector)
chi26A100<- unlist(chi6A100)
buenos<- which(chi6A100 <= qchisq(.95,25))
length(buenos)

transiz6A100<- transiz6A100[buenos]

resid6A100<- llply(transiz6A100, function(x) resid_ln(x))

resid6A100v<- ldply(resid6A100,as.vector)
buenos<- complete.cases(resid6A100v)
resid6A100v<- resid6A100v[buenos,]
resid6A100v<- resid6A100v[1:5000,]


resid6A100v$c11 <- sapply(X = resid6A100v$V1,zsignc)
resid6A100v$c12 <- sapply(X = resid6A100v$V7,zsignc)
resid6A100v$c13 <- sapply(X = resid6A100v$V13,zsignc)
resid6A100v$c14 <- sapply(X = resid6A100v$V19,zsignc)
resid6A100v$c15 <- sapply(X = resid6A100v$V25,zsignc)
resid6A100v$c16 <- sapply(X = resid6A100v$V31,zsignc)


resid6A100v$c21 <- sapply(X = resid6A100v$V2,zsignc)
resid6A100v$c22 <- sapply(X = resid6A100v$V8,zsignc)
resid6A100v$c23 <- sapply(X = resid6A100v$V14,zsignc)
resid6A100v$c24 <- sapply(X = resid6A100v$V20,zsignc)
resid6A100v$c25 <- sapply(X = resid6A100v$V26,zsignc)
resid6A100v$c26 <- sapply(X = resid6A100v$V28,zsignc)


resid6A100v$c31 <- sapply(X = resid6A100v$V3,zsignc)
resid6A100v$c32 <- sapply(X = resid6A100v$V9,zsignc)
resid6A100v$c33 <- sapply(X = resid6A100v$V15,zsignc)
resid6A100v$c34 <- sapply(X = resid6A100v$V21,zsignc)
resid6A100v$c35 <- sapply(X = resid6A100v$V27,zsignc)
resid6A100v$c36 <- sapply(X = resid6A100v$V33,zsignc)


resid6A100v$c41 <- sapply(X = resid6A100v$V4,zsignc)
resid6A100v$c42 <- sapply(X = resid6A100v$V10,zsignc)
resid6A100v$c43 <- sapply(X = resid6A100v$V16,zsignc)
resid6A100v$c44 <- sapply(X = resid6A100v$V22,zsignc)
resid6A100v$c45 <- sapply(X = resid6A100v$V28,zsignc)
resid6A100v$c46 <- sapply(X = resid6A100v$V34,zsignc)


resid6A100v$c51 <- sapply(X = resid6A100v$V5,zsignc)
resid6A100v$c52 <- sapply(X = resid6A100v$V11,zsignc)
resid6A100v$c53 <- sapply(X = resid6A100v$V17,zsignc)
resid6A100v$c54 <- sapply(X = resid6A100v$V23,zsignc)
resid6A100v$c55 <- sapply(X = resid6A100v$V29,zsignc)
resid6A100v$c56 <- sapply(X = resid6A100v$V32,zsignc)



resid6A100v$c61 <- sapply(X = resid6A100v$V6,zsignc)
resid6A100v$c62 <- sapply(X = resid6A100v$V12,zsignc)
resid6A100v$c63 <- sapply(X = resid6A100v$V18,zsignc)
resid6A100v$c64 <- sapply(X = resid6A100v$V24,zsignc)
resid6A100v$c65 <- sapply(X = resid6A100v$V30,zsignc)
resid6A100v$c66 <- sapply(X = resid6A100v$V36,zsignc)



zsig6A100<- resid6A100v[,38:73]

datosLong6A100<-pivot_longer(zsig6A100, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6A100$N <- gl(1,180000, labels = "100")

datosLong6A100$CATEGORY<- rep(nombres6,5000)

tabla6A100<-table(datosLong6A100$CATEGORY,datosLong6A100$MEASURE)
tabla6A100p<- as.data.frame(prop.table(tabla6A100,1))
colnames(tabla6A100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6A100p$CATEGORY<- rep(6,72)
tabla6A100p$MATRIX<- rep("A",72)
tabla6A100p$N<- rep(100,72)




## N = 100, matrix B

transiz6B100<-alply(seq6B100,1,transi)
dim6<- ldply(transiz6B100,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

#transiz6B100 <- transiz6B100[-falsos]
chisq6B100<- llply(transiz6B100,function(x) {loglm(~1+2, data=x)$pearson})
chi6B100<- llply(chisq6B100, as.vector)
chi26B100<- unlist(chi6B100)
buenos<- which(chi6B100 > qchisq(.95,25))
length(buenos)

transiz6B100<- transiz6B100[buenos]


resid6B100<- llply(transiz6B100, function(x) resid_ln(x))

#resid6B100 <- resid6B100[-falsos]
resid6B100v<- ldply(resid6B100,as.vector)
buenos<- complete.cases(resid6B100v)
resid6B100v<- resid6B100v[buenos,]
resid6B100v<- resid6B100v[1:5000,]


resid6B100v$c11 <- sapply(X = resid6B100v$V1,zsignc)
resid6B100v$c12 <- sapply(X = resid6B100v$V7,zsignc)
resid6B100v$c13 <- sapply(X = resid6B100v$V13,zsignc)
resid6B100v$c14 <- sapply(X = resid6B100v$V19,zsignc)
resid6B100v$c15 <- sapply(X = resid6B100v$V25,zsignc)
resid6B100v$c16 <- sapply(X = resid6B100v$V31,zsignc)


resid6B100v$c21 <- sapply(X = resid6B100v$V2,zsignc)
resid6B100v$c22 <- sapply(X = resid6B100v$V8,zsignc)
resid6B100v$c23 <- sapply(X = resid6B100v$V14,zsignc)
resid6B100v$c24 <- sapply(X = resid6B100v$V20,zsignc)
resid6B100v$c25 <- sapply(X = resid6B100v$V26,zsignc)
resid6B100v$c26 <- sapply(X = resid6B100v$V28,zsignc)


resid6B100v$c31 <- sapply(X = resid6B100v$V3,zsignc)
resid6B100v$c32 <- sapply(X = resid6B100v$V9,zsignc)
resid6B100v$c33 <- sapply(X = resid6B100v$V15,zsignc)
resid6B100v$c34 <- sapply(X = resid6B100v$V21,zsignc)
resid6B100v$c35 <- sapply(X = resid6B100v$V27,zsignc)
resid6B100v$c36 <- sapply(X = resid6B100v$V33,zsignc)


resid6B100v$c41 <- sapply(X = resid6B100v$V4,zsignc)
resid6B100v$c42 <- sapply(X = resid6B100v$V10,zsignc)
resid6B100v$c43 <- sapply(X = resid6B100v$V16,zsignc)
resid6B100v$c44 <- sapply(X = resid6B100v$V22,zsignc)
resid6B100v$c45 <- sapply(X = resid6B100v$V28,zsignc)
resid6B100v$c46 <- sapply(X = resid6B100v$V34,zsignc)


resid6B100v$c51 <- sapply(X = resid6B100v$V5,zsignc)
resid6B100v$c52 <- sapply(X = resid6B100v$V11,zsignc)
resid6B100v$c53 <- sapply(X = resid6B100v$V17,zsignc)
resid6B100v$c54 <- sapply(X = resid6B100v$V23,zsignc)
resid6B100v$c55 <- sapply(X = resid6B100v$V29,zsignc)
resid6B100v$c56 <- sapply(X = resid6B100v$V32,zsignc)



resid6B100v$c61 <- sapply(X = resid6B100v$V6,zsignc)
resid6B100v$c62 <- sapply(X = resid6B100v$V12,zsignc)
resid6B100v$c63 <- sapply(X = resid6B100v$V18,zsignc)
resid6B100v$c64 <- sapply(X = resid6B100v$V24,zsignc)
resid6B100v$c65 <- sapply(X = resid6B100v$V30,zsignc)
resid6B100v$c66 <- sapply(X = resid6B100v$V36,zsignc)

zsig6B100<- resid6B100v[,38:73]

datosLong6B100<-pivot_longer(zsig6B100, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6B100$N <- gl(1,180000, labels = "100")
datosLong6B100$CATEGORY<- rep(nombres6,5000)

tabla6B100<-table(datosLong6B100$CATEGORY,datosLong6B100$MEASURE)
tabla6B100p<- as.data.frame(prop.table(tabla6B100,1))
colnames(tabla6B100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6B100p$CATEGORY<- rep(6,108)
tabla6B100p$MATRIX<- rep("B",108)
tabla6B100p$N<- rep(100,108)


## N = 100, matrix C

transiz6C100<-alply(seq6C100,1,transi)

dim6<- ldply(transiz6C100,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

#transiz6C100 <- transiz6C100[-falsos]
chisq6C100<- llply(transiz6C100,function(x) {loglm(~1+2, data=x)$pearson})
chi6C100<- llply(chisq6C100, as.vector)
chi26C100<- unlist(chi6C100)
buenos<- which(chi6C100 > qchisq(.95,25))
length(buenos)

transiz6C100<- transiz6C100[buenos]
transiz6C100<- transiz6C100[1:5000]
resid6C100<- llply(transiz6C100, function(x) resid_ln(x))


resid6C100v<- ldply(resid6C100,as.vector)
buenos<- complete.cases(resid6C100v)
resid6C100v<- resid6C100v[buenos,]
resid6C100v<- resid6C100v[1:5000,]


resid6C100v$c11 <- sapply(X = resid6C100v$V1,zsignc)
resid6C100v$c12 <- sapply(X = resid6C100v$V7,zsignc)
resid6C100v$c13 <- sapply(X = resid6C100v$V13,zsignc)
resid6C100v$c14 <- sapply(X = resid6C100v$V19,zsignc)
resid6C100v$c15 <- sapply(X = resid6C100v$V25,zsignc)
resid6C100v$c16 <- sapply(X = resid6C100v$V31,zsignc)


resid6C100v$c21 <- sapply(X = resid6C100v$V2,zsignc)
resid6C100v$c22 <- sapply(X = resid6C100v$V8,zsignc)
resid6C100v$c23 <- sapply(X = resid6C100v$V14,zsignc)
resid6C100v$c24 <- sapply(X = resid6C100v$V20,zsignc)
resid6C100v$c25 <- sapply(X = resid6C100v$V26,zsignc)
resid6C100v$c26 <- sapply(X = resid6C100v$V28,zsignc)


resid6C100v$c31 <- sapply(X = resid6C100v$V3,zsignc)
resid6C100v$c32 <- sapply(X = resid6C100v$V9,zsignc)
resid6C100v$c33 <- sapply(X = resid6C100v$V15,zsignc)
resid6C100v$c34 <- sapply(X = resid6C100v$V21,zsignc)
resid6C100v$c35 <- sapply(X = resid6C100v$V27,zsignc)
resid6C100v$c36 <- sapply(X = resid6C100v$V33,zsignc)


resid6C100v$c41 <- sapply(X = resid6C100v$V4,zsignc)
resid6C100v$c42 <- sapply(X = resid6C100v$V10,zsignc)
resid6C100v$c43 <- sapply(X = resid6C100v$V16,zsignc)
resid6C100v$c44 <- sapply(X = resid6C100v$V22,zsignc)
resid6C100v$c45 <- sapply(X = resid6C100v$V28,zsignc)
resid6C100v$c46 <- sapply(X = resid6C100v$V34,zsignc)


resid6C100v$c51 <- sapply(X = resid6C100v$V5,zsignc)
resid6C100v$c52 <- sapply(X = resid6C100v$V11,zsignc)
resid6C100v$c53 <- sapply(X = resid6C100v$V17,zsignc)
resid6C100v$c54 <- sapply(X = resid6C100v$V23,zsignc)
resid6C100v$c55 <- sapply(X = resid6C100v$V29,zsignc)
resid6C100v$c56 <- sapply(X = resid6C100v$V32,zsignc)



resid6C100v$c61 <- sapply(X = resid6C100v$V6,zsignc)
resid6C100v$c62 <- sapply(X = resid6C100v$V12,zsignc)
resid6C100v$c63 <- sapply(X = resid6C100v$V18,zsignc)
resid6C100v$c64 <- sapply(X = resid6C100v$V24,zsignc)
resid6C100v$c65 <- sapply(X = resid6C100v$V30,zsignc)
resid6C100v$c66 <- sapply(X = resid6C100v$V36,zsignc)



zsig6C100<- resid6C100v[,38:73]

datosLong6C100<-pivot_longer(zsig6C100, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6C100$N <- gl(1,180000, labels = "100")

datosLong6C100$CATEGORY<- rep(nombres6,5000)

tabla6C100<-table(datosLong6C100$CATEGORY,datosLong6C100$MEASURE)
tabla6C100p<- as.data.frame(prop.table(tabla6C100,1))
colnames(tabla6C100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6C100p$CATEGORY<- rep(6,108)
tabla6C100p$MATRIX<- rep("C",108)
tabla6C100p$N<- rep(100,108)



## N = 150, matrix A

transiz6A150<-alply(seq6A150,1,transi)
dim6<- ldply(transiz6A150,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

#transiz6A150 <- transiz6A150[-falsos]
chisq6A150<- llply(transiz6A150,function(x) {loglm(~1+2, data=x)$pearson})
chi6A150<- llply(chisq6A150, as.vector)
chi26A150<- unlist(chi6A150)
buenos<- which(chi6A150 <= qchisq(.95,25))
length(buenos)

transiz6A150<- transiz6A150[buenos]

resid6A150<- llply(transiz6A150, function(x) resid_ln(x))

resid6A150v<- ldply(resid6A150,as.vector)
buenos<- complete.cases(resid6A150v)
resid6A150v<- resid6A150v[buenos,]
resid6A150v<- resid6A150v[1:5000,]


resid6A150v$c11 <- sapply(X = resid6A150v$V1,zsignc)
resid6A150v$c12 <- sapply(X = resid6A150v$V7,zsignc)
resid6A150v$c13 <- sapply(X = resid6A150v$V13,zsignc)
resid6A150v$c14 <- sapply(X = resid6A150v$V19,zsignc)
resid6A150v$c15 <- sapply(X = resid6A150v$V25,zsignc)
resid6A150v$c16 <- sapply(X = resid6A150v$V31,zsignc)


resid6A150v$c21 <- sapply(X = resid6A150v$V2,zsignc)
resid6A150v$c22 <- sapply(X = resid6A150v$V8,zsignc)
resid6A150v$c23 <- sapply(X = resid6A150v$V14,zsignc)
resid6A150v$c24 <- sapply(X = resid6A150v$V20,zsignc)
resid6A150v$c25 <- sapply(X = resid6A150v$V26,zsignc)
resid6A150v$c26 <- sapply(X = resid6A150v$V28,zsignc)


resid6A150v$c31 <- sapply(X = resid6A150v$V3,zsignc)
resid6A150v$c32 <- sapply(X = resid6A150v$V9,zsignc)
resid6A150v$c33 <- sapply(X = resid6A150v$V15,zsignc)
resid6A150v$c34 <- sapply(X = resid6A150v$V21,zsignc)
resid6A150v$c35 <- sapply(X = resid6A150v$V27,zsignc)
resid6A150v$c36 <- sapply(X = resid6A150v$V33,zsignc)


resid6A150v$c41 <- sapply(X = resid6A150v$V4,zsignc)
resid6A150v$c42 <- sapply(X = resid6A150v$V10,zsignc)
resid6A150v$c43 <- sapply(X = resid6A150v$V16,zsignc)
resid6A150v$c44 <- sapply(X = resid6A150v$V22,zsignc)
resid6A150v$c45 <- sapply(X = resid6A150v$V28,zsignc)
resid6A150v$c46 <- sapply(X = resid6A150v$V34,zsignc)


resid6A150v$c51 <- sapply(X = resid6A150v$V5,zsignc)
resid6A150v$c52 <- sapply(X = resid6A150v$V11,zsignc)
resid6A150v$c53 <- sapply(X = resid6A150v$V17,zsignc)
resid6A150v$c54 <- sapply(X = resid6A150v$V23,zsignc)
resid6A150v$c55 <- sapply(X = resid6A150v$V29,zsignc)
resid6A150v$c56 <- sapply(X = resid6A150v$V32,zsignc)



resid6A150v$c61 <- sapply(X = resid6A150v$V6,zsignc)
resid6A150v$c62 <- sapply(X = resid6A150v$V12,zsignc)
resid6A150v$c63 <- sapply(X = resid6A150v$V18,zsignc)
resid6A150v$c64 <- sapply(X = resid6A150v$V24,zsignc)
resid6A150v$c65 <- sapply(X = resid6A150v$V30,zsignc)
resid6A150v$c66 <- sapply(X = resid6A150v$V36,zsignc)



zsig6A150<- resid6A150v[,38:73]

datosLong6A150<-pivot_longer(zsig6A150, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong6A150$N <- gl(1,180000, labels = "150")

datosLong6A150$CATEGORY<- rep(nombres6,5000)

tabla6A150<-table(datosLong6A150$CATEGORY,datosLong6A150$MEASURE)
tabla6A150p<- as.data.frame(prop.table(tabla6A150,1))
colnames(tabla6A150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6A150p$CATEGORY<- rep(6,108)
tabla6A150p$MATRIX<- rep("A",108)
tabla6A150p$N<- rep(150,108)


## N = 150, matrix B

transiz6B150<-alply(seq6B150,1,transi)
dim6<- ldply(transiz6B150,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

#transiz6B150 <- transiz6B150[-falsos]
chisq6B150<- llply(transiz6B150,function(x) {loglm(~1+2, data=x)$pearson})
chi6B150<- llply(chisq6B150, as.vector)
chi26B150<- unlist(chi6B150)
buenos<- which(chi6B150 > qchisq(.95,25))
length(buenos)

transiz6B150<- transiz6B150[buenos]

resid6B150<- llply(transiz6B150, function(x) resid_ln(x))

resid6B150v<- ldply(resid6B150,as.vector)
buenos<- complete.cases(resid6B150v)
resid6B150v<- resid6B150v[buenos,]
resid6B150v<- resid6B150v[1:5000,]

resid6B150v$c11 <- sapply(X = resid6B150v$V1,zsignc)
resid6B150v$c12 <- sapply(X = resid6B150v$V7,zsignc)
resid6B150v$c13 <- sapply(X = resid6B150v$V13,zsignc)
resid6B150v$c14 <- sapply(X = resid6B150v$V19,zsignc)
resid6B150v$c15 <- sapply(X = resid6B150v$V25,zsignc)
resid6B150v$c16 <- sapply(X = resid6B150v$V31,zsignc)


resid6B150v$c21 <- sapply(X = resid6B150v$V2,zsignc)
resid6B150v$c22 <- sapply(X = resid6B150v$V8,zsignc)
resid6B150v$c23 <- sapply(X = resid6B150v$V14,zsignc)
resid6B150v$c24 <- sapply(X = resid6B150v$V20,zsignc)
resid6B150v$c25 <- sapply(X = resid6B150v$V26,zsignc)
resid6B150v$c26 <- sapply(X = resid6B150v$V28,zsignc)


resid6B150v$c31 <- sapply(X = resid6B150v$V3,zsignc)
resid6B150v$c32 <- sapply(X = resid6B150v$V9,zsignc)
resid6B150v$c33 <- sapply(X = resid6B150v$V15,zsignc)
resid6B150v$c34 <- sapply(X = resid6B150v$V21,zsignc)
resid6B150v$c35 <- sapply(X = resid6B150v$V27,zsignc)
resid6B150v$c36 <- sapply(X = resid6B150v$V33,zsignc)


resid6B150v$c41 <- sapply(X = resid6B150v$V4,zsignc)
resid6B150v$c42 <- sapply(X = resid6B150v$V10,zsignc)
resid6B150v$c43 <- sapply(X = resid6B150v$V16,zsignc)
resid6B150v$c44 <- sapply(X = resid6B150v$V22,zsignc)
resid6B150v$c45 <- sapply(X = resid6B150v$V28,zsignc)
resid6B150v$c46 <- sapply(X = resid6B150v$V34,zsignc)


resid6B150v$c51 <- sapply(X = resid6B150v$V5,zsignc)
resid6B150v$c52 <- sapply(X = resid6B150v$V11,zsignc)
resid6B150v$c53 <- sapply(X = resid6B150v$V17,zsignc)
resid6B150v$c54 <- sapply(X = resid6B150v$V23,zsignc)
resid6B150v$c55 <- sapply(X = resid6B150v$V29,zsignc)
resid6B150v$c56 <- sapply(X = resid6B150v$V32,zsignc)



resid6B150v$c61 <- sapply(X = resid6B150v$V6,zsignc)
resid6B150v$c62 <- sapply(X = resid6B150v$V12,zsignc)
resid6B150v$c63 <- sapply(X = resid6B150v$V18,zsignc)
resid6B150v$c64 <- sapply(X = resid6B150v$V24,zsignc)
resid6B150v$c65 <- sapply(X = resid6B150v$V30,zsignc)
resid6B150v$c66 <- sapply(X = resid6B150v$V36,zsignc)



zsig6B150<- resid6B150v[,38:73]

datosLong6B150<-pivot_longer(zsig6B150, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong6B150$N <- gl(1,180000, labels = "150")

datosLong6B150$CATEGORY<- rep(nombres6,5000)

tabla6B150<-table(datosLong6B150$CATEGORY,datosLong6B150$MEASURE)
tabla6B150p<- as.data.frame(prop.table(tabla6B150,1))
colnames(tabla6B150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6B150p$CATEGORY<- rep(6,108)
tabla6B150p$MATRIX<- rep("B",108)
tabla6B150p$N<- rep(150,108)


# N = 150 matrix C

transiz6C150 <-alply(seq6C150,1,transi)

dim6<- ldply(transiz6C150,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 != 6)
length(falsos)

#transiz6C150 <- transiz6C150[-falsos]
chisq6C150<- llply(transiz6C150,function(x) {loglm(~1+2, data=x)$pearson})
chi6C150<- llply(chisq6C150, as.vector)
chi26C150<- unlist(chi6C150)
buenos<- which(chi6C150 > qchisq(.95,25))
length(buenos)

transiz6C150<- transiz6C150[buenos]

resid6C150<- llply(transiz6C150, function(x) resid_ln(x))

resid6C150v<- ldply(resid6C150,as.vector)
buenos<- complete.cases(resid6C150v)
length(buenos)
resid6C150v<- resid6C150v[buenos,]
resid6C150v<- resid6C150v[1:5000,]


resid6C150v$c11 <- sapply(X = resid6C150v$V1,zsignc)
resid6C150v$c12 <- sapply(X = resid6C150v$V7,zsignc)
resid6C150v$c13 <- sapply(X = resid6C150v$V13,zsignc)
resid6C150v$c14 <- sapply(X = resid6C150v$V19,zsignc)
resid6C150v$c15 <- sapply(X = resid6C150v$V25,zsignc)
resid6C150v$c16 <- sapply(X = resid6C150v$V31,zsignc)


resid6C150v$c21 <- sapply(X = resid6C150v$V2,zsignc)
resid6C150v$c22 <- sapply(X = resid6C150v$V8,zsignc)
resid6C150v$c23 <- sapply(X = resid6C150v$V14,zsignc)
resid6C150v$c24 <- sapply(X = resid6C150v$V20,zsignc)
resid6C150v$c25 <- sapply(X = resid6C150v$V26,zsignc)
resid6C150v$c26 <- sapply(X = resid6C150v$V28,zsignc)


resid6C150v$c31 <- sapply(X = resid6C150v$V3,zsignc)
resid6C150v$c32 <- sapply(X = resid6C150v$V9,zsignc)
resid6C150v$c33 <- sapply(X = resid6C150v$V15,zsignc)
resid6C150v$c34 <- sapply(X = resid6C150v$V21,zsignc)
resid6C150v$c35 <- sapply(X = resid6C150v$V27,zsignc)
resid6C150v$c36 <- sapply(X = resid6C150v$V33,zsignc)


resid6C150v$c41 <- sapply(X = resid6C150v$V4,zsignc)
resid6C150v$c42 <- sapply(X = resid6C150v$V10,zsignc)
resid6C150v$c43 <- sapply(X = resid6C150v$V16,zsignc)
resid6C150v$c44 <- sapply(X = resid6C150v$V22,zsignc)
resid6C150v$c45 <- sapply(X = resid6C150v$V28,zsignc)
resid6C150v$c46 <- sapply(X = resid6C150v$V34,zsignc)


resid6C150v$c51 <- sapply(X = resid6C150v$V5,zsignc)
resid6C150v$c52 <- sapply(X = resid6C150v$V11,zsignc)
resid6C150v$c53 <- sapply(X = resid6C150v$V17,zsignc)
resid6C150v$c54 <- sapply(X = resid6C150v$V23,zsignc)
resid6C150v$c55 <- sapply(X = resid6C150v$V29,zsignc)
resid6C150v$c56 <- sapply(X = resid6C150v$V32,zsignc)

resid6C150v$c61 <- sapply(X = resid6C150v$V6,zsignc)
resid6C150v$c62 <- sapply(X = resid6C150v$V12,zsignc)
resid6C150v$c63 <- sapply(X = resid6C150v$V18,zsignc)
resid6C150v$c64 <- sapply(X = resid6C150v$V24,zsignc)
resid6C150v$c65 <- sapply(X = resid6C150v$V30,zsignc)
resid6C150v$c66 <- sapply(X = resid6C150v$V36,zsignc)



zsig6C150<- resid6C150v[,38:73]

datosLong6C150<-pivot_longer(zsig6C150, cols = 1:36,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong6C150$N <- gl(1,180000, labels = "150")

datosLong6C150$CATEGORY<- rep(nombres6,5000)

tabla6C150<-table(datosLong6C150$CATEGORY,datosLong6C150$MEASURE)
tabla6C150p<- as.data.frame(prop.table(tabla6C150,1))
colnames(tabla6C150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla6C150p$CATEGORY<- rep(6,108)
tabla6C150p$MATRIX<- rep("C",108)
tabla6C150p$N<- rep(150,108)



# Data summary table

DEPENDENCY<- gl(1,36,labels= "Inh")
PATTERN <- nombres6
P <- rep(0,36)
CATEGORY <- rep(6,36)
MATRIX<- gl(1,36,labels= "A")
N<- gl(1,36, labels = "20")

d20<- data.frame(PATTERN,DEPENDENCY,P,CATEGORY,MATRIX,N)

N<- gl(1,36, labels = "50")
d50<- data.frame(PATTERN,DEPENDENCY,P,CATEGORY,MATRIX,N)

d100<-data.frame(PATTERN,DEPENDENCY,P,CATEGORY,MATRIX,N)
d100$N<- gl(1,36, labels = "100")

d150<-data.frame(PATTERN,DEPENDENCY,P,CATEGORY,MATRIX,N)
d150$N<- gl(1,36, labels = "150")

# tabla6A20p<- rbind(tabla6A20p,d20)
 tabla6A50p<- rbind(tabla6A50p,d50)
 tabla6A100p<- rbind(tabla6A100p,d100)
# tabla6A150p<- rbind(tabla6A150p,d150)
# tabla6B20p<- rbind(tabla6B20p,d20)
# tabla6C20p<- rbind(tabla6C20p,d20)




tabla6_lm_bon<- rbind(tabla6A50p,tabla6A100p,tabla6A150p,
                      tabla6B50p,tabla6B100p,tabla6B150p, 
                      tabla6C50p,tabla6C100p,tabla6C150p)



# Activation patterns

tab6bc5<- subset(tabla6_lm_bon,PATTERN == "BC" & DEPENDENCY == "Act" )
tab6dc5<- subset(tabla6_lm_bon,PATTERN == "DC" & DEPENDENCY == "Act" )
tab6_lm_bon<- rbind(tab6bc5,tab6dc5)
tab6_lm_bon$MODEL<- gl(1,18, labels = "BON-LOGL")
tab6_lm_bon$N<- rep(gl(3,1, labels = c("50","100","150")),6)

# Inhibition patterns


tab6ba5<- subset(tabla6_lm_bon,PATTERN == "BA" & DEPENDENCY == "Inh" )
tab6bb5<- subset(tabla6_lm_bon,PATTERN == "BB" & DEPENDENCY == "Inh" )
tab6bd5<- subset(tabla6_lm_bon,PATTERN == "BD" & DEPENDENCY == "Inh" )
tab6be5<- subset(tabla6_lm_bon,PATTERN == "BE" & DEPENDENCY == "Inh" )
tab6bf5<- subset(tabla6_lm_bon,PATTERN == "BF" & DEPENDENCY == "Inh" )


tab6da5<- subset(tabla6_lm_bon,PATTERN == "DA" & DEPENDENCY == "Inh" )
tab6db5<- subset(tabla6_lm_bon,PATTERN == "DB" & DEPENDENCY == "Inh" )
tab6dd5<- subset(tabla6_lm_bon,PATTERN == "DD" & DEPENDENCY == "Inh" )
tab6de5<- subset(tabla6_lm_bon,PATTERN == "DE" & DEPENDENCY == "Inh" )
tab6df5<- subset(tabla6_lm_bon,PATTERN == "DF" & DEPENDENCY == "Inh" )

tab6_inh_lm_bon<- rbind(tab6ba5, tab6bb5,tab6bd5,tab6be5,tab6bf5,
                        tab6da5,tab6db5,tab6dd5,tab6de5,tab6df5)
tab6_inh_lm_bon$MODEL<- gl(1,90, labels = "BON-LOGL")
tab6_inh_lm_bon$N<- rep(gl(3,1, labels = c("50","100","150")),30)

# Spurious pattern

tab6ac1<- subset(tabla6_lm_bon,PATTERN == "AC" & DEPENDENCY == "Inh" )
tab6cc1<- subset(tabla6_lm_bon,PATTERN == "CC" & DEPENDENCY == "Inh" )

tab6_sp_lm_bon<- rbind(tab6ac1,tab6cc1)
tab6_sp_lm_bon$MODEL<- gl(1,18, labels = "BON-LOGL")
tab6_sp_lm_bon$N<- rep(gl(3,1, labels = c("50","100","150")),6)


end <- Sys.time()
print(end-start)
