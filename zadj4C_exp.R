start <- Sys.time()

library(simstudy)
library(MCTM)
library(plyr)
library(ggplot2)
#library(dplyr)
library(tidyr)
library(misty)

set.seed(123456)

# Functions 
 
## Calculation of transition frequencies

transi<- function(x) {
  TransMatrix(as.numeric(x),order= 1,probs= F)
}



# Significance of the residuals

zsign<-function(x){
  ifelse(abs(x) > 1.96,1,0) 
}


zsignc<-function(x){
  ifelse(x > 1.96, "Act",ifelse(x < -1.96, "Inh", "Ind"))
}



# Matrices with 4 states


## Matrix   A


m4A <- matrix(c(1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4), 4, 4, 
              byrow = T)

## Matrix B p =3/4 states = 4, 1 significant (BC)

m4B <- matrix(c(1/4, 1/4, 1/4,1/4,
                1/12, 1/12, 9/12,1/12,
                1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4), 4, 4, 
              byrow = T)


## Matrix C p =3/4 states = 4, 2 significant (BC y DC)

m4C <- matrix(c( 1/4, 1/4, 1/4,1/4,
                 1/12, 1/12, 9/12,1/12,
                 1/4, 1/4, 1/4,1/4,
                 1/12, 1/12, 9/12,1/12), 4, 4, 
              byrow = T)



## names4

nombres4<- c("AA","AB","AC","AD",
            "BA","BB","BC","BD",
            "CA","CB","CC","CD",
            "DA","DB","DC","DD")
            

# Generation of the sequences

## Sequences with 4 states

seq4A20 <- genMarkov(n = 10000, transMat = m4A, 
                     chainLen = 20, wide = TRUE)
seq4B20 <- genMarkov(n = 15000, transMat = m4B, 
                     chainLen = 20, wide = TRUE)
seq4C20 <- genMarkov(n = 15000, transMat = m4C, 
                     chainLen = 20, wide = TRUE)


seq4A50 <- genMarkov(n = 15000, transMat = m4A, 
                     chainLen = 50, wide = TRUE)
seq4B50 <- genMarkov(n = 15000, transMat = m4B, 
                     chainLen = 50, wide = TRUE)
seq4C50 <- genMarkov(n = 15000, transMat = m4C, 
                     chainLen = 50, wide = TRUE)


seq4A100 <- genMarkov(n = 10000, transMat = m4A, 
                      chainLen = 100, wide = TRUE)
seq4B100 <- genMarkov(n = 10000, transMat = m4B, 
                      chainLen = 100, wide = TRUE)
seq4C100 <- genMarkov(n = 10000, transMat = m4C, 
                      chainLen = 100, wide = TRUE)


seq4A150 <- genMarkov(n = 10000, transMat = m4A, 
                      chainLen = 150, wide = TRUE)
seq4B150 <- genMarkov(n = 10000, transMat = m4B, 
                      chainLen = 150, wide = TRUE)
seq4C150 <- genMarkov(n = 10000, transMat = m4C, 
                      chainLen = 150, wide = TRUE)


seq4A20<- seq4A20[,-1]
seq4B20<- seq4B20[,-1]
seq4C20<- seq4C20[,-1]

seq4A50<- seq4A50[,-1]
seq4B50<- seq4B50[,-1]
seq4C50<- seq4C50[,-1]

seq4A100<- seq4A100[,-1]
seq4B100<- seq4B100[,-1]
seq4C100<- seq4C100[,-1]


seq4A150<- seq4A150[,-1]
seq4B150<- seq4B150[,-1]
seq4C150<- seq4C150[,-1]


# Residuals of the 4-category condition

## N = 20, matrix A

transiz4A20<-alply(seq4A20,1,transi)
resid4A20<- llply(transiz4A20, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4A20,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)


resid4A20 <- resid4A20[-falsos]
resid4A20v<- ldply(resid4A20,as.vector)
buenos<- complete.cases(resid4A20v)
length(buenos)
resid4A20v<- resid4A20v[buenos,]
resid4A20v<- resid4A20v[1:5000,]


resid4A20v$c11 <- sapply(X = resid4A20v$V1,zsignc)
resid4A20v$c12 <- sapply(X = resid4A20v$V5,zsignc)
resid4A20v$c13 <- sapply(X = resid4A20v$V9,zsignc)
resid4A20v$c14 <- sapply(X = resid4A20v$V13,zsignc)

resid4A20v$c21 <- sapply(X = resid4A20v$V2,zsignc)
resid4A20v$c22 <- sapply(X = resid4A20v$V6,zsignc)
resid4A20v$c23 <- sapply(X = resid4A20v$V10,zsignc)
resid4A20v$c24 <- sapply(X = resid4A20v$V14,zsignc)

resid4A20v$c31 <- sapply(X = resid4A20v$V3,zsignc)
resid4A20v$c32 <- sapply(X = resid4A20v$V7,zsignc)
resid4A20v$c33 <- sapply(X = resid4A20v$V11,zsignc)
resid4A20v$c34 <- sapply(X = resid4A20v$V15,zsignc)

resid4A20v$c41 <- sapply(X = resid4A20v$V4,zsignc)
resid4A20v$c42 <- sapply(X = resid4A20v$V8,zsignc)
resid4A20v$c43 <- sapply(X = resid4A20v$V12,zsignc)
resid4A20v$c44 <- sapply(X = resid4A20v$V16,zsignc)




zsig4A20<- resid4A20v[,18:33]

datosLong4A20<-pivot_longer(zsig4A20, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong4A20$N <- gl(1,80000, labels = "20")
datosLong4A20$CATEGORY<- rep(nombres4,5000)

tabla4A20<-table(datosLong4A20$CATEGORY,datosLong4A20$MEASURE)
tabla4A20p<- as.data.frame(prop.table(tabla4A20,1))
colnames(tabla4A20p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4A20p$CATEGORY<- rep(4,48)
tabla4A20p$MATRIX<- rep("A",48)
tabla4A20p$N<- rep(20,48)


# N = 20 MATRIX B

transiz4B20<-alply(seq4B20,1,transi)
resid4B20<- llply(transiz4B20, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4B20,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4B20 <- resid4B20[-falsos]
resid4B20v<- ldply(resid4B20,as.vector)
buenos<- complete.cases(resid4B20v)
resid4B20v<- resid4B20v[buenos,]
resid4B20v<- resid4B20v[1:5000,]



resid4B20v$c11 <- sapply(X = resid4B20v$V1,zsignc)
resid4B20v$c12 <- sapply(X = resid4B20v$V5,zsignc)
resid4B20v$c13 <- sapply(X = resid4B20v$V9,zsignc)
resid4B20v$c14 <- sapply(X = resid4B20v$V13,zsignc)

resid4B20v$c21 <- sapply(X = resid4B20v$V2,zsignc)
resid4B20v$c22 <- sapply(X = resid4B20v$V6,zsignc)
resid4B20v$c23 <- sapply(X = resid4B20v$V10,zsignc)
resid4B20v$c24 <- sapply(X = resid4B20v$V14,zsignc)

resid4B20v$c31 <- sapply(X = resid4B20v$V3,zsignc)
resid4B20v$c32 <- sapply(X = resid4B20v$V7,zsignc)
resid4B20v$c33 <- sapply(X = resid4B20v$V11,zsignc)
resid4B20v$c34 <- sapply(X = resid4B20v$V15,zsignc)

resid4B20v$c41 <- sapply(X = resid4B20v$V4,zsignc)
resid4B20v$c42 <- sapply(X = resid4B20v$V8,zsignc)
resid4B20v$c43 <- sapply(X = resid4B20v$V12,zsignc)
resid4B20v$c44 <- sapply(X = resid4B20v$V16,zsignc)




zsig4B20<- resid4B20v[,18:33]

datosLong4B20<-pivot_longer(zsig4B20, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong4B20$N <- gl(1,80000, labels = "20")

datosLong4B20$CATEGORY<- rep(nombres4,5000)



tabla4B20<-table(datosLong4B20$CATEGORY,datosLong4B20$MEASURE)
tabla4B20p<- as.data.frame(prop.table(tabla4B20,1))
colnames(tabla4B20p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4B20p$CATEGORY<- rep(4,48)
tabla4B20p$MATRIX<- rep("B",48)
tabla4B20p$N<- rep(20,48)


tabla4_20p<- rbind(tabla4A20p,tabla4B20p,tabla4B20p)



# N = 20 Matrix C

transiz4C20<-alply(seq4C20,1,transi)
resid4C20<- llply(transiz4C20, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4C20,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4C20 <- resid4C20[-falsos]
resid4C20v<- ldply(resid4C20,as.vector)
buenos<- complete.cases(resid4C20v)
length(buenos)
resid4C20v<- resid4C20v[buenos,]
resid4C20v<- resid4C20v[1:5000,]


resid4C20v$c11 <- sapply(X = resid4C20v$V1,zsignc)
resid4C20v$c12 <- sapply(X = resid4C20v$V5,zsignc)
resid4C20v$c13 <- sapply(X = resid4C20v$V9,zsignc)
resid4C20v$c14 <- sapply(X = resid4C20v$V13,zsignc)

resid4C20v$c21 <- sapply(X = resid4C20v$V2,zsignc)
resid4C20v$c22 <- sapply(X = resid4C20v$V6,zsignc)
resid4C20v$c23 <- sapply(X = resid4C20v$V10,zsignc)
resid4C20v$c24 <- sapply(X = resid4C20v$V14,zsignc)

resid4C20v$c31 <- sapply(X = resid4C20v$V3,zsignc)
resid4C20v$c32 <- sapply(X = resid4C20v$V7,zsignc)
resid4C20v$c33 <- sapply(X = resid4C20v$V11,zsignc)
resid4C20v$c34 <- sapply(X = resid4C20v$V15,zsignc)

resid4C20v$c41 <- sapply(X = resid4C20v$V4,zsignc)
resid4C20v$c42 <- sapply(X = resid4C20v$V8,zsignc)
resid4C20v$c43 <- sapply(X = resid4C20v$V12,zsignc)
resid4C20v$c44 <- sapply(X = resid4C20v$V16,zsignc)



zsig4C20<- resid4C20v[,18:33]


datosLong4C20<-pivot_longer(zsig4C20, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong4C20$N <- gl(1,80000, labels = "20")

datosLong4C20$CATEGORY<- rep(nombres4,5000)



tabla4C20<-table(datosLong4C20$CATEGORY,datosLong4C20$MEASURE)
tabla4C20p<- as.data.frame(prop.table(tabla4C20,1))
colnames(tabla4C20p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4C20p$CATEGORY<- rep(4,48)
tabla4C20p$MATRIX<- rep("C",48)
tabla4C20p$N<- rep(20,48)


tabla4_20p<- rbind(tabla4A20p,tabla4B20p,tabla4C20p)


# N = 50 Matrix A

transiz4A50<-alply(seq4A50,1,transi)
resid4A50<- llply(transiz4A50, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4A50,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)


resid4A50v<- ldply(resid4A50,as.vector)
buenos<- complete.cases(resid4A50v)
resid4A50v<- resid4A50v[buenos,]
resid4A50v<- resid4A50v[1:5000,]


resid4A50v$c11 <- sapply(X = resid4A50v$V1,zsignc)
resid4A50v$c12 <- sapply(X = resid4A50v$V5,zsignc)
resid4A50v$c13 <- sapply(X = resid4A50v$V9,zsignc)
resid4A50v$c14 <- sapply(X = resid4A50v$V13,zsignc)

resid4A50v$c21 <- sapply(X = resid4A50v$V2,zsignc)
resid4A50v$c22 <- sapply(X = resid4A50v$V6,zsignc)
resid4A50v$c23 <- sapply(X = resid4A50v$V10,zsignc)
resid4A50v$c24 <- sapply(X = resid4A50v$V14,zsignc)

resid4A50v$c31 <- sapply(X = resid4A50v$V3,zsignc)
resid4A50v$c32 <- sapply(X = resid4A50v$V7,zsignc)
resid4A50v$c33 <- sapply(X = resid4A50v$V11,zsignc)
resid4A50v$c34 <- sapply(X = resid4A50v$V15,zsignc)

resid4A50v$c41 <- sapply(X = resid4A50v$V4,zsignc)
resid4A50v$c42 <- sapply(X = resid4A50v$V8,zsignc)
resid4A50v$c43 <- sapply(X = resid4A50v$V12,zsignc)
resid4A50v$c44 <- sapply(X = resid4A50v$V16,zsignc)


zsig4A50<- resid4A50v[,18:33]

datosLong4A50<-pivot_longer(zsig4A50, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong4A50$N <- gl(1,80000, labels = "20")

datosLong4A50$CATEGORY<- rep(nombres4,5000)

tabla4A50<-table(datosLong4A50$CATEGORY,datosLong4A50$MEASURE)
tabla4A50p<- as.data.frame(prop.table(tabla4A50,1))
colnames(tabla4A50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4A50p$CATEGORY<- rep(4,48)
tabla4A50p$MATRIX<- rep("A",48)
tabla4A50p$N<- rep(50,48)



# N = 50 Matrix B

transiz4B50<-alply(seq4B50,1,transi)
resid4B50<- llply(transiz4B50, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4B50,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4B50v<- ldply(resid4B50,as.vector)
buenos<- complete.cases(resid4B50v)
resid4B50v<- resid4B50v[buenos,]
resid4B50v<- resid4B50v[1:5000,]



resid4B50v$c11 <- sapply(X = resid4B50v$V1,zsignc)
resid4B50v$c12 <- sapply(X = resid4B50v$V5,zsignc)
resid4B50v$c13 <- sapply(X = resid4B50v$V9,zsignc)
resid4B50v$c14 <- sapply(X = resid4B50v$V13,zsignc)

resid4B50v$c21 <- sapply(X = resid4B50v$V2,zsignc)
resid4B50v$c22 <- sapply(X = resid4B50v$V6,zsignc)
resid4B50v$c23 <- sapply(X = resid4B50v$V10,zsignc)
resid4B50v$c24 <- sapply(X = resid4B50v$V14,zsignc)

resid4B50v$c31 <- sapply(X = resid4B50v$V3,zsignc)
resid4B50v$c32 <- sapply(X = resid4B50v$V7,zsignc)
resid4B50v$c33 <- sapply(X = resid4B50v$V11,zsignc)
resid4B50v$c34 <- sapply(X = resid4B50v$V15,zsignc)

resid4B50v$c41 <- sapply(X = resid4B50v$V4,zsignc)
resid4B50v$c42 <- sapply(X = resid4B50v$V8,zsignc)
resid4B50v$c43 <- sapply(X = resid4B50v$V12,zsignc)
resid4B50v$c44 <- sapply(X = resid4B50v$V16,zsignc)




zsig4B50<- resid4B50v[,18:33]

datosLong4B50<-pivot_longer(zsig4B50, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong4B50$N <- gl(1,80000, labels = "20")
datosLong4B50$CATEGORY<- rep(nombres4,5000)

tabla4B50<-table(datosLong4B50$CATEGORY,datosLong4B50$MEASURE)
tabla4B50p<- as.data.frame(prop.table(tabla4B50,1))
colnames(tabla4B50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4B50p$CATEGORY<- rep(4,48)
tabla4B50p$MATRIX<- rep("B",48)
tabla4B50p$N<- rep(50,48)




# N = 50 Matrix C

transiz4C50<-alply(seq4C50,1,transi)
resid4C50<- llply(transiz4C50, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4C50,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4C50v<- ldply(resid4C50,as.vector)
buenos<- complete.cases(resid4C50v)
resid4C50v<- resid4C50v[buenos,]
resid4C50v<- resid4C50v[1:5000,]


resid4C50v$c11 <- sapply(X = resid4C50v$V1,zsignc)
resid4C50v$c12 <- sapply(X = resid4C50v$V5,zsignc)
resid4C50v$c13 <- sapply(X = resid4C50v$V9,zsignc)
resid4C50v$c14 <- sapply(X = resid4C50v$V13,zsignc)

resid4C50v$c21 <- sapply(X = resid4C50v$V2,zsignc)
resid4C50v$c22 <- sapply(X = resid4C50v$V6,zsignc)
resid4C50v$c23 <- sapply(X = resid4C50v$V10,zsignc)
resid4C50v$c24 <- sapply(X = resid4C50v$V14,zsignc)

resid4C50v$c31 <- sapply(X = resid4C50v$V3,zsignc)
resid4C50v$c32 <- sapply(X = resid4C50v$V7,zsignc)
resid4C50v$c33 <- sapply(X = resid4C50v$V11,zsignc)
resid4C50v$c34 <- sapply(X = resid4C50v$V15,zsignc)

resid4C50v$c41 <- sapply(X = resid4C50v$V4,zsignc)
resid4C50v$c42 <- sapply(X = resid4C50v$V8,zsignc)
resid4C50v$c43 <- sapply(X = resid4C50v$V12,zsignc)
resid4C50v$c44 <- sapply(X = resid4C50v$V16,zsignc)



zsig4C50<- resid4C50v[,18:33]

datosLong4C50<-pivot_longer(zsig4C50, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "CATEGORY"))

datosLong4C50$N <- gl(1,80000, labels = "20")
datosLong4C50$CATEGORY<- rep(nombres4,5000)

tabla4C50<-table(datosLong4C50$CATEGORY,datosLong4C50$MEASURE)
tabla4C50p<- as.data.frame(prop.table(tabla4C50,1))
colnames(tabla4C50p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4C50p$CATEGORY<- rep(4,48)
tabla4C50p$MATRIX<- rep("C",48)
tabla4C50p$N<- rep(50,48)

tabla4_50p<- rbind(tabla4A50p,tabla4B50p,tabla4C50p)


# N = 100 Matrix A

transiz4A100<-alply(seq4A100,1,transi)
resid4A100<- llply(transiz4A100, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4A100,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)


resid4A100v<- ldply(resid4A100,as.vector)
buenos<- complete.cases(resid4A100v)
resid4A100v<- resid4A100v[buenos,]
resid4A100v<- resid4A100v[1:5000,]


resid4A100v$c11 <- sapply(X = resid4A100v$V1,zsignc)
resid4A100v$c12 <- sapply(X = resid4A100v$V5,zsignc)
resid4A100v$c13 <- sapply(X = resid4A100v$V9,zsignc)
resid4A100v$c14 <- sapply(X = resid4A100v$V13,zsignc)

resid4A100v$c21 <- sapply(X = resid4A100v$V2,zsignc)
resid4A100v$c22 <- sapply(X = resid4A100v$V6,zsignc)
resid4A100v$c23 <- sapply(X = resid4A100v$V10,zsignc)
resid4A100v$c24 <- sapply(X = resid4A100v$V14,zsignc)

resid4A100v$c31 <- sapply(X = resid4A100v$V3,zsignc)
resid4A100v$c32 <- sapply(X = resid4A100v$V7,zsignc)
resid4A100v$c33 <- sapply(X = resid4A100v$V11,zsignc)
resid4A100v$c34 <- sapply(X = resid4A100v$V15,zsignc)

resid4A100v$c41 <- sapply(X = resid4A100v$V4,zsignc)
resid4A100v$c42 <- sapply(X = resid4A100v$V8,zsignc)
resid4A100v$c43 <- sapply(X = resid4A100v$V12,zsignc)
resid4A100v$c44 <- sapply(X = resid4A100v$V16,zsignc)




zsig4A100<- resid4A100v[,18:33]

datosLong4A100<-pivot_longer(zsig4A100, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong4A100$N <- gl(1,80000, labels = "20")
datosLong4A100$CATEGORY<- rep(nombres4,5000)

tabla4A100<-table(datosLong4A100$CATEGORY,datosLong4A100$MEASURE)
tabla4A100p<- as.data.frame(prop.table(tabla4A100,1))
colnames(tabla4A100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4A100p$CATEGORY<- rep(4,48)
tabla4A100p$MATRIX<- rep("A",48)
tabla4A100p$N<- rep(100,48)



# N = 100 Matrix B

transiz4B100<-alply(seq4B100,1,transi)
resid4B100<- llply(transiz4B100, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4B100,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4B100v<- ldply(resid4B100,as.vector)
buenos<- complete.cases(resid4B100v)
resid4B100v<- resid4B100v[buenos,]
resid4B100v<- resid4B100v[1:5000,]



resid4B100v$c11 <- sapply(X = resid4B100v$V1,zsignc)
resid4B100v$c12 <- sapply(X = resid4B100v$V5,zsignc)
resid4B100v$c13 <- sapply(X = resid4B100v$V9,zsignc)
resid4B100v$c14 <- sapply(X = resid4B100v$V13,zsignc)

resid4B100v$c21 <- sapply(X = resid4B100v$V2,zsignc)
resid4B100v$c22 <- sapply(X = resid4B100v$V6,zsignc)
resid4B100v$c23 <- sapply(X = resid4B100v$V10,zsignc)
resid4B100v$c24 <- sapply(X = resid4B100v$V14,zsignc)

resid4B100v$c31 <- sapply(X = resid4B100v$V3,zsignc)
resid4B100v$c32 <- sapply(X = resid4B100v$V7,zsignc)
resid4B100v$c33 <- sapply(X = resid4B100v$V11,zsignc)
resid4B100v$c34 <- sapply(X = resid4B100v$V15,zsignc)

resid4B100v$c41 <- sapply(X = resid4B100v$V4,zsignc)
resid4B100v$c42 <- sapply(X = resid4B100v$V8,zsignc)
resid4B100v$c43 <- sapply(X = resid4B100v$V12,zsignc)
resid4B100v$c44 <- sapply(X = resid4B100v$V16,zsignc)




zsig4B100<- resid4B100v[,18:33]

datosLong4B100<-pivot_longer(zsig4B100, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))
datosLong4B100$N <- gl(1,80000, labels = "20")
datosLong4B100$CATEGORY<- rep(nombres4,5000)

tabla4B100<-table(datosLong4B100$CATEGORY,datosLong4B100$MEASURE)
tabla4B100p<- as.data.frame(prop.table(tabla4B100,1))
colnames(tabla4B100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4B100p$CATEGORY<- rep(4,48)
tabla4B100p$MATRIX<- rep("B",48)
tabla4B100p$N<- rep(100,48)


# N = 100 Matrix C

transiz4C100<-alply(seq4C100,1,transi)
resid4C100<- llply(transiz4C100, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4C100,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4C100v<- ldply(resid4C100,as.vector)
buenos<- complete.cases(resid4C100v)
resid4C100v<- resid4C100v[buenos,]
resid4C100v<- resid4C100v[1:5000,]


resid4C100v$c11 <- sapply(X = resid4C100v$V1,zsignc)
resid4C100v$c12 <- sapply(X = resid4C100v$V5,zsignc)
resid4C100v$c13 <- sapply(X = resid4C100v$V9,zsignc)
resid4C100v$c14 <- sapply(X = resid4C100v$V13,zsignc)

resid4C100v$c21 <- sapply(X = resid4C100v$V2,zsignc)
resid4C100v$c22 <- sapply(X = resid4C100v$V6,zsignc)
resid4C100v$c23 <- sapply(X = resid4C100v$V10,zsignc)
resid4C100v$c24 <- sapply(X = resid4C100v$V14,zsignc)

resid4C100v$c31 <- sapply(X = resid4C100v$V3,zsignc)
resid4C100v$c32 <- sapply(X = resid4C100v$V7,zsignc)
resid4C100v$c33 <- sapply(X = resid4C100v$V11,zsignc)
resid4C100v$c34 <- sapply(X = resid4C100v$V15,zsignc)

resid4C100v$c41 <- sapply(X = resid4C100v$V4,zsignc)
resid4C100v$c42 <- sapply(X = resid4C100v$V8,zsignc)
resid4C100v$c43 <- sapply(X = resid4C100v$V12,zsignc)
resid4C100v$c44 <- sapply(X = resid4C100v$V16,zsignc)



zsig4C100<- resid4C100v[,18:33]

datosLong4C100<-pivot_longer(zsig4C100, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong4C100$N <- gl(1,80000, labels = "20")
datosLong4C100$CATEGORY<- rep(nombres4,5000)

tabla4C100<-table(datosLong4C100$CATEGORY,datosLong4C100$MEASURE)
tabla4C100p<- as.data.frame(prop.table(tabla4C100,1))
colnames(tabla4C100p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4C100p$CATEGORY<- rep(4,48)
tabla4C100p$MATRIX<- rep("C",48)
tabla4C100p$N<- rep(100,48)


tabla4_100p<- rbind(tabla4A100p,tabla4B100p,tabla4C100p)



# N = 150 Matrix A

transiz4A150<-alply(seq4A150,1,transi)
resid4A150<- llply(transiz4A150, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4A150,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)


resid4A150v<- ldply(resid4A150,as.vector)
buenos<- complete.cases(resid4A150v)
resid4A150v<- resid4A150v[buenos,]
resid4A150v<- resid4A150v[1:5000,]


resid4A150v$c11 <- sapply(X = resid4A150v$V1,zsignc)
resid4A150v$c12 <- sapply(X = resid4A150v$V5,zsignc)
resid4A150v$c13 <- sapply(X = resid4A150v$V9,zsignc)
resid4A150v$c14 <- sapply(X = resid4A150v$V13,zsignc)

resid4A150v$c21 <- sapply(X = resid4A150v$V2,zsignc)
resid4A150v$c22 <- sapply(X = resid4A150v$V6,zsignc)
resid4A150v$c23 <- sapply(X = resid4A150v$V10,zsignc)
resid4A150v$c24 <- sapply(X = resid4A150v$V14,zsignc)

resid4A150v$c31 <- sapply(X = resid4A150v$V3,zsignc)
resid4A150v$c32 <- sapply(X = resid4A150v$V7,zsignc)
resid4A150v$c33 <- sapply(X = resid4A150v$V11,zsignc)
resid4A150v$c34 <- sapply(X = resid4A150v$V15,zsignc)

resid4A150v$c41 <- sapply(X = resid4A150v$V4,zsignc)
resid4A150v$c42 <- sapply(X = resid4A150v$V8,zsignc)
resid4A150v$c43 <- sapply(X = resid4A150v$V12,zsignc)
resid4A150v$c44 <- sapply(X = resid4A150v$V16,zsignc)




zsig4A150<- resid4A150v[,18:33]

datosLong4A150<-pivot_longer(zsig4A150, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong4A150$N <- gl(1,80000, labels = "20")
datosLong4A150$CATEGORY<- rep(nombres4,5000)

tabla4A150<-table(datosLong4A150$CATEGORY,datosLong4A150$MEASURE)
tabla4A150p<- as.data.frame(prop.table(tabla4A150,1))
colnames(tabla4A150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4A150p$CATEGORY<- rep(4,48)
tabla4A150p$MATRIX<- rep("A",48)
tabla4A150p$N<- rep(150,48)



# N = 150 Matrix B

transiz4B150<-alply(seq4B150,1,transi)
resid4B150<- llply(transiz4B150, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4B150,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4B150v<- ldply(resid4B150,as.vector)
buenos<- complete.cases(resid4B150v)
resid4B150v<- resid4B150v[buenos,]
resid4B150v<- resid4B150v[1:5000,]



resid4B150v$c11 <- sapply(X = resid4B150v$V1,zsignc)
resid4B150v$c12 <- sapply(X = resid4B150v$V5,zsignc)
resid4B150v$c13 <- sapply(X = resid4B150v$V9,zsignc)
resid4B150v$c14 <- sapply(X = resid4B150v$V13,zsignc)

resid4B150v$c21 <- sapply(X = resid4B150v$V2,zsignc)
resid4B150v$c22 <- sapply(X = resid4B150v$V6,zsignc)
resid4B150v$c23 <- sapply(X = resid4B150v$V10,zsignc)
resid4B150v$c24 <- sapply(X = resid4B150v$V14,zsignc)

resid4B150v$c31 <- sapply(X = resid4B150v$V3,zsignc)
resid4B150v$c32 <- sapply(X = resid4B150v$V7,zsignc)
resid4B150v$c33 <- sapply(X = resid4B150v$V11,zsignc)
resid4B150v$c34 <- sapply(X = resid4B150v$V15,zsignc)

resid4B150v$c41 <- sapply(X = resid4B150v$V4,zsignc)
resid4B150v$c42 <- sapply(X = resid4B150v$V8,zsignc)
resid4B150v$c43 <- sapply(X = resid4B150v$V12,zsignc)
resid4B150v$c44 <- sapply(X = resid4B150v$V16,zsignc)



zsig4B150<- resid4B150v[,18:33]

datosLong4B150<-pivot_longer(zsig4B150, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))
datosLong4B150$N <- gl(1,80000, labels = "20")
datosLong4B150$CATEGORY<- rep(nombres4,5000)

tabla4B150<-table(datosLong4B150$CATEGORY,datosLong4B150$MEASURE)
tabla4B150p<- as.data.frame(prop.table(tabla4B150,1))
colnames(tabla4B150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4B150p$CATEGORY<- rep(4,48)
tabla4B150p$MATRIX<- rep("B",48)
tabla4B150p$N<- rep(150,48)


# N = 150 Matrix C

transiz4C150<-alply(seq4C150,1,transi)
resid4C150<- llply(transiz4C150, function(x) chisq.test(x)$stdres)

dim4<- ldply(resid4C150,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 != 4)
length(falsos)

resid4C150v<- ldply(resid4C150,as.vector)
buenos<- complete.cases(resid4C150v)
resid4C150v<- resid4C150v[buenos,]
resid4C150v<- resid4C150v[1:5000,]


resid4C150v$c11 <- sapply(X = resid4C150v$V1,zsignc)
resid4C150v$c12 <- sapply(X = resid4C150v$V5,zsignc)
resid4C150v$c13 <- sapply(X = resid4C150v$V9,zsignc)
resid4C150v$c14 <- sapply(X = resid4C150v$V13,zsignc)

resid4C150v$c21 <- sapply(X = resid4C150v$V2,zsignc)
resid4C150v$c22 <- sapply(X = resid4C150v$V6,zsignc)
resid4C150v$c23 <- sapply(X = resid4C150v$V10,zsignc)
resid4C150v$c24 <- sapply(X = resid4C150v$V14,zsignc)

resid4C150v$c31 <- sapply(X = resid4C150v$V3,zsignc)
resid4C150v$c32 <- sapply(X = resid4C150v$V7,zsignc)
resid4C150v$c33 <- sapply(X = resid4C150v$V11,zsignc)
resid4C150v$c34 <- sapply(X = resid4C150v$V15,zsignc)

resid4C150v$c41 <- sapply(X = resid4C150v$V4,zsignc)
resid4C150v$c42 <- sapply(X = resid4C150v$V8,zsignc)
resid4C150v$c43 <- sapply(X = resid4C150v$V12,zsignc)
resid4C150v$c44 <- sapply(X = resid4C150v$V16,zsignc)



zsig4C150<- resid4C150v[,18:33]

datosLong4C150<-pivot_longer(zsig4C150, cols = 1:16,values_to = "MEASURE", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "CATEGORY"))

datosLong4C150$N <- gl(1,80000, labels = "20")
datosLong4C150$CATEGORY<- rep(nombres4,5000)

tabla4C150<-table(datosLong4C150$CATEGORY,datosLong4C150$MEASURE)
tabla4C150p<- as.data.frame(prop.table(tabla4C150,1))
colnames(tabla4C150p)<- c("PATTERN", "DEPENDENCY", "P")
tabla4C150p$CATEGORY<- rep(4,48)
tabla4C150p$MATRIX<- rep("C",48)
tabla4C150p$N<- rep(150,48)


tabla4_150p<- rbind(tabla4A150p,tabla4B150p,tabla4C150p)


# Data summary table


tabla4_exp<- rbind(tabla4_20p,tabla4_50p,tabla4_100p,tabla4_150p)


# Activation patterns

tab4bc1<- subset(tabla4_exp,PATTERN == "BC" & DEPENDENCY == "Act" )
tab4dc1<- subset(tabla4_exp,PATTERN == "DC" & DEPENDENCY == "Act" )
tab4_act_exp<- rbind(tab4bc1,tab4dc1)
tab4_act_exp$MODEL<- gl(1,24, labels = "EXPL")

write.xlsx(tab4_act_exp, "Tabla4C_act_exp.xlsx")


# Inhibition patterns


tab4ba1<- subset(tabla4_exp,PATTERN == "BA" & DEPENDENCY == "Inh" )
tab4bb1<- subset(tabla4_exp,PATTERN == "BB" & DEPENDENCY == "Inh" )
tab4bd1<- subset(tabla4_exp,PATTERN == "BD" & DEPENDENCY == "Inh" )


tab4da1<- subset(tabla4_exp,PATTERN == "DA" & DEPENDENCY == "Inh" )
tab4db1<- subset(tabla4_exp,PATTERN == "DB" & DEPENDENCY == "Inh" )
tab4dd1<- subset(tabla4_exp,PATTERN == "DD" & DEPENDENCY == "Inh" )


tab4_inh_exp<- rbind(tab4ba1, tab4bb1,tab4bd1, tab4da1,tab4db1,tab4dd1)
tab4_inh_exp$MODEL<- gl(1,72, labels = "EXPLORATORY")

write.xlsx(tab4_inh_exp, "Tabla4C_inh_exp.xlsx")



# Spurious patterns

tab4ac1<- subset(tabla4_exp,PATTERN == "AC" & DEPENDENCY == "Inh" )
tab4cc1<- subset(tabla4_exp,PATTERN == "CC" & DEPENDENCY == "Inh" )


tab4_sp_exp<- rbind(tab4ac1, tab4cc1)
tab4_sp_exp$MODEL<- gl(1,24, labels = "EXPL")

write.xlsx(tab4_sp_exp, "Tabla4C_sp_exp.xlsx")


end <- Sys.time()
print(end-start)