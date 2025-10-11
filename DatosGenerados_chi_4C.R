library(SMM)
library(smmR)
library(simstudy)
library(markovchain)
library(march)
library(MCTM)
library(plyr)
library(DescTools)
library(DescToolsAddIns)
library(ggplot2)
library(plotrix)
library(plotly)
library(psych)
library(dplyr)
library(tidyr)

## MCTM library is removed from R. 
## Formerly available versions can be obtained from the archive.

set.seed(123456)

# Funciones 
 
## Cálculo de las frecuencias de transición

transi<- function(x) {
  TransMatrix(as.numeric(x),order= 1,probs= F)
}


# Cálculo de residuales

resid2<- function(x){
  RS <- rowSums(x)
  CS <- colSums(x) 
  GT <- sum(x) 
  CST<- chisq.test(x)
  ASR <- (CST$observed - CST$expected) / sqrt(CST$expected * ((1 - RS / GT) %*% t(1 - CS / GT))) 
  ASR
}


# Significación de los residuales

zsign<-function(x){
  ifelse(abs(x) > 1.96,1,0) 
}

zsignb<-function(x){
  ifelse(abs(x) > 2.73,1,0) 
}


zsignc<-function(x){
  ifelse(x > 1.96, "Act",ifelse(x < -1.96, "Inh", "Ind"))
}



# Matrices con 4 estados


## Matriz nula   A


m4A <- matrix(c(1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4), 4, 4, 
              byrow = T)

## Matriz B p =3/4 estados = 4, 1 significativo (BC)

m4B <- matrix(c(1/4, 1/4, 1/4,1/4,
                1/12, 1/12, 9/12,1/12,
                1/4, 1/4, 1/4,1/4,
                1/4, 1/4, 1/4,1/4), 4, 4, 
              byrow = T)


## Matriz C p =3/4 estados = 4, 2 significativos (BC y DC)

m4C <- matrix(c( 1/4, 1/4, 1/4,1/4,
                 1/12, 1/12, 9/12,1/12,
                 1/4, 1/4, 1/4,1/4,
                 1/12, 1/12, 9/12,1/12), 4, 4, 
              byrow = T)



## nombres4

nombres4<- c("AA","AB","AC","AD",
            "BA","BB","BC","BD",
            "CA","CB","CC","CD",
            "DA","DB","DC","DD")
            

# Generación de las secuencias

## Secuencias con 4 estados

seq4A20 <- genMarkov(n = 25000, transMat = m4A, 
                     chainLen = 20, wide = TRUE)
seq4B20 <- genMarkov(n = 45000, transMat = m4B, 
                     chainLen = 20, wide = TRUE)
seq4C20 <- genMarkov(n = 45000, transMat = m4C, 
                     chainLen = 20, wide = TRUE)


seq4A50 <- genMarkov(n = 25000, transMat = m4A, 
                     chainLen = 50, wide = TRUE)
seq4B50 <- genMarkov(n = 25000, transMat = m4B, 
                     chainLen = 50, wide = TRUE)
seq4C50 <- genMarkov(n = 25000, transMat = m4C, 
                     chainLen = 50, wide = TRUE)


seq4A100 <- genMarkov(n = 25000, transMat = m4A, 
                      chainLen = 100, wide = TRUE)
seq4B100 <- genMarkov(n = 25000, transMat = m4B, 
                      chainLen = 100, wide = TRUE)
seq4C100 <- genMarkov(n = 25000, transMat = m4C, 
                      chainLen = 100, wide = TRUE)


seq4A20<- seq4A20[,-1]
seq4B20<- seq4B20[,-1]
seq4C20<- seq4C20[,-1]

seq4A50<- seq4A50[,-1]
seq4B50<- seq4B50[,-1]
seq4C50<- seq4C50[,-1]

seq4A100<- seq4A100[,-1]
seq4B100<- seq4B100[,-1]
seq4C100<- seq4C100[,-1]



# chi cuadrado de la condición de 4 categorías

## N = 20, matriz A

transiz4A20<-alply(seq4A20,1,transi)
dim4<- ldply(transiz4A20,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

transiz4A20 <- transiz4A20[-falsos]
chisq4A20<- llply(transiz4A20,function(x) chisq.test(x)$statistic)
chi4A20<- llply(chisq4A20, as.vector)
chi24A20<- unlist(chi4A20)
buenos<- which(chi4A20 <= qchisq(.95,6))
quantile(chi24A20[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)

transiz4A20<- transiz4A20[buenos]
resid4A20<- llply(transiz4A20, function(x) resid2(x))

resid4A20<- llply(transiz4A20, function(x) resid2(x))

resid4A20v<- ldply(resid4A20,as.vector)
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

datosLong4A20<-pivot_longer(zsig4A20, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong4A20$N <- gl(1,80000, labels = "20")

datosLong4A20$Categoria<- rep(nombres4,5000)

tabla4A20<-table(datosLong4A20$Categoria,datosLong4A20$Medida)
tabla4A20p<- as.data.frame(prop.table(tabla4A20,1))
colnames(tabla4A20p)<- c("Patrón", "Dependencia", "P")
tabla4A20p$Cat<- rep(4,48)
tabla4A20p$Matriz<- rep("A",48)
tabla4A20p$N<- rep(20,48)




# N = 20 Matriz B

transiz4B20<-alply(seq4B20,1,transi)
dim4<- ldply(transiz4B20,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

transiz4B20 <- transiz4B20[-falsos]
chisq4B20<- llply(transiz4B20,function(x) chisq.test(x)$statistic)
chi4B20<- llply(chisq4B20, as.vector)
chi24B20<- unlist(chi4B20)
buenos<- which(chi4B20 > qchisq(.95,6))
quantile(chi24B20[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)

transiz4B20<- transiz4B20[buenos]
resid4B20<- llply(transiz4B20, function(x) resid2(x))

resid4B20<- llply(transiz4B20, function(x) resid2(x))

resid4B20v<- ldply(resid4B20,as.vector)
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

datosLong4B20<-pivot_longer(zsig4B20, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong4B20$N <- gl(1,80000, labels = "20")

datosLong4B20$Categoria<- rep(nombres4,5000)

tabla4B20<-table(datosLong4B20$Categoria,datosLong4B20$Medida)
tabla4B20p<- as.data.frame(prop.table(tabla4B20,1))
colnames(tabla4B20p)<- c("Patrón", "Dependencia", "P")
tabla4B20p$Cat<- rep(4,48)
tabla4B20p$Matriz<- rep("B",48)
tabla4B20p$N<- rep(20,48)



# N = 20 Matriz C

transiz4C20<-alply(seq4C20,1,transi)
dim4<- ldply(transiz4C20,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

transiz4C20 <- transiz4C20[-falsos]
chisq4C20<- llply(transiz4C20,function(x) chisq.test(x)$statistic)
chi4C20<- llply(chisq4C20, as.vector)
chi24C20<- unlist(chi4C20)
buenos<- which(chi4C20 > qchisq(.95,6))
length(buenos)
quantile(chi24C20[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)

transiz4C20<- transiz4C20[buenos]
resid4C20<- llply(transiz4C20, function(x) resid2(x))

resid4C20<- llply(transiz4C20, function(x) resid2(x))

resid4C20v<- ldply(resid4C20,as.vector)
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

datosLong4C20<-pivot_longer(zsig4C20, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong4C20$N <- gl(1,80000, labels = "20")

datosLong4C20$Categoria<- rep(nombres4,5000)

tabla4C20<-table(datosLong4C20$Categoria,datosLong4C20$Medida)
tabla4C20p<- as.data.frame(prop.table(tabla4C20,1))
colnames(tabla4C20p)<- c("Patrón", "Dependencia", "P")
tabla4C20p$Cat<- rep(4,48)
tabla4C20p$Matriz<- rep("C",48)
tabla4C20p$N<- rep(20,48)


tabla4_20p<- rbind(tabla4A20p,tabla4B20p,tabla4C20p)

# N = 50 Matriz A

transiz4A50<-alply(seq4A50,1,transi)
dim4<- ldply(transiz4A50,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

chisq4A50<- llply(transiz4A50,function(x) chisq.test(x)$statistic)
chi4A50<- llply(chisq4A50, as.vector)
chi24A50<- unlist(chi4A50)
buenos<- which(chi4A50 <= qchisq(.95,6))
length(buenos)
quantile(chi24A50[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)
resid4A50<- llply(transiz4A50, function(x) resid2(x))
resid4A50<- llply(transiz4A50, function(x) resid2(x))
resid4A50v<- ldply(resid4A50,as.vector)
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

datosLong4A50<-pivot_longer(zsig4A50, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong4A50$N <- gl(1,80000, labels = "20")

datosLong4A50$Categoria<- rep(nombres4,5000)

tabla4A50<-table(datosLong4A50$Categoria,datosLong4A50$Medida)
tabla4A50p<- as.data.frame(prop.table(tabla4A50,1))
colnames(tabla4A50p)<- c("Patrón", "Dependencia", "P")
tabla4A50p$Cat<- rep(4,48)
tabla4A50p$Matriz<- rep("A",48)
tabla4A50p$N<- rep(50,48)


# N = 50 Matriz B

transiz4B50<-alply(seq4B50,1,transi)
dim4<- ldply(transiz4B50,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

chisq4B50<- llply(transiz4B50,function(x) chisq.test(x)$statistic)
chi4B50<- llply(chisq4B50, as.vector)
chi24B50<- unlist(chi4B50)
buenos<- which(chi4B50 > qchisq(.95,6))
length(buenos)
quantile(chi24B50[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)
resid4B50<- llply(transiz4B50, function(x) resid2(x))
resid4B50<- llply(transiz4B50, function(x) resid2(x))
resid4B50v<- ldply(resid4B50,as.vector)
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

datosLong4B50<-pivot_longer(zsig4B50, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong4B50$N <- gl(1,80000, labels = "20")
datosLong4B50$Categoria<- rep(nombres4,5000)

tabla4B50<-table(datosLong4B50$Categoria,datosLong4B50$Medida)
tabla4B50p<- as.data.frame(prop.table(tabla4B50,1))
colnames(tabla4B50p)<- c("Patrón", "Dependencia", "P")
tabla4B50p$Cat<- rep(4,48)
tabla4B50p$Matriz<- rep("B",48)
tabla4B50p$N<- rep(50,48)




# N = 50 Matriz C

transiz4C50<-alply(seq4C50,1,transi)
dim4<- ldply(transiz4C50,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

transiz4C50 <- transiz4C50[-falsos]
chisq4C50<- llply(transiz4C50,function(x) chisq.test(x)$statistic)
chi4C50<- llply(chisq4C50, as.vector)
chi24C50<- unlist(chi4C50)
buenos<- which(chi4C50 > qchisq(.95,6))
length(buenos)
quantile(chi24C50[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)
resid4C50<- llply(transiz4C50, function(x) resid2(x))
resid4C50<- llply(transiz4C50, function(x) resid2(x))
resid4C50v<- ldply(resid4C50,as.vector)
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

datosLong4C50<-pivot_longer(zsig4C50, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong4C50$N <- gl(1,80000, labels = "20")
datosLong4C50$Categoria<- rep(nombres4,5000)

tabla4C50<-table(datosLong4C50$Categoria,datosLong4C50$Medida)
tabla4C50p<- as.data.frame(prop.table(tabla4C50,1))
colnames(tabla4C50p)<- c("Patrón", "Dependencia", "P")
tabla4C50p$Cat<- rep(4,48)
tabla4C50p$Matriz<- rep("C",48)
tabla4C50p$N<- rep(50,48)

tabla4_50p<- rbind(tabla4A50p,tabla4B50p,tabla4C50p)




# N = 100 Matriz A

transiz4A100<-alply(seq4A100,1,transi)
dim4<- ldply(transiz4A100,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

chisq4A100<- llply(transiz4A100,function(x) chisq.test(x)$statistic)
chi4A100<- llply(chisq4A100, as.vector)
chi24A100<- unlist(chi4A100)
buenos<- which(chi4A100 < qchisq(.95,6))
length(buenos)
quantile(chi24A100[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)
resid4A100<- llply(transiz4A100, function(x) resid2(x))
resid4A100<- llply(transiz4A100, function(x) resid2(x))
resid4A100v<- ldply(resid4A100,as.vector)
resid4A100v<- resid4A50v[1:5000,]

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

datosLong4A100<-pivot_longer(zsig4A100, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong4A100$N <- gl(1,80000, labels = "20")
datosLong4A100$Categoria<- rep(nombres4,5000)

tabla4A100<-table(datosLong4A100$Categoria,datosLong4A100$Medida)
tabla4A100p<- as.data.frame(prop.table(tabla4A100,1))
colnames(tabla4A100p)<- c("Patrón", "Dependencia", "P")
tabla4A100p$Cat<- rep(4,48)
tabla4A100p$Matriz<- rep("A",48)
tabla4A100p$N<- rep(100,48)



# N = 100 Matriz B

transiz4B100<-alply(seq4B100,1,transi)
dim4<- ldply(transiz4B100,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

chisq4B100<- llply(transiz4B100,function(x) chisq.test(x)$statistic)
chi4B100<- llply(chisq4B100, as.vector)
chi24B100<- unlist(chi4B100)
buenos<- which(chi4B100 > qchisq(.95,6))
length(buenos)
quantile(chi24B100[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)
resid4B100<- llply(transiz4B100, function(x) resid2(x))
resid4B100<- llply(transiz4B100, function(x) resid2(x))
resid4B100v<- ldply(resid4B100,as.vector)
resid4B100v<- resid4B50v[1:5000,]



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

datosLong4B100<-pivot_longer(zsig4B100, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))
datosLong4B100$N <- gl(1,80000, labels = "20")
datosLong4B100$Categoria<- rep(nombres4,5000)

tabla4B100<-table(datosLong4B100$Categoria,datosLong4B100$Medida)
tabla4B100p<- as.data.frame(prop.table(tabla4B100,1))
colnames(tabla4B100p)<- c("Patrón", "Dependencia", "P")
tabla4B100p$Cat<- rep(4,48)
tabla4B100p$Matriz<- rep("B",48)
tabla4B100p$N<- rep(100,48)


# N = 100 Matriz C

transiz4C100<-alply(seq4C100,1,transi)
dim4<- ldply(transiz4C100,dim)
falsos<-which(dim4$V1 != 4 | dim4$V2 !=4)
length(falsos)

chisq4C100<- llply(transiz4C100,function(x) chisq.test(x)$statistic)
chi4C100<- llply(chisq4C100, as.vector)
chi24C100<- unlist(chi4C100)
buenos<- which(chi4C100 > qchisq(.95,6))
length(buenos)
quantile(chi24C100[buenos],probs=c(.05,.25,.50, .75,.95), na.rm =T)
resid4C100<- llply(transiz4C100, function(x) resid2(x))
resid4C100<- llply(transiz4C100, function(x) resid2(x))
resid4C100v<- ldply(resid4C100,as.vector)
resid4C100v<- resid4C50v[1:5000,]



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

datosLong4C100<-pivot_longer(zsig4C100, cols = 1:16,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong4C100$N <- gl(1,80000, labels = "20")
datosLong4C100$Categoria<- rep(nombres4,5000)

tabla4C100<-table(datosLong4C100$Categoria,datosLong4C100$Medida)
tabla4C100p<- as.data.frame(prop.table(tabla4C100,1))
colnames(tabla4C100p)<- c("Patrón", "Dependencia", "P")
tabla4C100p$Cat<- rep(4,48)
tabla4C100p$Matriz<- rep("C",48)
tabla4C100p$N<- rep(100,48)


tabla4_100p<- rbind(tabla4A100p,tabla4B100p,tabla4C100p)


# Tabla resumen de los datos


tabla4<- rbind(tabla4_20p,tabla4_50p,tabla4_100p)

patron_BC4<- tabla4[tabla4$Patrón== "BC",]

# Representación gráfica

plot4<-ggplot(patron_BC4,aes(x = factor(N), y = P,group = Dependencia,shape = Dependencia))+
  geom_line(aes(linetype = Dependencia),size = 1)+
  geom_point(size =2)+
  facet_grid(~Matriz)+
  theme_gray()

print(plot4)


tab4b<- pivot_wider(tabla4, names_from = c(Matriz,Dependencia), values_from = P)
write.csv2(tab4b, "Tabla4_chi2.csv")

