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

## Cálculo de la matriz de transición

transi<- function(x) {
  TransMatrix(as.numeric(x),order= 1,probs= F)
}



## Cálculo de residuales

resid2<- function(x){
  RS <- rowSums(x)
  CS <- colSums(x) 
  GT <- sum(x) 
  CST<- chisq.test(x)
  ASR <- (CST$observed - CST$expected) / sqrt(CST$expected * ((1 - RS / GT) %*% t(1 - CS / GT))) 
  ASR
}


## Significación de los residuales

zsign<-function(x){
  ifelse(abs(x) > 1.96,1,0) 
}

zsignb<-function(x){
  ifelse(abs(x) > 2.73,1,0) 
}


zsignc<-function(x){
  ifelse(x > 1.96, "Act",ifelse(x < -1.96, "Inh", "Ind"))
}



# Matrices con 6 estados


## Matriz nula N

m6A <-matrix(c(1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6,
               1/6, 1/6, 1/6,1/6,1/6, 1/6), 
             6, 6, 
             byrow = T)

## Matriz B con p = .75, 1 casilla significativa (2,4)


m6B <-matrix(c(1/6,  1/6,  1/6,  1/6,   1/6,  1/6,
               .05,  .05,  .75,  .05,   .05, .05,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6), 
             6, 6, 
             byrow = T)


## Matriz C con p = .75, 2 casillas significativas (2,4) y (4,4)

m6C <-matrix(c(1/6,  1/6,  1/6,  1/6,   1/6,  1/6,
               .05,  .05,  .75,  .05,   .05, .05,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               .05,  .05,  .75,  .05,   .05, .05,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6,
               1/6,  1/6,  1/6,  1/6,   1/6, 1/6), 
             6, 6, 
             byrow = T)

## nombres

nombres6<- c("AA","AB","AC","AD","AE","AF",
            "BA","BB","BC","BD","BE","BF",
            "CA","CB","CC","CD","CE","CF",
            "DA","DB","DC","DD","DE","DF",
            "EA","EB","EC","ED","EE","EF",
            "FA","FB","FC","FD","FE","FF")


## Secuencias con 6 estados



seq6A20 <- genMarkov(n = 15000, transMat = m6A, 
                     chainLen = 20, wide = TRUE)
seq6B20 <- genMarkov(n = 15000, transMat = m6B, 
                     chainLen = 20, wide = TRUE)
seq6C20 <- genMarkov(n = 15000, transMat = m6C, 
                     chainLen = 20, wide = TRUE)


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


seq6A20<- seq6A20[,-1]
seq6B20<- seq6B20[,-1]
seq6C20<- seq6C20[,-1]

seq6A50<- seq6A50[,-1]
seq6B50<- seq6B50[,-1]
seq6C50<- seq6C50[,-1]

seq6A100<- seq6A100[,-1]
seq6B100<- seq6B100[,-1]
seq6C100<- seq6C100[,-1]




# Residuales de la condición de 6 categorías

## N = 20, matriz A

transiz6A20<-alply(seq6A20,1,transi)
resid6A20<- llply(transiz6A20, function(x) resid2(x))

dim6<- ldply(resid6A20,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

resid6A20 <- resid6A20[-falsos]
resid6A20v<- ldply(resid6A20,as.vector)
buenos<- complete.cases(resid6A20v)
resid6A20v<- resid6A20v[buenos,]
resid6A20v<- resid6A20v[1:5000,]


resid6A20v$c11 <- sapply(X = resid6A20v$V1,zsignc)
resid6A20v$c12 <- sapply(X = resid6A20v$V7,zsignc)
resid6A20v$c13 <- sapply(X = resid6A20v$V13,zsignc)
resid6A20v$c14 <- sapply(X = resid6A20v$V19,zsignc)
resid6A20v$c15 <- sapply(X = resid6A20v$V25,zsignc)
resid6A20v$c16 <- sapply(X = resid6A20v$V31,zsignc)


resid6A20v$c21 <- sapply(X = resid6A20v$V2,zsignc)
resid6A20v$c22 <- sapply(X = resid6A20v$V8,zsignc)
resid6A20v$c23 <- sapply(X = resid6A20v$V14,zsignc)
resid6A20v$c24 <- sapply(X = resid6A20v$V20,zsignc)
resid6A20v$c25 <- sapply(X = resid6A20v$V26,zsignc)
resid6A20v$c26 <- sapply(X = resid6A20v$V28,zsignc)


resid6A20v$c31 <- sapply(X = resid6A20v$V3,zsignc)
resid6A20v$c32 <- sapply(X = resid6A20v$V9,zsignc)
resid6A20v$c33 <- sapply(X = resid6A20v$V15,zsignc)
resid6A20v$c34 <- sapply(X = resid6A20v$V21,zsignc)
resid6A20v$c35 <- sapply(X = resid6A20v$V27,zsignc)
resid6A20v$c36 <- sapply(X = resid6A20v$V33,zsignc)


resid6A20v$c41 <- sapply(X = resid6A20v$V4,zsignc)
resid6A20v$c42 <- sapply(X = resid6A20v$V10,zsignc)
resid6A20v$c43 <- sapply(X = resid6A20v$V16,zsignc)
resid6A20v$c44 <- sapply(X = resid6A20v$V22,zsignc)
resid6A20v$c45 <- sapply(X = resid6A20v$V28,zsignc)
resid6A20v$c46 <- sapply(X = resid6A20v$V34,zsignc)


resid6A20v$c51 <- sapply(X = resid6A20v$V5,zsignc)
resid6A20v$c52 <- sapply(X = resid6A20v$V11,zsignc)
resid6A20v$c53 <- sapply(X = resid6A20v$V17,zsignc)
resid6A20v$c54 <- sapply(X = resid6A20v$V23,zsignc)
resid6A20v$c55 <- sapply(X = resid6A20v$V29,zsignc)
resid6A20v$c56 <- sapply(X = resid6A20v$V32,zsignc)



resid6A20v$c61 <- sapply(X = resid6A20v$V6,zsignc)
resid6A20v$c62 <- sapply(X = resid6A20v$V12,zsignc)
resid6A20v$c63 <- sapply(X = resid6A20v$V18,zsignc)
resid6A20v$c64 <- sapply(X = resid6A20v$V24,zsignc)
resid6A20v$c65 <- sapply(X = resid6A20v$V30,zsignc)
resid6A20v$c66 <- sapply(X = resid6A20v$V36,zsignc)

# Cálculo de los patrones  de activación, inhibición e independencia 

zsig6A20<- resid6A20v[,38:73]
 
datosLong6A20<-pivot_longer(zsig6A20, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong6A20$N <- gl(1,180000, labels = "20")
datosLong6A20$Categoria<- rep(nombres6,5000)
 
tabla6A20<-table(datosLong6A20$Categoria,datosLong6A20$Medida)
tabla6A20p<- as.data.frame(prop.table(tabla6A20,1))
colnames(tabla6A20p)<- c("Patrón", "Dependencia", "P")
tabla6A20p$Cat<- rep(6,108)
tabla6A20p$Matriz<- rep("A",108)
tabla6A20p$N<- rep(20,108)


## N = 20, matriz B

transiz6B20<-alply(seq6B20,1,transi)
resid6B20<- llply(transiz6B20, function(x) resid2(x))

dim6<- ldply(resid6B20,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

resid6B20 <- resid6B20[-falsos]
resid6B20v<- ldply(resid6B20,as.vector)
buenos<- complete.cases(resid6B20v)
resid6B20v<- resid6B20v[buenos,]
resid6B20v<- resid6B20v[1:5000,]


resid6B20v$c11 <- sapply(X = resid6B20v$V1,zsignc)
resid6B20v$c12 <- sapply(X = resid6B20v$V7,zsignc)
resid6B20v$c13 <- sapply(X = resid6B20v$V13,zsignc)
resid6B20v$c14 <- sapply(X = resid6B20v$V19,zsignc)
resid6B20v$c15 <- sapply(X = resid6B20v$V25,zsignc)
resid6B20v$c16 <- sapply(X = resid6B20v$V31,zsignc)


resid6B20v$c21 <- sapply(X = resid6B20v$V2,zsignc)
resid6B20v$c22 <- sapply(X = resid6B20v$V8,zsignc)
resid6B20v$c23 <- sapply(X = resid6B20v$V14,zsignc)
resid6B20v$c24 <- sapply(X = resid6B20v$V20,zsignc)
resid6B20v$c25 <- sapply(X = resid6B20v$V26,zsignc)
resid6B20v$c26 <- sapply(X = resid6B20v$V28,zsignc)


resid6B20v$c31 <- sapply(X = resid6B20v$V3,zsignc)
resid6B20v$c32 <- sapply(X = resid6B20v$V9,zsignc)
resid6B20v$c33 <- sapply(X = resid6B20v$V15,zsignc)
resid6B20v$c34 <- sapply(X = resid6B20v$V21,zsignc)
resid6B20v$c35 <- sapply(X = resid6B20v$V27,zsignc)
resid6B20v$c36 <- sapply(X = resid6B20v$V33,zsignc)


resid6B20v$c41 <- sapply(X = resid6B20v$V4,zsignc)
resid6B20v$c42 <- sapply(X = resid6B20v$V10,zsignc)
resid6B20v$c43 <- sapply(X = resid6B20v$V16,zsignc)
resid6B20v$c44 <- sapply(X = resid6B20v$V22,zsignc)
resid6B20v$c45 <- sapply(X = resid6B20v$V28,zsignc)
resid6B20v$c46 <- sapply(X = resid6B20v$V34,zsignc)


resid6B20v$c51 <- sapply(X = resid6B20v$V5,zsignc)
resid6B20v$c52 <- sapply(X = resid6B20v$V11,zsignc)
resid6B20v$c53 <- sapply(X = resid6B20v$V17,zsignc)
resid6B20v$c54 <- sapply(X = resid6B20v$V23,zsignc)
resid6B20v$c55 <- sapply(X = resid6B20v$V29,zsignc)
resid6B20v$c56 <- sapply(X = resid6B20v$V32,zsignc)



resid6B20v$c61 <- sapply(X = resid6B20v$V6,zsignc)
resid6B20v$c62 <- sapply(X = resid6B20v$V12,zsignc)
resid6B20v$c63 <- sapply(X = resid6B20v$V18,zsignc)
resid6B20v$c64 <- sapply(X = resid6B20v$V24,zsignc)
resid6B20v$c65 <- sapply(X = resid6B20v$V30,zsignc)
resid6B20v$c66 <- sapply(X = resid6B20v$V36,zsignc)



zsig6B20<- resid6B20v[,38:73]

datosLong6B20<-pivot_longer(zsig6B20, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6B20$N <- gl(1,180000, labels = "20")
datosLong6B20$Categoria<- rep(nombres6,5000)

tabla6B20<-table(datosLong6B20$Categoria,datosLong6B20$Medida)
tabla6B20p<- as.data.frame(prop.table(tabla6B20,1))
colnames(tabla6B20p)<- c("Patrón", "Dependencia", "P")
tabla6B20p$Cat<- rep(6,108)
tabla6B20p$Matriz<- rep("B",108)
tabla6B20p$N<- rep(20,108)



## N = 20, matriz C

transiz6C20<-alply(seq6C20,1,transi)
resid6C20<- llply(transiz6C20, function(x) resid2(x))

dim6<- ldply(resid6C20,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

resid6C20 <- resid6C20[-falsos]
resid6C20v<- ldply(resid6C20,as.vector)
buenos<- complete.cases(resid6C20v)
resid6C20v<- resid6C20v[buenos,]
resid6C20v<- resid6C20v[1:5000,]


resid6C20v$c11 <- sapply(X = resid6C20v$V1,zsignc)
resid6C20v$c12 <- sapply(X = resid6C20v$V7,zsignc)
resid6C20v$c13 <- sapply(X = resid6C20v$V13,zsignc)
resid6C20v$c14 <- sapply(X = resid6C20v$V19,zsignc)
resid6C20v$c15 <- sapply(X = resid6C20v$V25,zsignc)
resid6C20v$c16 <- sapply(X = resid6C20v$V31,zsignc)


resid6C20v$c21 <- sapply(X = resid6C20v$V2,zsignc)
resid6C20v$c22 <- sapply(X = resid6C20v$V8,zsignc)
resid6C20v$c23 <- sapply(X = resid6C20v$V14,zsignc)
resid6C20v$c24 <- sapply(X = resid6C20v$V20,zsignc)
resid6C20v$c25 <- sapply(X = resid6C20v$V26,zsignc)
resid6C20v$c26 <- sapply(X = resid6C20v$V28,zsignc)


resid6C20v$c31 <- sapply(X = resid6C20v$V3,zsignc)
resid6C20v$c32 <- sapply(X = resid6C20v$V9,zsignc)
resid6C20v$c33 <- sapply(X = resid6C20v$V15,zsignc)
resid6C20v$c34 <- sapply(X = resid6C20v$V21,zsignc)
resid6C20v$c35 <- sapply(X = resid6C20v$V27,zsignc)
resid6C20v$c36 <- sapply(X = resid6C20v$V33,zsignc)


resid6C20v$c41 <- sapply(X = resid6C20v$V4,zsignc)
resid6C20v$c42 <- sapply(X = resid6C20v$V10,zsignc)
resid6C20v$c43 <- sapply(X = resid6C20v$V16,zsignc)
resid6C20v$c44 <- sapply(X = resid6C20v$V22,zsignc)
resid6C20v$c45 <- sapply(X = resid6C20v$V28,zsignc)
resid6C20v$c46 <- sapply(X = resid6C20v$V34,zsignc)


resid6C20v$c51 <- sapply(X = resid6C20v$V5,zsignc)
resid6C20v$c52 <- sapply(X = resid6C20v$V11,zsignc)
resid6C20v$c53 <- sapply(X = resid6C20v$V17,zsignc)
resid6C20v$c54 <- sapply(X = resid6C20v$V23,zsignc)
resid6C20v$c55 <- sapply(X = resid6C20v$V29,zsignc)
resid6C20v$c56 <- sapply(X = resid6C20v$V32,zsignc)



resid6C20v$c61 <- sapply(X = resid6C20v$V6,zsignc)
resid6C20v$c62 <- sapply(X = resid6C20v$V12,zsignc)
resid6C20v$c63 <- sapply(X = resid6C20v$V18,zsignc)
resid6C20v$c64 <- sapply(X = resid6C20v$V24,zsignc)
resid6C20v$c65 <- sapply(X = resid6C20v$V30,zsignc)
resid6C20v$c66 <- sapply(X = resid6C20v$V36,zsignc)



zsig6C20<- resid6C20v[,38:73]

datosLong6C20<-pivot_longer(zsig6C20, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6C20$N <- gl(1,180000, labels = "20")
datosLong6C20$Categoria<- rep(nombres6,5000)

tabla6C20<-table(datosLong6C20$Categoria,datosLong6C20$Medida)
tabla6C20p<- as.data.frame(prop.table(tabla6C20,1))
colnames(tabla6C20p)<- c("Patrón", "Dependencia", "P")
tabla6C20p$Cat<- rep(6,108)
tabla6C20p$Matriz<- rep("C",108)
tabla6C20p$N<- rep(20,108)

######


## N = 50, matriz A

transiz6A50<-alply(seq6A50,1,transi)
resid6A50<- llply(transiz6A50, function(x) resid2(x))

dim6<- ldply(resid6A50,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

resid6A50 <- resid6A50[-falsos]
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

datosLong6A50<-pivot_longer(zsig6A50, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6A50$N <- gl(1,180000, labels = "50")

datosLong6A50$Categoria<- rep(nombres6,5000)

tabla6A50<-table(datosLong6A50$Categoria,datosLong6A50$Medida)
tabla6A50p<- as.data.frame(prop.table(tabla6A50,1))
colnames(tabla6A50p)<- c("Patrón", "Dependencia", "P")
tabla6A50p$Cat<- rep(6,108)
tabla6A50p$Matriz<- rep("A",108)
tabla6A50p$N<- rep(50,108)




## N = 50, matriz B

transiz6B50<-alply(seq6B50,1,transi)
resid6B50<- llply(transiz6B50, function(x) resid2(x))

dim6<- ldply(resid6B50,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

resid6B50 <- resid6B50[-falsos]
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

datosLong6B50<-pivot_longer(zsig6B50, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6B50$N <- gl(1,180000, labels = "50")

datosLong6B50$Categoria<- rep(nombres6,5000)

tabla6B50<-table(datosLong6B50$Categoria,datosLong6B50$Medida)
tabla6B50p<- as.data.frame(prop.table(tabla6B50,1))
colnames(tabla6B50p)<- c("Patrón", "Dependencia", "P")
tabla6B50p$Cat<- rep(6,108)
tabla6B50p$Matriz<- rep("B",108)
tabla6B50p$N<- rep(50,108)



## N = 50, matriz C

transiz6C50<-alply(seq6C50,1,transi)
resid6C50<- llply(transiz6C50, function(x) resid2(x))

dim6<- ldply(resid6C50,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

resid6C50 <- resid6C50[-falsos]
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

datosLong6C50<-pivot_longer(zsig6C50, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6C50$N <- gl(1,180000, labels = "50")

datosLong6C50$Categoria<- rep(nombres6,5000)

tabla6C50<-table(datosLong6C50$Categoria,datosLong6C50$Medida)
tabla6C50p<- as.data.frame(prop.table(tabla6C50,1))
colnames(tabla6C50p)<- c("Patrón", "Dependencia", "P")
tabla6C50p$Cat<- rep(6,108)
tabla6C50p$Matriz<- rep("C",108)
tabla6C50p$N<- rep(50,108)



########



## N = 100, matriz A

transiz6A100<-alply(seq6A100,1,transi)
resid6A100<- llply(transiz6A100, function(x) resid2(x))

dim6<- ldply(resid6A100,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

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

datosLong6A100<-pivot_longer(zsig6A100, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6A100$N <- gl(1,180000, labels = "100")

datosLong6A100$Categoria<- rep(nombres6,5000)

tabla6A100<-table(datosLong6A100$Categoria,datosLong6A100$Medida)
tabla6A100p<- as.data.frame(prop.table(tabla6A100,1))
colnames(tabla6A100p)<- c("Patrón", "Dependencia", "P")
tabla6A100p$Cat<- rep(6,108)
tabla6A100p$Matriz<- rep("A",108)
tabla6A100p$N<- rep(100,108)




## N = 100, matriz B

transiz6B100<-alply(seq6B100,1,transi)
resid6B100<- llply(transiz6B100, function(x) resid2(x))

dim6<- ldply(resid6B100,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

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

datosLong6B100<-pivot_longer(zsig6B100, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6B100$N <- gl(1,180000, labels = "100")
datosLong6B100$Categoria<- rep(nombres6,5000)

tabla6B100<-table(datosLong6B100$Categoria,datosLong6B100$Medida)
tabla6B100p<- as.data.frame(prop.table(tabla6B100,1))
colnames(tabla6B100p)<- c("Patrón", "Dependencia", "P")
tabla6B100p$Cat<- rep(6,108)
tabla6B100p$Matriz<- rep("B",108)
tabla6B100p$N<- rep(100,108)


## N = 100, matriz C

transiz6C100<-alply(seq6C100,1,transi)
resid6C100<- llply(transiz6C100, function(x) resid2(x))

dim6<- ldply(resid6C100,dim)
falsos<-which(dim6$V1 != 6 | dim6$V2 !=6)
length(falsos)

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

datosLong6C100<-pivot_longer(zsig6C100, cols = 1:36,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong6C100$N <- gl(1,180000, labels = "100")

datosLong6C100$Categoria<- rep(nombres6,5000)

tabla6C100<-table(datosLong6C100$Categoria,datosLong6C100$Medida)
tabla6C100p<- as.data.frame(prop.table(tabla6C100,1))
colnames(tabla6C100p)<- c("Patrón", "Dependencia", "P")
tabla6C100p$Cat<- rep(6,108)
tabla6C100p$Matriz<- rep("C",108)
tabla6C100p$N<- rep(100,108)



# Tabla resumen de los datos


tabla6<- rbind(tabla6A20p,tabla6A50p,tabla6A100p,tabla6B20p,tabla6B50p,tabla6B100p, 
               tabla6C20p,tabla6C50p,tabla6C100p)

patron_BC6<- tabla6[tabla6$Patrón== "BC",]


# Representación gráfica

plot6<-ggplot(patron_BC6,aes(x = factor(N), y = P,group = Dependencia,shape = Dependencia))+
  geom_line(aes(linetype = Dependencia),size = 1)+
  geom_point(size =2)+
  facet_grid(~Matriz)+
  theme_gray()

print(plot6)



tab6c<- pivot_wider(tabla6, names_from = c(Matriz,Dependencia), values_from = P)
write.csv2(tab6c, "Tabla6.csv")


