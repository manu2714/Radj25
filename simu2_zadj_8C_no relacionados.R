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



# Definición  de las matrices con 8 estados


# Matrices con 8 estados 

# Matriz A casillas significativas (1,2) y (3,4)

m82A <- matrix(c(1/28, 6/8, 1/28,1/28,1/28, 1/28, 1/28,1/28,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 1/28,6/8,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
                8, 8, byrow = T)

# Matriz B con p = 6/8, 4 casillas significativas  (1,2), (3,4), (5,6) y (7,8)


m82B <- matrix(c(1/28, 6/8, 1/28,1/28,1/28, 1/28, 1/28,1/28,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                 1/28, 1/28, 1/28,6/8,1/28, 1/28, 1/28,1/28,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                 1/28, 1/28, 1/28,1/28,1/28, 6/8, 1/28,1/28,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                 1/28, 1/28, 1/28,1/28,1/28, 1/28, 1/28,6/8,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
               8, 8, byrow = T)


# Matriz C con p = 6/8, 6 casillas significativas (1,2), (3,4), (5,6), (7,8), (2,5) y (4,8)


m82C <- matrix(c(1/28, 6/8, 1/28,1/28,1/28, 1/28, 1/28,1/28,
                 1/28, 1/28, 1/28,1/28,6/8, 1/28, 1/28,1/28,
                 1/28, 1/28, 1/28,6/8,1/28, 1/28, 1/28,1/28,
                 1/28, 1/28, 1/28,1/28,1/28, 1/28, 1/28,6/8,
                 1/28, 1/28, 1/28,1/28,1/28, 6/8, 1/28,1/28,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                 1/28, 1/28, 1/28,1/28,1/28, 1/28, 1/28,6/8,
                 1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
                 8, 8, byrow = T)



# Condición con 8 estados

## Matriz 2A

seq82A100 <- genMarkov(n = 15000, transMat = m82A, 
                     chainLen = 100, wide = TRUE)
seq82B100 <- genMarkov(n = 15000, transMat = m82B, 
                     chainLen = 100, wide = TRUE)
seq82C100 <- genMarkov(n = 15000, transMat = m82C, 
                     chainLen = 100, wide = TRUE)


seq82A150 <- genMarkov(n = 15000, transMat = m82A, 
                     chainLen = 150, wide = TRUE)
seq82B150 <- genMarkov(n = 15000, transMat = m82B, 
                     chainLen = 150, wide = TRUE)
seq82C150 <- genMarkov(n = 15000, transMat = m82C, 
                     chainLen = 150, wide = TRUE)


seq82A200 <- genMarkov(n = 15000, transMat = m82A, 
                      chainLen = 200, wide = TRUE)
seq82B200 <- genMarkov(n = 15000, transMat = m82B, 
                      chainLen = 200, wide = TRUE)
seq82C200 <- genMarkov(n = 15000, transMat = m82C, 
                      chainLen = 200, wide = TRUE)


seq82A100<- seq82A100[,-1]
seq82B100<- seq82B100[,-1]
seq82C100<- seq82C100[,-1]

seq82A150<- seq82A150[,-1]
seq82B150<- seq82B150[,-1]
seq82C150<- seq82C150[,-1]

seq82A200<- seq82A200[,-1]
seq82B200<- seq82B200[,-1]
seq82C200<- seq82C200[,-1]

## nombres8 

nombres8<- c("AA","AB","AC","AD","AE","AF","AG","AH",
            "BA","BB","BC","BD","BE","BF","BG","BH",
            "CA","CB","CC","CD","CE","CF","CG","CH",
            "DA","DB","DC","DD","DE","DF","DG","DH",
            "EA","EB","EC","ED","EE","EF","EG","EH",
            "FA","FB","FC","FD","FE","FF","FG","FH",
            "GA","GB","GC","GD","GE","GF","GG","GH",
            "HA","HB","HC","HD","HE","HF","HG","HH")


# Residuales de la condición de 8 categorías

## N = 100, matriz A2

transiz82A100<-alply(seq82A100,1,transi)
resid82A100<- llply(transiz82A100, function(x) resid2(x))

dim8<- ldply(resid82A100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

#resid82A100 <- resid82A100[-falsos]
resid82A100v<- ldply(resid82A100,as.vector)
buenos<- complete.cases(resid82A100v)
resid82A100v<- resid82A100v[buenos,]
resid82A100v<- resid82A100v[1:5000,]


resid82A100v$c11 <- sapply(X = resid82A100v$V1,zsignc)
resid82A100v$c12 <- sapply(X = resid82A100v$V9,zsignc)
resid82A100v$c13 <- sapply(X = resid82A100v$V17,zsignc)
resid82A100v$c14 <- sapply(X = resid82A100v$V25,zsignc)
resid82A100v$c15 <- sapply(X = resid82A100v$V33,zsignc)
resid82A100v$c16 <- sapply(X = resid82A100v$V41,zsignc)
resid82A100v$c17 <- sapply(X = resid82A100v$V49,zsignc)
resid82A100v$c18 <- sapply(X = resid82A100v$V57,zsignc)


resid82A100v$c21 <- sapply(X = resid82A100v$V2,zsignc)
resid82A100v$c22 <- sapply(X = resid82A100v$V10,zsignc)
resid82A100v$c23 <- sapply(X = resid82A100v$V18,zsignc)
resid82A100v$c24 <- sapply(X = resid82A100v$V26,zsignc)
resid82A100v$c25 <- sapply(X = resid82A100v$V34,zsignc)
resid82A100v$c26 <- sapply(X = resid82A100v$V42,zsignc)
resid82A100v$c27 <- sapply(X = resid82A100v$V50,zsignc)
resid82A100v$c28 <- sapply(X = resid82A100v$V58,zsignc)


resid82A100v$c31 <- sapply(X = resid82A100v$V3,zsignc)
resid82A100v$c32 <- sapply(X = resid82A100v$V11,zsignc)
resid82A100v$c33 <- sapply(X = resid82A100v$V19,zsignc)
resid82A100v$c34 <- sapply(X = resid82A100v$V27,zsignc)
resid82A100v$c35 <- sapply(X = resid82A100v$V35,zsignc)
resid82A100v$c36 <- sapply(X = resid82A100v$V43,zsignc)
resid82A100v$c37 <- sapply(X = resid82A100v$V51,zsignc)
resid82A100v$c38 <- sapply(X = resid82A100v$V59,zsignc)


resid82A100v$c41 <- sapply(X = resid82A100v$V4,zsignc)
resid82A100v$c42 <- sapply(X = resid82A100v$V12,zsignc)
resid82A100v$c43 <- sapply(X = resid82A100v$V20,zsignc)
resid82A100v$c44 <- sapply(X = resid82A100v$V28,zsignc)
resid82A100v$c45 <- sapply(X = resid82A100v$V36,zsignc)
resid82A100v$c46 <- sapply(X = resid82A100v$V44,zsignc)
resid82A100v$c47 <- sapply(X = resid82A100v$V52,zsignc)
resid82A100v$c48 <- sapply(X = resid82A100v$V60,zsignc)


resid82A100v$c51 <- sapply(X = resid82A100v$V5,zsignc)
resid82A100v$c52 <- sapply(X = resid82A100v$V13,zsignc)
resid82A100v$c53 <- sapply(X = resid82A100v$V21,zsignc)
resid82A100v$c54 <- sapply(X = resid82A100v$V29,zsignc)
resid82A100v$c55 <- sapply(X = resid82A100v$V37,zsignc)
resid82A100v$c56 <- sapply(X = resid82A100v$V45,zsignc)
resid82A100v$c57 <- sapply(X = resid82A100v$V53,zsignc)
resid82A100v$c58 <- sapply(X = resid82A100v$V61,zsignc)


resid82A100v$c61 <- sapply(X = resid82A100v$V6,zsignc)
resid82A100v$c62 <- sapply(X = resid82A100v$V14,zsignc)
resid82A100v$c63 <- sapply(X = resid82A100v$V22,zsignc)
resid82A100v$c64 <- sapply(X = resid82A100v$V30,zsignc)
resid82A100v$c65 <- sapply(X = resid82A100v$V38,zsignc)
resid82A100v$c66 <- sapply(X = resid82A100v$V46,zsignc)
resid82A100v$c67 <- sapply(X = resid82A100v$V54,zsignc)
resid82A100v$c68 <- sapply(X = resid82A100v$V62,zsignc)


resid82A100v$c71 <- sapply(X = resid82A100v$V7,zsignc)
resid82A100v$c72 <- sapply(X = resid82A100v$V15,zsignc)
resid82A100v$c73 <- sapply(X = resid82A100v$V23,zsignc)
resid82A100v$c74 <- sapply(X = resid82A100v$V31,zsignc)
resid82A100v$c75 <- sapply(X = resid82A100v$V39,zsignc)
resid82A100v$c76 <- sapply(X = resid82A100v$V47,zsignc)
resid82A100v$c77 <- sapply(X = resid82A100v$V55,zsignc)
resid82A100v$c78 <- sapply(X = resid82A100v$V63,zsignc)


resid82A100v$c81 <- sapply(X = resid82A100v$V8,zsignc)
resid82A100v$c82 <- sapply(X = resid82A100v$V16,zsignc)
resid82A100v$c83 <- sapply(X = resid82A100v$V24,zsignc)
resid82A100v$c84 <- sapply(X = resid82A100v$V32,zsignc)
resid82A100v$c85 <- sapply(X = resid82A100v$V40,zsignc)
resid82A100v$c86 <- sapply(X = resid82A100v$V48,zsignc)
resid82A100v$c87 <- sapply(X = resid82A100v$V56,zsignc)
resid82A100v$c88 <- sapply(X = resid82A100v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82A100<- resid82A100v[,66:129]

datosLong82A100<-pivot_longer(zsig82A100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82A100$N <- gl(1,320000, labels = 100)
datosLong82A100$Categoria<- rep(nombres8,5000)

tabla82A100<-table(datosLong82A100$Categoria,datosLong82A100$Medida)
tabla82A100p<- as.data.frame(prop.table(tabla82A100,1))
colnames(tabla82A100p)<- c("Patrón", "Dependencia", "P")
tabla82A100p$Cat<- rep(8,192)
tabla82A100p$Matriz<- rep("A",192)
tabla82A100p$N<- rep(100,192)

####

## N = 150, matriz 2A

transiz82A150<-alply(seq82A150,1,transi)
resid82A150<- llply(transiz82A150, function(x) resid2(x))

dim8<- ldply(resid82A150,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid82A150v<- ldply(resid82A150,as.vector)
buenos<- complete.cases(resid82A150v)
resid82A150v<- resid82A150v[buenos,]
resid82A150v<- resid82A150v[1:5000,]


resid82A150v$c11 <- sapply(X = resid82A150v$V1,zsignc)
resid82A150v$c12 <- sapply(X = resid82A150v$V9,zsignc)
resid82A150v$c13 <- sapply(X = resid82A150v$V17,zsignc)
resid82A150v$c14 <- sapply(X = resid82A150v$V25,zsignc)
resid82A150v$c15 <- sapply(X = resid82A150v$V33,zsignc)
resid82A150v$c16 <- sapply(X = resid82A150v$V41,zsignc)
resid82A150v$c17 <- sapply(X = resid82A150v$V49,zsignc)
resid82A150v$c18 <- sapply(X = resid82A150v$V57,zsignc)


resid82A150v$c21 <- sapply(X = resid82A150v$V2,zsignc)
resid82A150v$c22 <- sapply(X = resid82A150v$V10,zsignc)
resid82A150v$c23 <- sapply(X = resid82A150v$V18,zsignc)
resid82A150v$c24 <- sapply(X = resid82A150v$V26,zsignc)
resid82A150v$c25 <- sapply(X = resid82A150v$V34,zsignc)
resid82A150v$c26 <- sapply(X = resid82A150v$V42,zsignc)
resid82A150v$c27 <- sapply(X = resid82A150v$V50,zsignc)
resid82A150v$c28 <- sapply(X = resid82A150v$V58,zsignc)


resid82A150v$c31 <- sapply(X = resid82A150v$V3,zsignc)
resid82A150v$c32 <- sapply(X = resid82A150v$V11,zsignc)
resid82A150v$c33 <- sapply(X = resid82A150v$V19,zsignc)
resid82A150v$c34 <- sapply(X = resid82A150v$V27,zsignc)
resid82A150v$c35 <- sapply(X = resid82A150v$V35,zsignc)
resid82A150v$c36 <- sapply(X = resid82A150v$V43,zsignc)
resid82A150v$c37 <- sapply(X = resid82A150v$V51,zsignc)
resid82A150v$c38 <- sapply(X = resid82A150v$V59,zsignc)


resid82A150v$c41 <- sapply(X = resid82A150v$V4,zsignc)
resid82A150v$c42 <- sapply(X = resid82A150v$V12,zsignc)
resid82A150v$c43 <- sapply(X = resid82A150v$V20,zsignc)
resid82A150v$c44 <- sapply(X = resid82A150v$V28,zsignc)
resid82A150v$c45 <- sapply(X = resid82A150v$V36,zsignc)
resid82A150v$c46 <- sapply(X = resid82A150v$V44,zsignc)
resid82A150v$c47 <- sapply(X = resid82A150v$V52,zsignc)
resid82A150v$c48 <- sapply(X = resid82A150v$V60,zsignc)


resid82A150v$c51 <- sapply(X = resid82A150v$V5,zsignc)
resid82A150v$c52 <- sapply(X = resid82A150v$V13,zsignc)
resid82A150v$c53 <- sapply(X = resid82A150v$V21,zsignc)
resid82A150v$c54 <- sapply(X = resid82A150v$V29,zsignc)
resid82A150v$c55 <- sapply(X = resid82A150v$V37,zsignc)
resid82A150v$c56 <- sapply(X = resid82A150v$V45,zsignc)
resid82A150v$c57 <- sapply(X = resid82A150v$V53,zsignc)
resid82A150v$c58 <- sapply(X = resid82A150v$V61,zsignc)


resid82A150v$c61 <- sapply(X = resid82A150v$V6,zsignc)
resid82A150v$c62 <- sapply(X = resid82A150v$V14,zsignc)
resid82A150v$c63 <- sapply(X = resid82A150v$V22,zsignc)
resid82A150v$c64 <- sapply(X = resid82A150v$V30,zsignc)
resid82A150v$c65 <- sapply(X = resid82A150v$V38,zsignc)
resid82A150v$c66 <- sapply(X = resid82A150v$V46,zsignc)
resid82A150v$c67 <- sapply(X = resid82A150v$V54,zsignc)
resid82A150v$c68 <- sapply(X = resid82A150v$V62,zsignc)


resid82A150v$c71 <- sapply(X = resid82A150v$V7,zsignc)
resid82A150v$c72 <- sapply(X = resid82A150v$V15,zsignc)
resid82A150v$c73 <- sapply(X = resid82A150v$V23,zsignc)
resid82A150v$c74 <- sapply(X = resid82A150v$V31,zsignc)
resid82A150v$c75 <- sapply(X = resid82A150v$V39,zsignc)
resid82A150v$c76 <- sapply(X = resid82A150v$V47,zsignc)
resid82A150v$c77 <- sapply(X = resid82A150v$V55,zsignc)
resid82A150v$c78 <- sapply(X = resid82A150v$V63,zsignc)


resid82A150v$c81 <- sapply(X = resid82A150v$V8,zsignc)
resid82A150v$c82 <- sapply(X = resid82A150v$V16,zsignc)
resid82A150v$c83 <- sapply(X = resid82A150v$V24,zsignc)
resid82A150v$c84 <- sapply(X = resid82A150v$V32,zsignc)
resid82A150v$c85 <- sapply(X = resid82A150v$V40,zsignc)
resid82A150v$c86 <- sapply(X = resid82A150v$V48,zsignc)
resid82A150v$c87 <- sapply(X = resid82A150v$V56,zsignc)
resid82A150v$c88 <- sapply(X = resid82A150v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82A150<- resid82A150v[,66:129]

datosLong82A150<-pivot_longer(zsig82A150, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82A150$N <- gl(1,320000, labels = 150)
datosLong82A150$Categoria<- rep(nombres8,5000)

tabla82A150<-table(datosLong82A150$Categoria,datosLong82A150$Medida)
tabla82A150p<- as.data.frame(prop.table(tabla82A150,1))
colnames(tabla82A150p)<- c("Patrón", "Dependencia", "P")
tabla82A150p$Cat<- rep(8,192)
tabla82A150p$Matriz<- rep("A",192)
tabla82A150p$N<- rep(150,192)

####

## N = 200, matriz 2A

transiz82A200<-alply(seq82A200,1,transi)
resid82A200<- llply(transiz82A200, function(x) resid2(x))

dim8<- ldply(resid82A200,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid82A200v<- ldply(resid82A200,as.vector)
buenos<- complete.cases(resid82A200v)
resid82A200v<- resid82A200v[buenos,]
resid82A200v<- resid82A200v[1:5000,]


resid82A200v$c11 <- sapply(X = resid82A200v$V1,zsignc)
resid82A200v$c12 <- sapply(X = resid82A200v$V9,zsignc)
resid82A200v$c13 <- sapply(X = resid82A200v$V17,zsignc)
resid82A200v$c14 <- sapply(X = resid82A200v$V25,zsignc)
resid82A200v$c15 <- sapply(X = resid82A200v$V33,zsignc)
resid82A200v$c16 <- sapply(X = resid82A200v$V41,zsignc)
resid82A200v$c17 <- sapply(X = resid82A200v$V49,zsignc)
resid82A200v$c18 <- sapply(X = resid82A200v$V57,zsignc)


resid82A200v$c21 <- sapply(X = resid82A200v$V2,zsignc)
resid82A200v$c22 <- sapply(X = resid82A200v$V10,zsignc)
resid82A200v$c23 <- sapply(X = resid82A200v$V18,zsignc)
resid82A200v$c24 <- sapply(X = resid82A200v$V26,zsignc)
resid82A200v$c25 <- sapply(X = resid82A200v$V34,zsignc)
resid82A200v$c26 <- sapply(X = resid82A200v$V42,zsignc)
resid82A200v$c27 <- sapply(X = resid82A200v$V50,zsignc)
resid82A200v$c28 <- sapply(X = resid82A200v$V58,zsignc)


resid82A200v$c31 <- sapply(X = resid82A200v$V3,zsignc)
resid82A200v$c32 <- sapply(X = resid82A200v$V11,zsignc)
resid82A200v$c33 <- sapply(X = resid82A200v$V19,zsignc)
resid82A200v$c34 <- sapply(X = resid82A200v$V27,zsignc)
resid82A200v$c35 <- sapply(X = resid82A200v$V35,zsignc)
resid82A200v$c36 <- sapply(X = resid82A200v$V43,zsignc)
resid82A200v$c37 <- sapply(X = resid82A200v$V51,zsignc)
resid82A200v$c38 <- sapply(X = resid82A200v$V59,zsignc)


resid82A200v$c41 <- sapply(X = resid82A200v$V4,zsignc)
resid82A200v$c42 <- sapply(X = resid82A200v$V12,zsignc)
resid82A200v$c43 <- sapply(X = resid82A200v$V20,zsignc)
resid82A200v$c44 <- sapply(X = resid82A200v$V28,zsignc)
resid82A200v$c45 <- sapply(X = resid82A200v$V36,zsignc)
resid82A200v$c46 <- sapply(X = resid82A200v$V44,zsignc)
resid82A200v$c47 <- sapply(X = resid82A200v$V52,zsignc)
resid82A200v$c48 <- sapply(X = resid82A200v$V60,zsignc)


resid82A200v$c51 <- sapply(X = resid82A200v$V5,zsignc)
resid82A200v$c52 <- sapply(X = resid82A200v$V13,zsignc)
resid82A200v$c53 <- sapply(X = resid82A200v$V21,zsignc)
resid82A200v$c54 <- sapply(X = resid82A200v$V29,zsignc)
resid82A200v$c55 <- sapply(X = resid82A200v$V37,zsignc)
resid82A200v$c56 <- sapply(X = resid82A200v$V45,zsignc)
resid82A200v$c57 <- sapply(X = resid82A200v$V53,zsignc)
resid82A200v$c58 <- sapply(X = resid82A200v$V61,zsignc)


resid82A200v$c61 <- sapply(X = resid82A200v$V6,zsignc)
resid82A200v$c62 <- sapply(X = resid82A200v$V14,zsignc)
resid82A200v$c63 <- sapply(X = resid82A200v$V22,zsignc)
resid82A200v$c64 <- sapply(X = resid82A200v$V30,zsignc)
resid82A200v$c65 <- sapply(X = resid82A200v$V38,zsignc)
resid82A200v$c66 <- sapply(X = resid82A200v$V46,zsignc)
resid82A200v$c67 <- sapply(X = resid82A200v$V54,zsignc)
resid82A200v$c68 <- sapply(X = resid82A200v$V62,zsignc)


resid82A200v$c71 <- sapply(X = resid82A200v$V7,zsignc)
resid82A200v$c72 <- sapply(X = resid82A200v$V15,zsignc)
resid82A200v$c73 <- sapply(X = resid82A200v$V23,zsignc)
resid82A200v$c74 <- sapply(X = resid82A200v$V31,zsignc)
resid82A200v$c75 <- sapply(X = resid82A200v$V39,zsignc)
resid82A200v$c76 <- sapply(X = resid82A200v$V47,zsignc)
resid82A200v$c77 <- sapply(X = resid82A200v$V55,zsignc)
resid82A200v$c78 <- sapply(X = resid82A200v$V63,zsignc)


resid82A200v$c81 <- sapply(X = resid82A200v$V8,zsignc)
resid82A200v$c82 <- sapply(X = resid82A200v$V16,zsignc)
resid82A200v$c83 <- sapply(X = resid82A200v$V24,zsignc)
resid82A200v$c84 <- sapply(X = resid82A200v$V32,zsignc)
resid82A200v$c85 <- sapply(X = resid82A200v$V40,zsignc)
resid82A200v$c86 <- sapply(X = resid82A200v$V48,zsignc)
resid82A200v$c87 <- sapply(X = resid82A200v$V56,zsignc)
resid82A200v$c88 <- sapply(X = resid82A200v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82A200<- resid82A200v[,66:129]

datosLong82A200<-pivot_longer(zsig82A200, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82A200$N <- gl(1,320000, labels = 100)
datosLong82A200$Categoria<- rep(nombres8,5000)

tabla82A200<-table(datosLong82A200$Categoria,datosLong82A200$Medida)
tabla82A200p<- as.data.frame(prop.table(tabla82A200,1))
colnames(tabla82A200p)<- c("Patrón", "Dependencia", "P")
tabla82A200p$Cat<- rep(8,192)
tabla82A200p$Matriz<- rep("A",192)
tabla82A200p$N<- rep(200,192)


## N = 100, matriz 2B

transiz82B100<-alply(seq82B100,1,transi)
resid82B100<- llply(transiz82B100, function(x) resid2(x))

dim8<- ldply(resid82B100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

#resid82B100 <- resid82B100[-falsos]
resid82B100v<- ldply(resid82B100,as.vector)
buenos<- complete.cases(resid82B100v)
resid82B100v<- resid82B100v[buenos,]
resid82B100v<- resid82B100v[1:5000,]


resid82B100v$c11 <- sapply(X = resid82B100v$V1,zsignc)
resid82B100v$c12 <- sapply(X = resid82B100v$V9,zsignc)
resid82B100v$c13 <- sapply(X = resid82B100v$V17,zsignc)
resid82B100v$c14 <- sapply(X = resid82B100v$V25,zsignc)
resid82B100v$c15 <- sapply(X = resid82B100v$V33,zsignc)
resid82B100v$c16 <- sapply(X = resid82B100v$V41,zsignc)
resid82B100v$c17 <- sapply(X = resid82B100v$V49,zsignc)
resid82B100v$c18 <- sapply(X = resid82B100v$V57,zsignc)


resid82B100v$c21 <- sapply(X = resid82B100v$V2,zsignc)
resid82B100v$c22 <- sapply(X = resid82B100v$V10,zsignc)
resid82B100v$c23 <- sapply(X = resid82B100v$V18,zsignc)
resid82B100v$c24 <- sapply(X = resid82B100v$V26,zsignc)
resid82B100v$c25 <- sapply(X = resid82B100v$V34,zsignc)
resid82B100v$c26 <- sapply(X = resid82B100v$V42,zsignc)
resid82B100v$c27 <- sapply(X = resid82B100v$V50,zsignc)
resid82B100v$c28 <- sapply(X = resid82B100v$V58,zsignc)


resid82B100v$c31 <- sapply(X = resid82B100v$V3,zsignc)
resid82B100v$c32 <- sapply(X = resid82B100v$V11,zsignc)
resid82B100v$c33 <- sapply(X = resid82B100v$V19,zsignc)
resid82B100v$c34 <- sapply(X = resid82B100v$V27,zsignc)
resid82B100v$c35 <- sapply(X = resid82B100v$V35,zsignc)
resid82B100v$c36 <- sapply(X = resid82B100v$V43,zsignc)
resid82B100v$c37 <- sapply(X = resid82B100v$V51,zsignc)
resid82B100v$c38 <- sapply(X = resid82B100v$V59,zsignc)


resid82B100v$c41 <- sapply(X = resid82B100v$V4,zsignc)
resid82B100v$c42 <- sapply(X = resid82B100v$V12,zsignc)
resid82B100v$c43 <- sapply(X = resid82B100v$V20,zsignc)
resid82B100v$c44 <- sapply(X = resid82B100v$V28,zsignc)
resid82B100v$c45 <- sapply(X = resid82B100v$V36,zsignc)
resid82B100v$c46 <- sapply(X = resid82B100v$V44,zsignc)
resid82B100v$c47 <- sapply(X = resid82B100v$V52,zsignc)
resid82B100v$c48 <- sapply(X = resid82B100v$V60,zsignc)


resid82B100v$c51 <- sapply(X = resid82B100v$V5,zsignc)
resid82B100v$c52 <- sapply(X = resid82B100v$V13,zsignc)
resid82B100v$c53 <- sapply(X = resid82B100v$V21,zsignc)
resid82B100v$c54 <- sapply(X = resid82B100v$V29,zsignc)
resid82B100v$c55 <- sapply(X = resid82B100v$V37,zsignc)
resid82B100v$c56 <- sapply(X = resid82B100v$V45,zsignc)
resid82B100v$c57 <- sapply(X = resid82B100v$V53,zsignc)
resid82B100v$c58 <- sapply(X = resid82B100v$V61,zsignc)


resid82B100v$c61 <- sapply(X = resid82B100v$V6,zsignc)
resid82B100v$c62 <- sapply(X = resid82B100v$V14,zsignc)
resid82B100v$c63 <- sapply(X = resid82B100v$V22,zsignc)
resid82B100v$c64 <- sapply(X = resid82B100v$V30,zsignc)
resid82B100v$c65 <- sapply(X = resid82B100v$V38,zsignc)
resid82B100v$c66 <- sapply(X = resid82B100v$V46,zsignc)
resid82B100v$c67 <- sapply(X = resid82B100v$V54,zsignc)
resid82B100v$c68 <- sapply(X = resid82B100v$V62,zsignc)


resid82B100v$c71 <- sapply(X = resid82B100v$V7,zsignc)
resid82B100v$c72 <- sapply(X = resid82B100v$V15,zsignc)
resid82B100v$c73 <- sapply(X = resid82B100v$V23,zsignc)
resid82B100v$c74 <- sapply(X = resid82B100v$V31,zsignc)
resid82B100v$c75 <- sapply(X = resid82B100v$V39,zsignc)
resid82B100v$c76 <- sapply(X = resid82B100v$V47,zsignc)
resid82B100v$c77 <- sapply(X = resid82B100v$V55,zsignc)
resid82B100v$c78 <- sapply(X = resid82B100v$V63,zsignc)


resid82B100v$c81 <- sapply(X = resid82B100v$V8,zsignc)
resid82B100v$c82 <- sapply(X = resid82B100v$V16,zsignc)
resid82B100v$c83 <- sapply(X = resid82B100v$V24,zsignc)
resid82B100v$c84 <- sapply(X = resid82B100v$V32,zsignc)
resid82B100v$c85 <- sapply(X = resid82B100v$V40,zsignc)
resid82B100v$c86 <- sapply(X = resid82B100v$V48,zsignc)
resid82B100v$c87 <- sapply(X = resid82B100v$V56,zsignc)
resid82B100v$c88 <- sapply(X = resid82B100v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82B100<- resid82B100v[,66:129]

datosLong82B100<-pivot_longer(zsig82B100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82B100$N <- gl(1,320000, labels = 100)
datosLong82B100$Categoria<- rep(nombres8,5000)

tabla82B100<-table(datosLong82B100$Categoria,datosLong82B100$Medida)
tabla82B100p<- as.data.frame(prop.table(tabla82B100,1))
colnames(tabla82B100p)<- c("Patrón", "Dependencia", "P")
tabla82B100p$Cat<- rep(8,192)
tabla82B100p$Matriz<- rep("B",192)
tabla82B100p$N<- rep(100,192)

####

## N = 150, matriz 2B

transiz82B150<-alply(seq82B150,1,transi)
resid82B150<- llply(transiz82B150, function(x) resid2(x))

dim8<- ldply(resid82B150,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

#resid82B150<-resid82B150[-falsos]
resid82B150v<- ldply(resid82B150,as.vector)
buenos<- complete.cases(resid82B150v)
resid82B150v<- resid82B150v[buenos,]
resid82B150v<- resid82B150v[1:5000,]


resid82B150v$c11 <- sapply(X = resid82B150v$V1,zsignc)
resid82B150v$c12 <- sapply(X = resid82B150v$V9,zsignc)
resid82B150v$c13 <- sapply(X = resid82B150v$V17,zsignc)
resid82B150v$c14 <- sapply(X = resid82B150v$V25,zsignc)
resid82B150v$c15 <- sapply(X = resid82B150v$V33,zsignc)
resid82B150v$c16 <- sapply(X = resid82B150v$V41,zsignc)
resid82B150v$c17 <- sapply(X = resid82B150v$V49,zsignc)
resid82B150v$c18 <- sapply(X = resid82B150v$V57,zsignc)


resid82B150v$c21 <- sapply(X = resid82B150v$V2,zsignc)
resid82B150v$c22 <- sapply(X = resid82B150v$V10,zsignc)
resid82B150v$c23 <- sapply(X = resid82B150v$V18,zsignc)
resid82B150v$c24 <- sapply(X = resid82B150v$V26,zsignc)
resid82B150v$c25 <- sapply(X = resid82B150v$V34,zsignc)
resid82B150v$c26 <- sapply(X = resid82B150v$V42,zsignc)
resid82B150v$c27 <- sapply(X = resid82B150v$V50,zsignc)
resid82B150v$c28 <- sapply(X = resid82B150v$V58,zsignc)


resid82B150v$c31 <- sapply(X = resid82B150v$V3,zsignc)
resid82B150v$c32 <- sapply(X = resid82B150v$V11,zsignc)
resid82B150v$c33 <- sapply(X = resid82B150v$V19,zsignc)
resid82B150v$c34 <- sapply(X = resid82B150v$V27,zsignc)
resid82B150v$c35 <- sapply(X = resid82B150v$V35,zsignc)
resid82B150v$c36 <- sapply(X = resid82B150v$V43,zsignc)
resid82B150v$c37 <- sapply(X = resid82B150v$V51,zsignc)
resid82B150v$c38 <- sapply(X = resid82B150v$V59,zsignc)

resid82B150v$c41 <- sapply(X = resid82B150v$V4,zsignc)
resid82B150v$c42 <- sapply(X = resid82B150v$V12,zsignc)
resid82B150v$c43 <- sapply(X = resid82B150v$V20,zsignc)
resid82B150v$c44 <- sapply(X = resid82B150v$V28,zsignc)
resid82B150v$c45 <- sapply(X = resid82B150v$V36,zsignc)
resid82B150v$c46 <- sapply(X = resid82B150v$V44,zsignc)
resid82B150v$c47 <- sapply(X = resid82B150v$V52,zsignc)
resid82B150v$c48 <- sapply(X = resid82B150v$V60,zsignc)


resid82B150v$c51 <- sapply(X = resid82B150v$V5,zsignc)
resid82B150v$c52 <- sapply(X = resid82B150v$V13,zsignc)
resid82B150v$c53 <- sapply(X = resid82B150v$V21,zsignc)
resid82B150v$c54 <- sapply(X = resid82B150v$V29,zsignc)
resid82B150v$c55 <- sapply(X = resid82B150v$V37,zsignc)
resid82B150v$c56 <- sapply(X = resid82B150v$V45,zsignc)
resid82B150v$c57 <- sapply(X = resid82B150v$V53,zsignc)
resid82B150v$c58 <- sapply(X = resid82B150v$V61,zsignc)


resid82B150v$c61 <- sapply(X = resid82B150v$V6,zsignc)
resid82B150v$c62 <- sapply(X = resid82B150v$V14,zsignc)
resid82B150v$c63 <- sapply(X = resid82B150v$V22,zsignc)
resid82B150v$c64 <- sapply(X = resid82B150v$V30,zsignc)
resid82B150v$c65 <- sapply(X = resid82B150v$V38,zsignc)
resid82B150v$c66 <- sapply(X = resid82B150v$V46,zsignc)
resid82B150v$c67 <- sapply(X = resid82B150v$V54,zsignc)
resid82B150v$c68 <- sapply(X = resid82B150v$V62,zsignc)


resid82B150v$c71 <- sapply(X = resid82B150v$V7,zsignc)
resid82B150v$c72 <- sapply(X = resid82B150v$V15,zsignc)
resid82B150v$c73 <- sapply(X = resid82B150v$V23,zsignc)
resid82B150v$c74 <- sapply(X = resid82B150v$V31,zsignc)
resid82B150v$c75 <- sapply(X = resid82B150v$V39,zsignc)
resid82B150v$c76 <- sapply(X = resid82B150v$V47,zsignc)
resid82B150v$c77 <- sapply(X = resid82B150v$V55,zsignc)
resid82B150v$c78 <- sapply(X = resid82B150v$V63,zsignc)


resid82B150v$c81 <- sapply(X = resid82B150v$V8,zsignc)
resid82B150v$c82 <- sapply(X = resid82B150v$V16,zsignc)
resid82B150v$c83 <- sapply(X = resid82B150v$V24,zsignc)
resid82B150v$c84 <- sapply(X = resid82B150v$V32,zsignc)
resid82B150v$c85 <- sapply(X = resid82B150v$V40,zsignc)
resid82B150v$c86 <- sapply(X = resid82B150v$V48,zsignc)
resid82B150v$c87 <- sapply(X = resid82B150v$V56,zsignc)
resid82B150v$c88 <- sapply(X = resid82B150v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82B150<- resid82B150v[,66:129]

datosLong82B150<-pivot_longer(zsig82B150, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82B150$N <- gl(1,320000, labels = 150)

datosLong82B150$Categoria<- rep(nombres8,5000)

tabla82B150<-table(datosLong82B150$Categoria,datosLong82B150$Medida)
tabla82B150p<- as.data.frame(prop.table(tabla82B150,1))
colnames(tabla82B150p)<- c("Patrón", "Dependencia", "P")
tabla82B150p$Cat<- rep(8,192)
tabla82B150p$Matriz<- rep("B",192)
tabla82B150p$N<- rep(150,192)

####


## N = 200, matriz 2B

transiz82B200<-alply(seq82B200,1,transi)
resid82B200<- llply(transiz82B200, function(x) resid2(x))

dim8<- ldply(resid82B200,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)


resid82B200v<- ldply(resid82B200,as.vector)
buenos<- complete.cases(resid82B200v)
resid82B200v<- resid82B200v[buenos,]
resid82B200v<- resid82B200v[1:5000,]


resid82B200v$c11 <- sapply(X = resid82B200v$V1,zsignc)
resid82B200v$c12 <- sapply(X = resid82B200v$V9,zsignc)
resid82B200v$c13 <- sapply(X = resid82B200v$V17,zsignc)
resid82B200v$c14 <- sapply(X = resid82B200v$V25,zsignc)
resid82B200v$c15 <- sapply(X = resid82B200v$V33,zsignc)
resid82B200v$c16 <- sapply(X = resid82B200v$V41,zsignc)
resid82B200v$c17 <- sapply(X = resid82B200v$V49,zsignc)
resid82B200v$c18 <- sapply(X = resid82B200v$V57,zsignc)


resid82B200v$c21 <- sapply(X = resid82B200v$V2,zsignc)
resid82B200v$c22 <- sapply(X = resid82B200v$V10,zsignc)
resid82B200v$c23 <- sapply(X = resid82B200v$V18,zsignc)
resid82B200v$c24 <- sapply(X = resid82B200v$V26,zsignc)
resid82B200v$c25 <- sapply(X = resid82B200v$V34,zsignc)
resid82B200v$c26 <- sapply(X = resid82B200v$V42,zsignc)
resid82B200v$c27 <- sapply(X = resid82B200v$V50,zsignc)
resid82B200v$c28 <- sapply(X = resid82B200v$V58,zsignc)


resid82B200v$c31 <- sapply(X = resid82B200v$V3,zsignc)
resid82B200v$c32 <- sapply(X = resid82B200v$V11,zsignc)
resid82B200v$c33 <- sapply(X = resid82B200v$V19,zsignc)
resid82B200v$c34 <- sapply(X = resid82B200v$V27,zsignc)
resid82B200v$c35 <- sapply(X = resid82B200v$V35,zsignc)
resid82B200v$c36 <- sapply(X = resid82B200v$V43,zsignc)
resid82B200v$c37 <- sapply(X = resid82B200v$V51,zsignc)
resid82B200v$c38 <- sapply(X = resid82B200v$V59,zsignc)


resid82B200v$c41 <- sapply(X = resid82B200v$V4,zsignc)
resid82B200v$c42 <- sapply(X = resid82B200v$V12,zsignc)
resid82B200v$c43 <- sapply(X = resid82B200v$V20,zsignc)
resid82B200v$c44 <- sapply(X = resid82B200v$V28,zsignc)
resid82B200v$c45 <- sapply(X = resid82B200v$V36,zsignc)
resid82B200v$c46 <- sapply(X = resid82B200v$V44,zsignc)
resid82B200v$c47 <- sapply(X = resid82B200v$V52,zsignc)
resid82B200v$c48 <- sapply(X = resid82B200v$V60,zsignc)


resid82B200v$c51 <- sapply(X = resid82B200v$V5,zsignc)
resid82B200v$c52 <- sapply(X = resid82B200v$V13,zsignc)
resid82B200v$c53 <- sapply(X = resid82B200v$V21,zsignc)
resid82B200v$c54 <- sapply(X = resid82B200v$V29,zsignc)
resid82B200v$c55 <- sapply(X = resid82B200v$V37,zsignc)
resid82B200v$c56 <- sapply(X = resid82B200v$V45,zsignc)
resid82B200v$c57 <- sapply(X = resid82B200v$V53,zsignc)
resid82B200v$c58 <- sapply(X = resid82B200v$V61,zsignc)


resid82B200v$c61 <- sapply(X = resid82B200v$V6,zsignc)
resid82B200v$c62 <- sapply(X = resid82B200v$V14,zsignc)
resid82B200v$c63 <- sapply(X = resid82B200v$V22,zsignc)
resid82B200v$c64 <- sapply(X = resid82B200v$V30,zsignc)
resid82B200v$c65 <- sapply(X = resid82B200v$V38,zsignc)
resid82B200v$c66 <- sapply(X = resid82B200v$V46,zsignc)
resid82B200v$c67 <- sapply(X = resid82B200v$V54,zsignc)
resid82B200v$c68 <- sapply(X = resid82B200v$V62,zsignc)


resid82B200v$c71 <- sapply(X = resid82B200v$V7,zsignc)
resid82B200v$c72 <- sapply(X = resid82B200v$V15,zsignc)
resid82B200v$c73 <- sapply(X = resid82B200v$V23,zsignc)
resid82B200v$c74 <- sapply(X = resid82B200v$V31,zsignc)
resid82B200v$c75 <- sapply(X = resid82B200v$V39,zsignc)
resid82B200v$c76 <- sapply(X = resid82B200v$V47,zsignc)
resid82B200v$c77 <- sapply(X = resid82B200v$V55,zsignc)
resid82B200v$c78 <- sapply(X = resid82B200v$V63,zsignc)


resid82B200v$c81 <- sapply(X = resid82B200v$V8,zsignc)
resid82B200v$c82 <- sapply(X = resid82B200v$V16,zsignc)
resid82B200v$c83 <- sapply(X = resid82B200v$V24,zsignc)
resid82B200v$c84 <- sapply(X = resid82B200v$V32,zsignc)
resid82B200v$c85 <- sapply(X = resid82B200v$V40,zsignc)
resid82B200v$c86 <- sapply(X = resid82B200v$V48,zsignc)
resid82B200v$c87 <- sapply(X = resid82B200v$V56,zsignc)
resid82B200v$c88 <- sapply(X = resid82B200v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82B200<- resid82B200v[,66:129]

datosLong82B200<-pivot_longer(zsig82B200, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong82B200$N <- gl(1,320000, labels = 200)

datosLong82B200$Categoria<- rep(nombres8,5000)

tabla82B200<-table(datosLong82B200$Categoria,datosLong82B200$Medida)
tabla82B200p<- as.data.frame(prop.table(tabla82B200,1))
colnames(tabla82B200p)<- c("Patrón", "Dependencia", "P")
tabla82B200p$Cat<- rep(8,192)
tabla82B200p$Matriz<- rep("B",192)
tabla82B200p$N<- rep(200,192)


####



## N = 100, matriz 2C

transiz82C100<-alply(seq82C100,1,transi)
resid82C100<- llply(transiz82C100, function(x) resid2(x))

dim8<- ldply(resid82C100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid82C100 <- resid82C100[-falsos]
resid82C100v<- ldply(resid82C100,as.vector)
buenos<- complete.cases(resid82C100v)
resid82C100v<- resid82C100v[buenos,]
resid82C100v<- resid82C100v[1:5000,]


resid82C100v$c11 <- sapply(X = resid82C100v$V1,zsignc)
resid82C100v$c12 <- sapply(X = resid82C100v$V9,zsignc)
resid82C100v$c13 <- sapply(X = resid82C100v$V17,zsignc)
resid82C100v$c14 <- sapply(X = resid82C100v$V25,zsignc)
resid82C100v$c15 <- sapply(X = resid82C100v$V33,zsignc)
resid82C100v$c16 <- sapply(X = resid82C100v$V41,zsignc)
resid82C100v$c17 <- sapply(X = resid82C100v$V49,zsignc)
resid82C100v$c18 <- sapply(X = resid82C100v$V57,zsignc)


resid82C100v$c21 <- sapply(X = resid82C100v$V2,zsignc)
resid82C100v$c22 <- sapply(X = resid82C100v$V10,zsignc)
resid82C100v$c23 <- sapply(X = resid82C100v$V18,zsignc)
resid82C100v$c24 <- sapply(X = resid82C100v$V26,zsignc)
resid82C100v$c25 <- sapply(X = resid82C100v$V34,zsignc)
resid82C100v$c26 <- sapply(X = resid82C100v$V42,zsignc)
resid82C100v$c27 <- sapply(X = resid82C100v$V50,zsignc)
resid82C100v$c28 <- sapply(X = resid82C100v$V58,zsignc)


resid82C100v$c81 <- sapply(X = resid82C100v$V8,zsignc)
resid82C100v$c82 <- sapply(X = resid82C100v$V16,zsignc)
resid82C100v$c83 <- sapply(X = resid82C100v$V24,zsignc)
resid82C100v$c84 <- sapply(X = resid82C100v$V32,zsignc)
resid82C100v$c85 <- sapply(X = resid82C100v$V40,zsignc)
resid82C100v$c86 <- sapply(X = resid82C100v$V48,zsignc)
resid82C100v$c87 <- sapply(X = resid82C100v$V56,zsignc)
resid82C100v$c88 <- sapply(X = resid82C100v$V64,zsignc)

resid82C100v$c31 <- sapply(X = resid82C100v$V3,zsignc)
resid82C100v$c32 <- sapply(X = resid82C100v$V11,zsignc)
resid82C100v$c33 <- sapply(X = resid82C100v$V19,zsignc)
resid82C100v$c34 <- sapply(X = resid82C100v$V27,zsignc)
resid82C100v$c35 <- sapply(X = resid82C100v$V35,zsignc)
resid82C100v$c36 <- sapply(X = resid82C100v$V43,zsignc)
resid82C100v$c37 <- sapply(X = resid82C100v$V51,zsignc)
resid82C100v$c38 <- sapply(X = resid82C100v$V59,zsignc)


resid82C100v$c41 <- sapply(X = resid82C100v$V4,zsignc)
resid82C100v$c42 <- sapply(X = resid82C100v$V12,zsignc)
resid82C100v$c43 <- sapply(X = resid82C100v$V20,zsignc)
resid82C100v$c44 <- sapply(X = resid82C100v$V28,zsignc)
resid82C100v$c45 <- sapply(X = resid82C100v$V36,zsignc)
resid82C100v$c46 <- sapply(X = resid82C100v$V44,zsignc)
resid82C100v$c47 <- sapply(X = resid82C100v$V52,zsignc)
resid82C100v$c48 <- sapply(X = resid82C100v$V60,zsignc)


resid82C100v$c51 <- sapply(X = resid82C100v$V5,zsignc)
resid82C100v$c52 <- sapply(X = resid82C100v$V13,zsignc)
resid82C100v$c53 <- sapply(X = resid82C100v$V21,zsignc)
resid82C100v$c54 <- sapply(X = resid82C100v$V29,zsignc)
resid82C100v$c55 <- sapply(X = resid82C100v$V37,zsignc)
resid82C100v$c56 <- sapply(X = resid82C100v$V45,zsignc)
resid82C100v$c57 <- sapply(X = resid82C100v$V53,zsignc)
resid82C100v$c58 <- sapply(X = resid82C100v$V61,zsignc)


resid82C100v$c61 <- sapply(X = resid82C100v$V6,zsignc)
resid82C100v$c62 <- sapply(X = resid82C100v$V14,zsignc)
resid82C100v$c63 <- sapply(X = resid82C100v$V22,zsignc)
resid82C100v$c64 <- sapply(X = resid82C100v$V30,zsignc)
resid82C100v$c65 <- sapply(X = resid82C100v$V38,zsignc)
resid82C100v$c66 <- sapply(X = resid82C100v$V46,zsignc)
resid82C100v$c67 <- sapply(X = resid82C100v$V54,zsignc)
resid82C100v$c68 <- sapply(X = resid82C100v$V62,zsignc)


resid82C100v$c71 <- sapply(X = resid82C100v$V7,zsignc)
resid82C100v$c72 <- sapply(X = resid82C100v$V15,zsignc)
resid82C100v$c73 <- sapply(X = resid82C100v$V23,zsignc)
resid82C100v$c74 <- sapply(X = resid82C100v$V31,zsignc)
resid82C100v$c75 <- sapply(X = resid82C100v$V39,zsignc)
resid82C100v$c76 <- sapply(X = resid82C100v$V47,zsignc)
resid82C100v$c77 <- sapply(X = resid82C100v$V55,zsignc)
resid82C100v$c78 <- sapply(X = resid82C100v$V63,zsignc)


resid82C100v$c81 <- sapply(X = resid82C100v$V8,zsignc)
resid82C100v$c82 <- sapply(X = resid82C100v$V16,zsignc)
resid82C100v$c83 <- sapply(X = resid82C100v$V24,zsignc)
resid82C100v$c84 <- sapply(X = resid82C100v$V32,zsignc)
resid82C100v$c85 <- sapply(X = resid82C100v$V40,zsignc)
resid82C100v$c86 <- sapply(X = resid82C100v$V48,zsignc)
resid82C100v$c87 <- sapply(X = resid82C100v$V56,zsignc)
resid82C100v$c88 <- sapply(X = resid82C100v$V64,zsignc)




# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82C100<- resid82C100v[,66:129]

datosLong82C100<-pivot_longer(zsig82C100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82C100$N <- gl(1,320000, labels = 100)

datosLong82C100$Categoria<- rep(nombres8,5000)

tabla82C100<-table(datosLong82C100$Categoria,datosLong82C100$Medida)
tabla82C100p<- as.data.frame(prop.table(tabla82C100,1))
colnames(tabla82C100p)<- c("Patrón", "Dependencia", "P")
tabla82C100p$Cat<- rep(8,192)
tabla82C100p$Matriz<- rep("C",192)
tabla82C100p$N<- rep(100,192)


####

## N = 150, matriz 2C

transiz82C150<-alply(seq82C150,1,transi)
resid82C150<- llply(transiz82C150, function(x) resid2(x))

dim8<- ldply(resid82C150,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid82C150 <- resid82C150[-falsos]
resid82C150v<- ldply(resid82C150,as.vector)
buenos<- complete.cases(resid82C150v)
resid82C150v<- resid82C150v[buenos,]
resid82C150v<- resid82C150v[1:5000,]


resid82C150v$c11 <- sapply(X = resid82C150v$V1,zsignc)
resid82C150v$c12 <- sapply(X = resid82C150v$V9,zsignc)
resid82C150v$c13 <- sapply(X = resid82C150v$V17,zsignc)
resid82C150v$c14 <- sapply(X = resid82C150v$V25,zsignc)
resid82C150v$c15 <- sapply(X = resid82C150v$V33,zsignc)
resid82C150v$c16 <- sapply(X = resid82C150v$V41,zsignc)
resid82C150v$c17 <- sapply(X = resid82C150v$V49,zsignc)
resid82C150v$c18 <- sapply(X = resid82C150v$V57,zsignc)


resid82C150v$c21 <- sapply(X = resid82C150v$V2,zsignc)
resid82C150v$c22 <- sapply(X = resid82C150v$V10,zsignc)
resid82C150v$c23 <- sapply(X = resid82C150v$V18,zsignc)
resid82C150v$c24 <- sapply(X = resid82C150v$V26,zsignc)
resid82C150v$c25 <- sapply(X = resid82C150v$V34,zsignc)
resid82C150v$c26 <- sapply(X = resid82C150v$V42,zsignc)
resid82C150v$c27 <- sapply(X = resid82C150v$V50,zsignc)
resid82C150v$c28 <- sapply(X = resid82C150v$V58,zsignc)


resid82C150v$c31 <- sapply(X = resid82C150v$V3,zsignc)
resid82C150v$c32 <- sapply(X = resid82C150v$V11,zsignc)
resid82C150v$c33 <- sapply(X = resid82C150v$V19,zsignc)
resid82C150v$c34 <- sapply(X = resid82C150v$V27,zsignc)
resid82C150v$c35 <- sapply(X = resid82C150v$V35,zsignc)
resid82C150v$c36 <- sapply(X = resid82C150v$V43,zsignc)
resid82C150v$c37 <- sapply(X = resid82C150v$V51,zsignc)
resid82C150v$c38 <- sapply(X = resid82C150v$V59,zsignc)


resid82C150v$c41 <- sapply(X = resid82C150v$V4,zsignc)
resid82C150v$c42 <- sapply(X = resid82C150v$V12,zsignc)
resid82C150v$c43 <- sapply(X = resid82C150v$V20,zsignc)
resid82C150v$c44 <- sapply(X = resid82C150v$V28,zsignc)
resid82C150v$c45 <- sapply(X = resid82C150v$V36,zsignc)
resid82C150v$c46 <- sapply(X = resid82C150v$V44,zsignc)
resid82C150v$c47 <- sapply(X = resid82C150v$V52,zsignc)
resid82C150v$c48 <- sapply(X = resid82C150v$V60,zsignc)


resid82C150v$c51 <- sapply(X = resid82C150v$V5,zsignc)
resid82C150v$c52 <- sapply(X = resid82C150v$V13,zsignc)
resid82C150v$c53 <- sapply(X = resid82C150v$V21,zsignc)
resid82C150v$c54 <- sapply(X = resid82C150v$V29,zsignc)
resid82C150v$c55 <- sapply(X = resid82C150v$V37,zsignc)
resid82C150v$c56 <- sapply(X = resid82C150v$V45,zsignc)
resid82C150v$c57 <- sapply(X = resid82C150v$V53,zsignc)
resid82C150v$c58 <- sapply(X = resid82C150v$V61,zsignc)


resid82C150v$c61 <- sapply(X = resid82C150v$V6,zsignc)
resid82C150v$c62 <- sapply(X = resid82C150v$V14,zsignc)
resid82C150v$c63 <- sapply(X = resid82C150v$V22,zsignc)
resid82C150v$c64 <- sapply(X = resid82C150v$V30,zsignc)
resid82C150v$c65 <- sapply(X = resid82C150v$V38,zsignc)
resid82C150v$c66 <- sapply(X = resid82C150v$V46,zsignc)
resid82C150v$c67 <- sapply(X = resid82C150v$V54,zsignc)
resid82C150v$c68 <- sapply(X = resid82C150v$V62,zsignc)


resid82C150v$c71 <- sapply(X = resid82C150v$V7,zsignc)
resid82C150v$c72 <- sapply(X = resid82C150v$V15,zsignc)
resid82C150v$c73 <- sapply(X = resid82C150v$V23,zsignc)
resid82C150v$c74 <- sapply(X = resid82C150v$V31,zsignc)
resid82C150v$c75 <- sapply(X = resid82C150v$V39,zsignc)
resid82C150v$c76 <- sapply(X = resid82C150v$V47,zsignc)
resid82C150v$c77 <- sapply(X = resid82C150v$V55,zsignc)
resid82C150v$c78 <- sapply(X = resid82C150v$V63,zsignc)


resid82C150v$c81 <- sapply(X = resid82C150v$V8,zsignc)
resid82C150v$c82 <- sapply(X = resid82C150v$V16,zsignc)
resid82C150v$c83 <- sapply(X = resid82C150v$V24,zsignc)
resid82C150v$c84 <- sapply(X = resid82C150v$V32,zsignc)
resid82C150v$c85 <- sapply(X = resid82C150v$V40,zsignc)
resid82C150v$c86 <- sapply(X = resid82C150v$V48,zsignc)
resid82C150v$c87 <- sapply(X = resid82C150v$V56,zsignc)
resid82C150v$c88 <- sapply(X = resid82C150v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82C150<- resid82C150v[,66:129]

datosLong82C150<-pivot_longer(zsig82C150, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong82C150$N <- gl(1,320000, labels = 150)

datosLong82C150$Categoria<- rep(nombres8,5000)

tabla82C150<-table(datosLong82C150$Categoria,datosLong82C150$Medida)
tabla82C150p<- as.data.frame(prop.table(tabla82C150,1))
colnames(tabla82C150p)<- c("Patrón", "Dependencia", "P")
tabla82C150p$Cat<- rep(8,192)
tabla82C150p$Matriz<- rep("C",192)
tabla82C150p$N<- rep(150,192)

####

## N = 200, matriz C

transiz82C200<-alply(seq82C200,1,transi)
resid82C200<- llply(transiz82C200, function(x) resid2(x))

dim8<- ldply(resid82C200,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

#resid82C200 <- resid82C200[-falsos]
resid82C200v<- ldply(resid82C200,as.vector)
buenos<- complete.cases(resid82C200v)
resid82C200v<- resid82C200v[buenos,]
resid82C200v<- resid82C200v[1:5000,]


resid82C200v$c11 <- sapply(X = resid82C200v$V1,zsignc)
resid82C200v$c12 <- sapply(X = resid82C200v$V9,zsignc)
resid82C200v$c13 <- sapply(X = resid82C200v$V17,zsignc)
resid82C200v$c14 <- sapply(X = resid82C200v$V25,zsignc)
resid82C200v$c15 <- sapply(X = resid82C200v$V33,zsignc)
resid82C200v$c16 <- sapply(X = resid82C200v$V41,zsignc)
resid82C200v$c17 <- sapply(X = resid82C200v$V49,zsignc)
resid82C200v$c18 <- sapply(X = resid82C200v$V57,zsignc)


resid82C200v$c21 <- sapply(X = resid82C200v$V2,zsignc)
resid82C200v$c22 <- sapply(X = resid82C200v$V10,zsignc)
resid82C200v$c23 <- sapply(X = resid82C200v$V18,zsignc)
resid82C200v$c24 <- sapply(X = resid82C200v$V26,zsignc)
resid82C200v$c25 <- sapply(X = resid82C200v$V34,zsignc)
resid82C200v$c26 <- sapply(X = resid82C200v$V42,zsignc)
resid82C200v$c27 <- sapply(X = resid82C200v$V50,zsignc)
resid82C200v$c28 <- sapply(X = resid82C200v$V58,zsignc)


resid82C200v$c31 <- sapply(X = resid82C200v$V3,zsignc)
resid82C200v$c32 <- sapply(X = resid82C200v$V11,zsignc)
resid82C200v$c33 <- sapply(X = resid82C200v$V19,zsignc)
resid82C200v$c34 <- sapply(X = resid82C200v$V27,zsignc)
resid82C200v$c35 <- sapply(X = resid82C200v$V35,zsignc)
resid82C200v$c36 <- sapply(X = resid82C200v$V43,zsignc)
resid82C200v$c37 <- sapply(X = resid82C200v$V51,zsignc)
resid82C200v$c38 <- sapply(X = resid82C200v$V59,zsignc)


resid82C200v$c41 <- sapply(X = resid82C200v$V4,zsignc)
resid82C200v$c42 <- sapply(X = resid82C200v$V12,zsignc)
resid82C200v$c43 <- sapply(X = resid82C200v$V20,zsignc)
resid82C200v$c44 <- sapply(X = resid82C200v$V28,zsignc)
resid82C200v$c45 <- sapply(X = resid82C200v$V36,zsignc)
resid82C200v$c46 <- sapply(X = resid82C200v$V44,zsignc)
resid82C200v$c47 <- sapply(X = resid82C200v$V52,zsignc)
resid82C200v$c48 <- sapply(X = resid82C200v$V60,zsignc)


resid82C200v$c51 <- sapply(X = resid82C200v$V5,zsignc)
resid82C200v$c52 <- sapply(X = resid82C200v$V13,zsignc)
resid82C200v$c53 <- sapply(X = resid82C200v$V21,zsignc)
resid82C200v$c54 <- sapply(X = resid82C200v$V29,zsignc)
resid82C200v$c55 <- sapply(X = resid82C200v$V37,zsignc)
resid82C200v$c56 <- sapply(X = resid82C200v$V45,zsignc)
resid82C200v$c57 <- sapply(X = resid82C200v$V53,zsignc)
resid82C200v$c58 <- sapply(X = resid82C200v$V61,zsignc)


resid82C200v$c61 <- sapply(X = resid82C200v$V6,zsignc)
resid82C200v$c62 <- sapply(X = resid82C200v$V14,zsignc)
resid82C200v$c63 <- sapply(X = resid82C200v$V22,zsignc)
resid82C200v$c64 <- sapply(X = resid82C200v$V30,zsignc)
resid82C200v$c65 <- sapply(X = resid82C200v$V38,zsignc)
resid82C200v$c66 <- sapply(X = resid82C200v$V46,zsignc)
resid82C200v$c67 <- sapply(X = resid82C200v$V54,zsignc)
resid82C200v$c68 <- sapply(X = resid82C200v$V62,zsignc)


resid82C200v$c71 <- sapply(X = resid82C200v$V7,zsignc)
resid82C200v$c72 <- sapply(X = resid82C200v$V15,zsignc)
resid82C200v$c73 <- sapply(X = resid82C200v$V23,zsignc)
resid82C200v$c74 <- sapply(X = resid82C200v$V31,zsignc)
resid82C200v$c75 <- sapply(X = resid82C200v$V39,zsignc)
resid82C200v$c76 <- sapply(X = resid82C200v$V47,zsignc)
resid82C200v$c77 <- sapply(X = resid82C200v$V55,zsignc)
resid82C200v$c78 <- sapply(X = resid82C200v$V63,zsignc)


resid82C200v$c81 <- sapply(X = resid82C200v$V8,zsignc)
resid82C200v$c82 <- sapply(X = resid82C200v$V16,zsignc)
resid82C200v$c83 <- sapply(X = resid82C200v$V24,zsignc)
resid82C200v$c84 <- sapply(X = resid82C200v$V32,zsignc)
resid82C200v$c85 <- sapply(X = resid82C200v$V40,zsignc)
resid82C200v$c86 <- sapply(X = resid82C200v$V48,zsignc)
resid82C200v$c87 <- sapply(X = resid82C200v$V56,zsignc)
resid82C200v$c88 <- sapply(X = resid82C200v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig82C200<- resid82C200v[,66:129]

datosLong82C200<-pivot_longer(zsig82C200, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong82C200$N <- gl(1,320000, labels = 100)

datosLong82C200$Categoria<- rep(nombres8,5000)

tabla82C200<-table(datosLong82C200$Categoria,datosLong82C200$Medida)
tabla82C200p<- as.data.frame(prop.table(tabla82C200,1))
colnames(tabla82C200p)<- c("Patrón", "Dependencia", "P")
tabla82C200p$Cat<- rep(8,192)
tabla82C200p$Matriz<- rep("C",192)
tabla82C200p$N<- rep(200,192)

# Tabla resumen de datos


tabla8C<- rbind(tabla82A100p,tabla82A150p,tabla82A200p,tabla82B100p,tabla82B150p,
               tabla82B200p, tabla82C100p,tabla82C150p,tabla82C200p)

patron_ab<- tabla8C[tabla8C$Patrón=="AB" & tabla8C$Dependencia== 'Act',]
patron_cd<- tabla8C[tabla8C$Patrón=="CD" & tabla8C$Dependencia== 'Act',]

patron_ef<- tabla8C[tabla8C$Patrón=="EF" & tabla8C$Dependencia== 'Act',]
patron_gh<- tabla8C[tabla8C$Patrón=="GH" & tabla8C$Dependencia== 'Act',]

patron_be<- tabla8C[tabla8C$Patrón=="BE" & tabla8C$Dependencia== 'Act',]
patron_dh<- tabla8C[tabla8C$Patrón=="DH" & tabla8C$Dependencia== 'Act',]




