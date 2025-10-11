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

set.seed(123456)


## MCTM library is removed from R. 
## Formerly available versions can be obtained from the archive.

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

# Matriz nula A

m8A <- matrix(c(1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
                8, 8, byrow = T)

# Matriz B con p = 6/8, 1 casilla significativa (2,3)


m8B <- matrix(c(1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 6/8,1/28,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
                8, 8, byrow = T)


# Matriz C con p = 6/8, 2 casillas significativa (2,3),(4,3)


m8C <- matrix(c(1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 6/8,1/28,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/28, 1/28, 6/8,1/28,1/28, 1/28, 1/28,1/28,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8,
                1/8, 1/8, 1/8,1/8,1/8, 1/8, 1/8,1/8), 
                8, 8, byrow = T)


# Condición con 8 estados


seq8A20 <- genMarkov(n = 15000, transMat = m8A, 
                     chainLen = 20, wide = TRUE)
seq8B20 <- genMarkov(n = 15000, transMat = m8B, 
                     chainLen = 20, wide = TRUE)
seq8C20 <- genMarkov(n = 15000, transMat = m8C, 
                     chainLen = 20, wide = TRUE)


seq8A50 <- genMarkov(n = 15000, transMat = m8A, 
                     chainLen = 50, wide = TRUE)
seq8B50 <- genMarkov(n = 15000, transMat = m8B, 
                     chainLen = 50, wide = TRUE)
seq8C50 <- genMarkov(n = 15000, transMat = m8C, 
                     chainLen = 50, wide = TRUE)


seq8A100 <- genMarkov(n = 15000, transMat = m8A, 
                      chainLen = 100, wide = TRUE)
seq8B100 <- genMarkov(n = 15000, transMat = m8B, 
                      chainLen = 100, wide = TRUE)
seq8C100 <- genMarkov(n = 15000, transMat = m8C, 
                      chainLen = 100, wide = TRUE)


seq8A20<- seq8A20[,-1]
seq8B20<- seq8B20[,-1]
seq8C20<- seq8C20[,-1]

seq8A50<- seq8A50[,-1]
seq8B50<- seq8B50[,-1]
seq8C50<- seq8C50[,-1]

seq8A100<- seq8A100[,-1]
seq8B100<- seq8B100[,-1]
seq8C100<- seq8C100[,-1]

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

## N = 20, matriz N

transiz8A20<-alply(seq8A20,1,transi)
resid8A20<- llply(transiz8A20, function(x) resid2(x))

dim8<- ldply(resid8A20,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8A20 <- resid8A20[-falsos]
resid8A20v<- ldply(resid8A20,as.vector)
buenos<- complete.cases(resid8A20v)
resid8A20v<- resid8A20v[buenos,]
resid8A20v<- resid8A20v[1:5000,]


resid8A20v$c11 <- sapply(X = resid8A20v$V1,zsignc)
resid8A20v$c12 <- sapply(X = resid8A20v$V9,zsignc)
resid8A20v$c13 <- sapply(X = resid8A20v$V17,zsignc)
resid8A20v$c14 <- sapply(X = resid8A20v$V25,zsignc)
resid8A20v$c15 <- sapply(X = resid8A20v$V33,zsignc)
resid8A20v$c16 <- sapply(X = resid8A20v$V41,zsignc)
resid8A20v$c17 <- sapply(X = resid8A20v$V49,zsignc)
resid8A20v$c18 <- sapply(X = resid8A20v$V57,zsignc)


resid8A20v$c21 <- sapply(X = resid8A20v$V2,zsignc)
resid8A20v$c22 <- sapply(X = resid8A20v$V10,zsignc)
resid8A20v$c23 <- sapply(X = resid8A20v$V18,zsignc)
resid8A20v$c24 <- sapply(X = resid8A20v$V26,zsignc)
resid8A20v$c25 <- sapply(X = resid8A20v$V34,zsignc)
resid8A20v$c26 <- sapply(X = resid8A20v$V42,zsignc)
resid8A20v$c27 <- sapply(X = resid8A20v$V50,zsignc)
resid8A20v$c28 <- sapply(X = resid8A20v$V58,zsignc)


resid8A20v$c31 <- sapply(X = resid8A20v$V3,zsignc)
resid8A20v$c32 <- sapply(X = resid8A20v$V11,zsignc)
resid8A20v$c33 <- sapply(X = resid8A20v$V19,zsignc)
resid8A20v$c34 <- sapply(X = resid8A20v$V27,zsignc)
resid8A20v$c35 <- sapply(X = resid8A20v$V35,zsignc)
resid8A20v$c36 <- sapply(X = resid8A20v$V43,zsignc)
resid8A20v$c37 <- sapply(X = resid8A20v$V51,zsignc)
resid8A20v$c38 <- sapply(X = resid8A20v$V59,zsignc)


resid8A20v$c41 <- sapply(X = resid8A20v$V4,zsignc)
resid8A20v$c42 <- sapply(X = resid8A20v$V12,zsignc)
resid8A20v$c43 <- sapply(X = resid8A20v$V20,zsignc)
resid8A20v$c44 <- sapply(X = resid8A20v$V28,zsignc)
resid8A20v$c45 <- sapply(X = resid8A20v$V36,zsignc)
resid8A20v$c46 <- sapply(X = resid8A20v$V44,zsignc)
resid8A20v$c47 <- sapply(X = resid8A20v$V52,zsignc)
resid8A20v$c48 <- sapply(X = resid8A20v$V60,zsignc)


resid8A20v$c51 <- sapply(X = resid8A20v$V5,zsignc)
resid8A20v$c52 <- sapply(X = resid8A20v$V13,zsignc)
resid8A20v$c53 <- sapply(X = resid8A20v$V21,zsignc)
resid8A20v$c54 <- sapply(X = resid8A20v$V29,zsignc)
resid8A20v$c55 <- sapply(X = resid8A20v$V37,zsignc)
resid8A20v$c56 <- sapply(X = resid8A20v$V45,zsignc)
resid8A20v$c57 <- sapply(X = resid8A20v$V53,zsignc)
resid8A20v$c58 <- sapply(X = resid8A20v$V61,zsignc)


resid8A20v$c61 <- sapply(X = resid8A20v$V6,zsignc)
resid8A20v$c62 <- sapply(X = resid8A20v$V14,zsignc)
resid8A20v$c63 <- sapply(X = resid8A20v$V22,zsignc)
resid8A20v$c64 <- sapply(X = resid8A20v$V30,zsignc)
resid8A20v$c65 <- sapply(X = resid8A20v$V38,zsignc)
resid8A20v$c66 <- sapply(X = resid8A20v$V46,zsignc)
resid8A20v$c67 <- sapply(X = resid8A20v$V54,zsignc)
resid8A20v$c68 <- sapply(X = resid8A20v$V62,zsignc)


resid8A20v$c71 <- sapply(X = resid8A20v$V7,zsignc)
resid8A20v$c72 <- sapply(X = resid8A20v$V15,zsignc)
resid8A20v$c73 <- sapply(X = resid8A20v$V23,zsignc)
resid8A20v$c74 <- sapply(X = resid8A20v$V31,zsignc)
resid8A20v$c75 <- sapply(X = resid8A20v$V39,zsignc)
resid8A20v$c76 <- sapply(X = resid8A20v$V47,zsignc)
resid8A20v$c77 <- sapply(X = resid8A20v$V55,zsignc)
resid8A20v$c78 <- sapply(X = resid8A20v$V63,zsignc)


resid8A20v$c81 <- sapply(X = resid8A20v$V8,zsignc)
resid8A20v$c82 <- sapply(X = resid8A20v$V16,zsignc)
resid8A20v$c83 <- sapply(X = resid8A20v$V24,zsignc)
resid8A20v$c84 <- sapply(X = resid8A20v$V32,zsignc)
resid8A20v$c85 <- sapply(X = resid8A20v$V40,zsignc)
resid8A20v$c86 <- sapply(X = resid8A20v$V48,zsignc)
resid8A20v$c87 <- sapply(X = resid8A20v$V56,zsignc)
resid8A20v$c88 <- sapply(X = resid8A20v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8A20<- resid8A20v[,66:129]

datosLong8A20<-pivot_longer(zsig8A20, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8A20$N <- gl(1,320000, labels = 20)
datosLong8A20$Categoria<- rep(nombres8,5000)

tabla8A20<-table(datosLong8A20$Categoria,datosLong8A20$Medida)
tabla8A20p<- as.data.frame(prop.table(tabla8A20,1))
colnames(tabla8A20p)<- c("Patrón", "Dependencia", "P")
tabla8A20p$Cat<- rep(8,192)
tabla8A20p$Matriz<- rep("A",192)
tabla8A20p$N<- rep(20,192)

####

## N = 50, matriz A

transiz8A50<-alply(seq8A50,1,transi)
resid8A50<- llply(transiz8A50, function(x) resid2(x))

dim8<- ldply(resid8A50,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8A50 <- resid8A50[-falsos]
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



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8A50<- resid8A50v[,66:129]

datosLong8A50<-pivot_longer(zsig8A50, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8A50$N <- gl(1,320000, labels = 50)
datosLong8A50$Categoria<- rep(nombres8,5000)

tabla8A50<-table(datosLong8A50$Categoria,datosLong8A50$Medida)
tabla8A50p<- as.data.frame(prop.table(tabla8A50,1))
colnames(tabla8A50p)<- c("Patrón", "Dependencia", "P")
tabla8A50p$Cat<- rep(8,192)
tabla8A50p$Matriz<- rep("A",192)
tabla8A50p$N<- rep(50,192)

####

## N = 100, matriz A

transiz8A100<-alply(seq8A100,1,transi)
resid8A100<- llply(transiz8A100, function(x) resid2(x))

dim8<- ldply(resid8A100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

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



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8A100<- resid8A100v[,66:129]

datosLong8A100<-pivot_longer(zsig8A100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8A100$N <- gl(1,320000, labels = 100)
datosLong8A100$Categoria<- rep(nombres8,5000)

tabla8A100<-table(datosLong8A100$Categoria,datosLong8A100$Medida)
tabla8A100p<- as.data.frame(prop.table(tabla8A100,1))
colnames(tabla8A100p)<- c("Patrón", "Dependencia", "P")
tabla8A100p$Cat<- rep(8,192)
tabla8A100p$Matriz<- rep("A",192)
tabla8A100p$N<- rep(100,192)


## N = 20, matriz B

transiz8B20<-alply(seq8B20,1,transi)
resid8B20<- llply(transiz8B20, function(x) resid2(x))

dim8<- ldply(resid8B20,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8B20 <- resid8B20[-falsos]
resid8B20v<- ldply(resid8B20,as.vector)
buenos<- complete.cases(resid8B20v)
resid8B20v<- resid8B20v[buenos,]
resid8B20v<- resid8B20v[1:5000,]


resid8B20v$c11 <- sapply(X = resid8B20v$V1,zsignc)
resid8B20v$c12 <- sapply(X = resid8B20v$V9,zsignc)
resid8B20v$c13 <- sapply(X = resid8B20v$V17,zsignc)
resid8B20v$c14 <- sapply(X = resid8B20v$V25,zsignc)
resid8B20v$c15 <- sapply(X = resid8B20v$V33,zsignc)
resid8B20v$c16 <- sapply(X = resid8B20v$V41,zsignc)
resid8B20v$c17 <- sapply(X = resid8B20v$V49,zsignc)
resid8B20v$c18 <- sapply(X = resid8B20v$V57,zsignc)


resid8B20v$c21 <- sapply(X = resid8B20v$V2,zsignc)
resid8B20v$c22 <- sapply(X = resid8B20v$V10,zsignc)
resid8B20v$c23 <- sapply(X = resid8B20v$V18,zsignc)
resid8B20v$c24 <- sapply(X = resid8B20v$V26,zsignc)
resid8B20v$c25 <- sapply(X = resid8B20v$V34,zsignc)
resid8B20v$c26 <- sapply(X = resid8B20v$V42,zsignc)
resid8B20v$c27 <- sapply(X = resid8B20v$V50,zsignc)
resid8B20v$c28 <- sapply(X = resid8B20v$V58,zsignc)


resid8B20v$c31 <- sapply(X = resid8B20v$V3,zsignc)
resid8B20v$c32 <- sapply(X = resid8B20v$V11,zsignc)
resid8B20v$c33 <- sapply(X = resid8B20v$V19,zsignc)
resid8B20v$c34 <- sapply(X = resid8B20v$V27,zsignc)
resid8B20v$c35 <- sapply(X = resid8B20v$V35,zsignc)
resid8B20v$c36 <- sapply(X = resid8B20v$V43,zsignc)
resid8B20v$c37 <- sapply(X = resid8B20v$V51,zsignc)
resid8B20v$c38 <- sapply(X = resid8B20v$V59,zsignc)


resid8B20v$c41 <- sapply(X = resid8B20v$V4,zsignc)
resid8B20v$c42 <- sapply(X = resid8B20v$V12,zsignc)
resid8B20v$c43 <- sapply(X = resid8B20v$V20,zsignc)
resid8B20v$c44 <- sapply(X = resid8B20v$V28,zsignc)
resid8B20v$c45 <- sapply(X = resid8B20v$V36,zsignc)
resid8B20v$c46 <- sapply(X = resid8B20v$V44,zsignc)
resid8B20v$c47 <- sapply(X = resid8B20v$V52,zsignc)
resid8B20v$c48 <- sapply(X = resid8B20v$V60,zsignc)


resid8B20v$c51 <- sapply(X = resid8B20v$V5,zsignc)
resid8B20v$c52 <- sapply(X = resid8B20v$V13,zsignc)
resid8B20v$c53 <- sapply(X = resid8B20v$V21,zsignc)
resid8B20v$c54 <- sapply(X = resid8B20v$V29,zsignc)
resid8B20v$c55 <- sapply(X = resid8B20v$V37,zsignc)
resid8B20v$c56 <- sapply(X = resid8B20v$V45,zsignc)
resid8B20v$c57 <- sapply(X = resid8B20v$V53,zsignc)
resid8B20v$c58 <- sapply(X = resid8B20v$V61,zsignc)


resid8B20v$c61 <- sapply(X = resid8B20v$V6,zsignc)
resid8B20v$c62 <- sapply(X = resid8B20v$V14,zsignc)
resid8B20v$c63 <- sapply(X = resid8B20v$V22,zsignc)
resid8B20v$c64 <- sapply(X = resid8B20v$V30,zsignc)
resid8B20v$c65 <- sapply(X = resid8B20v$V38,zsignc)
resid8B20v$c66 <- sapply(X = resid8B20v$V46,zsignc)
resid8B20v$c67 <- sapply(X = resid8B20v$V54,zsignc)
resid8B20v$c68 <- sapply(X = resid8B20v$V62,zsignc)


resid8B20v$c71 <- sapply(X = resid8B20v$V7,zsignc)
resid8B20v$c72 <- sapply(X = resid8B20v$V15,zsignc)
resid8B20v$c73 <- sapply(X = resid8B20v$V23,zsignc)
resid8B20v$c74 <- sapply(X = resid8B20v$V31,zsignc)
resid8B20v$c75 <- sapply(X = resid8B20v$V39,zsignc)
resid8B20v$c76 <- sapply(X = resid8B20v$V47,zsignc)
resid8B20v$c77 <- sapply(X = resid8B20v$V55,zsignc)
resid8B20v$c78 <- sapply(X = resid8B20v$V63,zsignc)


resid8B20v$c81 <- sapply(X = resid8B20v$V8,zsignc)
resid8B20v$c82 <- sapply(X = resid8B20v$V16,zsignc)
resid8B20v$c83 <- sapply(X = resid8B20v$V24,zsignc)
resid8B20v$c84 <- sapply(X = resid8B20v$V32,zsignc)
resid8B20v$c85 <- sapply(X = resid8B20v$V40,zsignc)
resid8B20v$c86 <- sapply(X = resid8B20v$V48,zsignc)
resid8B20v$c87 <- sapply(X = resid8B20v$V56,zsignc)
resid8B20v$c88 <- sapply(X = resid8B20v$V64,zsignc)



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8B20<- resid8B20v[,66:129]

datosLong8B20<-pivot_longer(zsig8B20, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8B20$N <- gl(1,320000, labels = 20)
datosLong8B20$Categoria<- rep(nombres8,5000)

tabla8B20<-table(datosLong8B20$Categoria,datosLong8B20$Medida)
tabla8B20p<- as.data.frame(prop.table(tabla8B20,1))
colnames(tabla8B20p)<- c("Patrón", "Dependencia", "P")
tabla8B20p$Cat<- rep(8,192)
tabla8B20p$Matriz<- rep("B",192)
tabla8B20p$N<- rep(20,192)

####

## N = 50, matriz B

transiz8B50<-alply(seq8B50,1,transi)
resid8B50<- llply(transiz8B50, function(x) resid2(x))

dim8<- ldply(resid8B50,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8B50<-resid8B50[-falsos]
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



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8B50<- resid8B50v[,66:129]

datosLong8B50<-pivot_longer(zsig8B50, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8B50$N <- gl(1,320000, labels = 50)

datosLong8B50$Categoria<- rep(nombres8,5000)

tabla8B50<-table(datosLong8B50$Categoria,datosLong8B50$Medida)
tabla8B50p<- as.data.frame(prop.table(tabla8B50,1))
colnames(tabla8B50p)<- c("Patrón", "Dependencia", "P")
tabla8B50p$Cat<- rep(8,192)
tabla8B50p$Matriz<- rep("B",192)
tabla8B50p$N<- rep(50,192)

####


## N = 100, matriz B

transiz8B100<-alply(seq8B100,1,transi)
resid8B100<- llply(transiz8B100, function(x) resid2(x))

dim8<- ldply(resid8B100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)


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



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8B100<- resid8B100v[,66:129]

datosLong8B100<-pivot_longer(zsig8B100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong8B100$N <- gl(1,320000, labels = 100)

datosLong8B100$Categoria<- rep(nombres8,5000)

tabla8B100<-table(datosLong8B100$Categoria,datosLong8B100$Medida)
tabla8B100p<- as.data.frame(prop.table(tabla8B100,1))
colnames(tabla8B100p)<- c("Patrón", "Dependencia", "P")
tabla8B100p$Cat<- rep(8,192)
tabla8B100p$Matriz<- rep("B",192)
tabla8B100p$N<- rep(100,192)


####



## N = 20, matriz C

transiz8C20<-alply(seq8C20,1,transi)
resid8C20<- llply(transiz8C20, function(x) resid2(x))

dim8<- ldply(resid8C20,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8C20 <- resid8C20[-falsos]
resid8C20v<- ldply(resid8C20,as.vector)
buenos<- complete.cases(resid8C20v)
resid8C20v<- resid8C20v[buenos,]
resid8C20v<- resid8C20v[1:5000,]


resid8C20v$c11 <- sapply(X = resid8C20v$V1,zsignc)
resid8C20v$c12 <- sapply(X = resid8C20v$V9,zsignc)
resid8C20v$c13 <- sapply(X = resid8C20v$V17,zsignc)
resid8C20v$c14 <- sapply(X = resid8C20v$V25,zsignc)
resid8C20v$c15 <- sapply(X = resid8C20v$V33,zsignc)
resid8C20v$c16 <- sapply(X = resid8C20v$V41,zsignc)
resid8C20v$c17 <- sapply(X = resid8C20v$V49,zsignc)
resid8C20v$c18 <- sapply(X = resid8C20v$V57,zsignc)


resid8C20v$c21 <- sapply(X = resid8C20v$V2,zsignc)
resid8C20v$c22 <- sapply(X = resid8C20v$V10,zsignc)
resid8C20v$c23 <- sapply(X = resid8C20v$V18,zsignc)
resid8C20v$c24 <- sapply(X = resid8C20v$V26,zsignc)
resid8C20v$c25 <- sapply(X = resid8C20v$V34,zsignc)
resid8C20v$c26 <- sapply(X = resid8C20v$V42,zsignc)
resid8C20v$c27 <- sapply(X = resid8C20v$V50,zsignc)
resid8C20v$c28 <- sapply(X = resid8C20v$V58,zsignc)


resid8C20v$c81 <- sapply(X = resid8C20v$V8,zsignc)
resid8C20v$c82 <- sapply(X = resid8C20v$V16,zsignc)
resid8C20v$c83 <- sapply(X = resid8C20v$V24,zsignc)
resid8C20v$c84 <- sapply(X = resid8C20v$V32,zsignc)
resid8C20v$c85 <- sapply(X = resid8C20v$V40,zsignc)
resid8C20v$c86 <- sapply(X = resid8C20v$V48,zsignc)
resid8C20v$c87 <- sapply(X = resid8C20v$V56,zsignc)
resid8C20v$c88 <- sapply(X = resid8C20v$V64,zsignc)

resid8C20v$c31 <- sapply(X = resid8C20v$V3,zsignc)
resid8C20v$c32 <- sapply(X = resid8C20v$V11,zsignc)
resid8C20v$c33 <- sapply(X = resid8C20v$V19,zsignc)
resid8C20v$c34 <- sapply(X = resid8C20v$V27,zsignc)
resid8C20v$c35 <- sapply(X = resid8C20v$V35,zsignc)
resid8C20v$c36 <- sapply(X = resid8C20v$V43,zsignc)
resid8C20v$c37 <- sapply(X = resid8C20v$V51,zsignc)
resid8C20v$c38 <- sapply(X = resid8C20v$V59,zsignc)


resid8C20v$c41 <- sapply(X = resid8C20v$V4,zsignc)
resid8C20v$c42 <- sapply(X = resid8C20v$V12,zsignc)
resid8C20v$c43 <- sapply(X = resid8C20v$V20,zsignc)
resid8C20v$c44 <- sapply(X = resid8C20v$V28,zsignc)
resid8C20v$c45 <- sapply(X = resid8C20v$V36,zsignc)
resid8C20v$c46 <- sapply(X = resid8C20v$V44,zsignc)
resid8C20v$c47 <- sapply(X = resid8C20v$V52,zsignc)
resid8C20v$c48 <- sapply(X = resid8C20v$V60,zsignc)


resid8C20v$c51 <- sapply(X = resid8C20v$V5,zsignc)
resid8C20v$c52 <- sapply(X = resid8C20v$V13,zsignc)
resid8C20v$c53 <- sapply(X = resid8C20v$V21,zsignc)
resid8C20v$c54 <- sapply(X = resid8C20v$V29,zsignc)
resid8C20v$c55 <- sapply(X = resid8C20v$V37,zsignc)
resid8C20v$c56 <- sapply(X = resid8C20v$V45,zsignc)
resid8C20v$c57 <- sapply(X = resid8C20v$V53,zsignc)
resid8C20v$c58 <- sapply(X = resid8C20v$V61,zsignc)


resid8C20v$c61 <- sapply(X = resid8C20v$V6,zsignc)
resid8C20v$c62 <- sapply(X = resid8C20v$V14,zsignc)
resid8C20v$c63 <- sapply(X = resid8C20v$V22,zsignc)
resid8C20v$c64 <- sapply(X = resid8C20v$V30,zsignc)
resid8C20v$c65 <- sapply(X = resid8C20v$V38,zsignc)
resid8C20v$c66 <- sapply(X = resid8C20v$V46,zsignc)
resid8C20v$c67 <- sapply(X = resid8C20v$V54,zsignc)
resid8C20v$c68 <- sapply(X = resid8C20v$V62,zsignc)


resid8C20v$c71 <- sapply(X = resid8C20v$V7,zsignc)
resid8C20v$c72 <- sapply(X = resid8C20v$V15,zsignc)
resid8C20v$c73 <- sapply(X = resid8C20v$V23,zsignc)
resid8C20v$c74 <- sapply(X = resid8C20v$V31,zsignc)
resid8C20v$c75 <- sapply(X = resid8C20v$V39,zsignc)
resid8C20v$c76 <- sapply(X = resid8C20v$V47,zsignc)
resid8C20v$c77 <- sapply(X = resid8C20v$V55,zsignc)
resid8C20v$c78 <- sapply(X = resid8C20v$V63,zsignc)


resid8C20v$c81 <- sapply(X = resid8C20v$V8,zsignc)
resid8C20v$c82 <- sapply(X = resid8C20v$V16,zsignc)
resid8C20v$c83 <- sapply(X = resid8C20v$V24,zsignc)
resid8C20v$c84 <- sapply(X = resid8C20v$V32,zsignc)
resid8C20v$c85 <- sapply(X = resid8C20v$V40,zsignc)
resid8C20v$c86 <- sapply(X = resid8C20v$V48,zsignc)
resid8C20v$c87 <- sapply(X = resid8C20v$V56,zsignc)
resid8C20v$c88 <- sapply(X = resid8C20v$V64,zsignc)




# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8C20<- resid8C20v[,66:129]

datosLong8C20<-pivot_longer(zsig8C20, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8C20$N <- gl(1,320000, labels = 20)

datosLong8C20$Categoria<- rep(nombres8,5000)

tabla8C20<-table(datosLong8C20$Categoria,datosLong8C20$Medida)
tabla8C20p<- as.data.frame(prop.table(tabla8C20,1))
colnames(tabla8C20p)<- c("Patrón", "Dependencia", "P")
tabla8C20p$Cat<- rep(8,192)
tabla8C20p$Matriz<- rep("C",192)
tabla8C20p$N<- rep(20,192)


####

## N = 50, matriz C

transiz8C50<-alply(seq8C50,1,transi)
resid8C50<- llply(transiz8C50, function(x) resid2(x))

dim8<- ldply(resid8C50,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8C50 <- resid8C50[-falsos]
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



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8C50<- resid8C50v[,66:129]

datosLong8C50<-pivot_longer(zsig8C50, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                            names_to = c("N", "Categoria"))

datosLong8C50$N <- gl(1,320000, labels = 50)

datosLong8C50$Categoria<- rep(nombres8,5000)

tabla8C50<-table(datosLong8C50$Categoria,datosLong8C50$Medida)
tabla8C50p<- as.data.frame(prop.table(tabla8C50,1))
colnames(tabla8C50p)<- c("Patrón", "Dependencia", "P")
tabla8C50p$Cat<- rep(8,192)
tabla8C50p$Matriz<- rep("C",192)
tabla8C50p$N<- rep(50,192)

####

## N = 100, matriz C

transiz8C100<-alply(seq8C100,1,transi)
resid8C100<- llply(transiz8C100, function(x) resid2(x))

dim8<- ldply(resid8C100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 !=8)
length(falsos)

resid8C100 <- resid8C100[-falsos]
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



# Cálculo de los patrones  de activación, inhibición e independencia 

zsig8C100<- resid8C100v[,66:129]

datosLong8C100<-pivot_longer(zsig8C100, cols = 1:64,values_to = "Medida", names_pattern = "(.*)_(.*)",
                             names_to = c("N", "Categoria"))

datosLong8C100$N <- gl(1,320000, labels = 100)

datosLong8C100$Categoria<- rep(nombres8,5000)

tabla8C100<-table(datosLong8C100$Categoria,datosLong8C100$Medida)
tabla8C100p<- as.data.frame(prop.table(tabla8C100,1))
colnames(tabla8C100p)<- c("Patrón", "Dependencia", "P")
tabla8C100p$Cat<- rep(8,192)
tabla8C100p$Matriz<- rep("C",192)
tabla8C100p$N<- rep(100,192)

# Tabla resumen de datos


tabla8<- rbind(tabla8A20p,tabla8A50p,tabla8A100p,tabla8B20p,tabla8B50p,tabla8B100p, 
               tabla8C20p,tabla8C50p,tabla8C100p)

patron_BC8<- tabla8[tabla8$Patrón== "BC",]


# Representación gráfica

plot8<-ggplot(patron_BC8,aes(x = factor(N), y = P,group = Dependencia,shape = Dependencia))+
  geom_line(aes(linetype = Dependencia),size = 1)+
  geom_point(size =2)+
  facet_grid(~Matriz)+
  theme_gray()

print(plot8)


tab8b<- pivot_wider(tabla8, names_from = c(N,Dependencia), values_from = P)
write.csv2(tab8b, "Tabla8.csv")

