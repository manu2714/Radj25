start<- Sys.time()

library(simstudy)
library(MCTM)
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(misty)
library(stats)
library(DescTools)
library(plotrix)

set.seed(123456)


## MCTM library is removed from R. 
## Formerly available versions can be obtained from the archive.

# Functions 

## Calculation of transition frequencies

transi<- function(x) {
  TransMatrix(as.numeric(x),order= 1,probs= F)
}



entropia2<- function(x){
  apply(x,1,Entropy)
}

  # Definition of matrices with 8 states


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

# nombres4<- c("AA","AB","AC","AD",
#              "BA","BB","BC","BD",
#              "CA","CB","CC","CD",
#              "DA","DB","DC","DD")


# Generation of the sequences

## Sequences with 4 states


seq4A100 <- genMarkov(n = 500, transMat = m4A, 
                      chainLen = 100, wide = TRUE)
seq4B100 <- genMarkov(n = 5000, transMat = m4B, 
                      chainLen = 100, wide = TRUE)
seq4C100 <- genMarkov(n = 5000, transMat = m4C, 
                      chainLen = 100, wide = TRUE)

seq4A100<- seq4A100[,-1]
seq4B100<- seq4B100[,-1]
seq4C100<- seq4C100[,-1]


# Entropy of the 4-category condition

## N = 100, matrix A

transiz4A100<-alply(seq4A100,1,transi)

dim8<- ldply(transiz4A100,dim)
falsos<-which(dim8$V1 != 4 | dim8$V2 != 4)
length(falsos)

transiz4A100 <- transiz4A100
chisq4A100<- llply(transiz4A100,function(x) chisq.test(x)$statistic)
chi4A100<- llply(chisq4A100, as.vector)
chi24A100<- unlist(chi4A100)
buenos<- which(chi4A100 <= qchisq(.95,9))
length(buenos)

transiz4A100<- transiz4A100[buenos]


entro4A100<- llply(transiz4A100, entropia2)

entro4A100v<- ldply(entro4A100,as.vector)
entro4A100r<- apply(entro4A100v,1,Rank)
entro4A100rm<- apply(entro4A100v[,-1],2,mean,na.rm=T)


## N = 100, matrix B

transiz4B100<-alply(seq4B100,1,transi)

dim8<- ldply(transiz4B100,dim)
falsos<-which(dim8$V1 != 4 | dim8$V2 != 4)
length(falsos)


transiz4B100 <- transiz4B100
chisq4B100<- llply(transiz4B100,function(x) chisq.test(x)$statistic)
chi4B100<- llply(chisq4B100, as.vector)
chi24B100<- unlist(chi4B100)
buenos<- which(chi4B100 <= qchisq(.95,9))


entro4B100<- llply(transiz4B100, entropia2)

entro4B100v<- ldply(entro4B100,as.vector)
entro4B100r<- apply(entro4B100v,1,Rank)
entro4B100rm<- apply(entro4B100v[,-1],2,mean)


## N = 100, matrix C

transiz4C100<-alply(seq4C100,1,transi)

dim8<- ldply(transiz4C100,dim)
falsos<-which(dim8$V1 != 4 | dim8$V2 != 4)
length(falsos)


transiz4C100 <- transiz4C100
chisq4C100<- llply(transiz4C100,function(x) chisq.test(x)$statistic)
chi4C100<- llply(chisq4C100, as.vector)
chi24C100<- unlist(chi4C100)
buenos<- which(chi4C100 <= qchisq(.95,9))


entro4C100<- llply(transiz4C100, entropia2)

entro4C100v<- ldply(entro4C100,as.vector)
entro4C100r<- apply(entro4C100v,1,Rank)
entro4C100rm<- apply(entro4C100v[,-1],2,mean)

nombres<-  c("A","B","C","D")

polar.plot(entro4A100rm,clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

polar.plot(entro4B100rm,clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

polar.plot(entro4C100rm,clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)



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

# nombres6<- c("AA","AB","AC","AD","AE","AF",
#              "BA","BB","BC","BD","BE","BF",
#              "CA","CB","CC","CD","CE","CF",
#              "DA","DB","DC","DD","DE","DF",
#              "EA","EB","EC","ED","EE","EF",
#              "FA","FB","FC","FD","FE","FF")


# Generation of the sequences

## Sequences with 6 states



seq6A100 <- genMarkov(n = 5000, transMat = m6A, 
                      chainLen = 100, wide = TRUE)
seq6B100 <- genMarkov(n = 5000, transMat = m6B, 
                      chainLen = 100, wide = TRUE)
seq6C100 <- genMarkov(n = 5000, transMat = m6C, 
                      chainLen = 100, wide = TRUE)
                    

seq6A100<- seq6A100[,-1]
seq6B100<- seq6B100[,-1]
seq6C100<- seq6C100[,-1]



# Entropy of the 6-category condition


## N = 100, matrix A

transiz6A100<-alply(seq6A100,1,transi)

dim8<- ldply(transiz6A100,dim)
falsos<-which(dim8$V1 != 6 | dim8$V2 != 6)
length(falsos)

transiz6A100 <- transiz6A100
chisq6A100<- llply(transiz6A100,function(x) chisq.test(x)$statistic)
chi6A100<- llply(chisq6A100, as.vector)
chi26A100<- unlist(chi6A100)
buenos<- which(chi6A100 <= qchisq(.95,9))
length(buenos)

transiz6A100<- transiz6A100[buenos]


entro6A100<- llply(transiz6A100, entropia2)

entro6A100v<- ldply(entro6A100,as.vector)
entro6A100r<- apply(entro6A100v,1,Rank)
entro6A100rm<- apply(entro6A100v[,-1],2,mean,na.rm=T)


## N = 100, matrix B

transiz6B100<-alply(seq6B100,1,transi)

dim8<- ldply(transiz6B100,dim)
falsos<-which(dim8$V1 != 6 | dim8$V2 != 6)
length(falsos)


transiz6B100 <- transiz6B100
chisq6B100<- llply(transiz6B100,function(x) chisq.test(x)$statistic)
chi6B100<- llply(chisq6B100, as.vector)
chi26B100<- unlist(chi6B100)
buenos<- which(chi6B100 <= qchisq(.95,9))


entro6B100<- llply(transiz6B100, entropia2)

entro6B100v<- ldply(entro6B100,as.vector)
entro6B100r<- apply(entro6B100v,1,Rank)
entro6B100rm<- apply(entro6B100v[,-1],2,mean)


## N = 100, matrix C

transiz6C100<-alply(seq6C100,1,transi)

dim8<- ldply(transiz6C100,dim)
falsos<-which(dim8$V1 != 6 | dim8$V2 != 6)
length(falsos)


transiz6C100 <- transiz6C100
chisq6C100<- llply(transiz6C100,function(x) chisq.test(x)$statistic)
chi6C100<- llply(chisq6C100, as.vector)
chi26C100<- unlist(chi6C100)
buenos<- which(chi6C100 <= qchisq(.95,9))


entro6C100<- llply(transiz6C100, entropia2)

entro6C100v<- ldply(entro6C100,as.vector)
entro6C100r<- apply(entro6C100v,1,Rank)
entro6C100rm<- apply(entro6C100v[,-1],2,mean)

nombres6<- c("A","B","C","D","E","F")

posic<- c(0,60,60*2,60*3, 60*4, 60*5,60*6)+30

polar.plot(entro6A100rm,clocwise=F, rp.type = "p",label=nombres6,
           label.pos = posic,
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

polar.plot(entro6B100rm,clocwise=T, rp.type = "p",label= nombres6,
           label.pos = posic,
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

polar.plot(entro6C100rm,clocwise=T, rp.type = "p",label=nombres6,
           label.pos = posic,
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)




## Matrix   A

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


seq8A100 <- genMarkov(n = 5000, transMat = m8A, 
                      chainLen = 100, wide = TRUE)
seq8B100 <- genMarkov(n = 5000, transMat = m8B, 
                      chainLen = 100, wide = TRUE)
seq8C100 <- genMarkov(n = 5000, transMat = m8C, 
                      chainLen = 100, wide = TRUE)


seq8A100<- seq8A100[,-1]
seq8B100<- seq8B100[,-1]
seq8C100<- seq8C100[,-1]


transiz8A100<-alply(seq8A100,1,transi)

dim8<- ldply(transiz8A100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)

transiz8A100 <- transiz8A100
chisq8A100<- llply(transiz8A100,function(x) chisq.test(x)$statistic)
chi8A100<- llply(chisq8A100, as.vector)
chi28A100<- unlist(chi8A100)
buenos<- which(chi8A100 <= qchisq(.95,9))
length(buenos)

transiz8A100<- transiz8A100[buenos]


entro8A100<- llply(transiz8A100, entropia2)

entro8A100v<- ldply(entro8A100,as.vector)
entro8A100r<- apply(entro8A100v,1,Rank)
entro8A100rm<- apply(entro8A100v[,-1],2,mean,na.rm=T)


## N = 100, matrix B

transiz8B100<-alply(seq8B100,1,transi)

dim8<- ldply(transiz8B100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)


transiz8B100 <- transiz8B100
chisq8B100<- llply(transiz8B100,function(x) chisq.test(x)$statistic)
chi8B100<- llply(chisq8B100, as.vector)
chi28B100<- unlist(chi8B100)
buenos<- which(chi8B100 <= qchisq(.95,9))


entro8B100<- llply(transiz8B100, entropia2)

entro8B100v<- ldply(entro8B100,as.vector)
entro8B100r<- apply(entro8B100v,1,Rank)
entro8B100rm<- apply(entro8B100v[,-1],2,mean)


## N = 100, matrix C

transiz8C100<-alply(seq8C100,1,transi)

dim8<- ldply(transiz8C100,dim)
falsos<-which(dim8$V1 != 8 | dim8$V2 != 8)
length(falsos)


transiz8C100 <- transiz8C100
chisq8C100<- llply(transiz8C100,function(x) chisq.test(x)$statistic)
chi8C100<- llply(chisq8C100, as.vector)
chi28C100<- unlist(chi8C100)
buenos<- which(chi8C100 <= qchisq(.95,9))


entro8C100<- llply(transiz8C100, entropia2)

entro8C100v<- ldply(entro8C100,as.vector)
entro8C100r<- apply(entro8C100v,1,Rank)
entro8C100rm<- apply(entro8C100v[,-1],2,mean)

nombres8<- c("A","B","C","D","E","F","G","H")


posic<- c(0,45,45*2,45*3, 45*4, 45*5,45*6,45*7,45*8)


polar.plot(entro6A100rm,clocwise=T, rp.type = "p",label=nombres8,
           label.pos = posic,polar.pos= posic,
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)


polar.plot(entro6B100rm,clocwise=T, rp.type = "p",label= nombres8,
           label.pos = posic,polar.pos= posic,
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)


polar.plot(entro6C100rm,clocwise=T, rp.type = "p",label=nombres8,
           label.pos = posic,polar.pos= posic,
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)






end<- Sys.time()
print(end-start)
