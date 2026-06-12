
# Loading libraries

library(plotrix)
library(DescTools)

# Function to calculate the entropy measure

entropia<- function(x){
  res<-  apply(x,1,Entropy)
  return(res)
}


# Lag frequency table

nombres<- c("H+","M+","H-","M-")
x<- c(13,6,7,9,9,25,5,6,6,7,6,14,7,8,14,46)
datos<- matrix(x,4,4,byrow=T, dimnames = list(nombres,nombres))

# Calculation of the entropy measure

entropia3<- entropia(datos)


# Calculation of expected frequencies and entropy for a random model

exp<- chisq.test(x)$expected

aleatorio1<- matrix(exp,4,4,byrow = T,dimnames = list(nombres,nombres))
aleatorio2<- matrix(11.75,4,4,byrow=T, dimnames = list(nombres,nombres))

entropia1<- entropia(aleatorio1)
entropia2<- entropia(aleatorio2)


resultado1<- scale(c(entropia1,entropia3))
resultado2<- scale(c(entropia2,entropia3))

# Plots

polar.plot(resultado1[1:4],clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

polar.plot(resultado1[5:8],clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)


polar.plot(resultado2[1:4],clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

polar.plot(resultado2[5:8],clocwise=T, rp.type = "p",label=nombres,
           label.pos = c(45,135,225,315),
           radial.lim=c(-4,4),start=90,lwd=3,line.col=1)

