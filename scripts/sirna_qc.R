##################################################
###                                            ###
###     siRNA Quality Control 15/11/16         ###
###                                            ###
##################################################

# Import libraries  ====
#source("http://bioconductor.org/biocLite.R")
library(stats)
library(ggplot2)
library(pROC)
library(FactoMineR)

# Import data
source("scripts/sirna_io.R")

# Raw data visualisation ====

ggplot(data = samp, aes(x=samp$wellindex, y=samp$p24.supt1))+
  ggtitle(label = "Raw transfer")+
  geom_point()

ggplot(data = samp, aes(x=samp$wellindex, y=samp$p24.supt1[order(samp$p24.supt1)]))+
  ggtitle(label = "Raw transfer")+
  geom_bar(stat = "identity", width=0.5)

# Normalisation to non-target siRNA ====

raw1$nt.norm = ((raw1$p24.supt1/raw1.c[which(raw1.c$sample=="Control_NT.fcs"),"p24.supt1"])-1)*100
raw2$nt.norm = ((raw2$p24.supt1/raw2.c[which(raw2.c$sample=="Control_NT.fcs"),"p24.supt1"])-1)*100
raw3$nt.norm = ((raw3$p24.supt1/raw3.c[which(raw3.c$sample=="Control_NT.fcs"),"p24.supt1"])-1)*100
raw4$nt.norm = ((raw4$p24.supt1/raw4.c[which(raw4.c$sample=="Control_NT.fcs"),"p24.supt1"])-1)*100

ggplot(data = raw1, aes(x=raw1$wellindex, y=raw1$nt.norm[order(raw1$nt.norm)]))+
  ggtitle(label = "Plate 1 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))

ggplot(data = raw2, aes(x=raw2$wellindex, y=raw2$nt.norm[order(raw2$nt.norm)]))+
  ggtitle(label = "Plate 2 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))

ggplot(data = raw3, aes(x=raw3$wellindex, y=raw3$nt.norm[order(raw3$nt.norm)]))+
  ggtitle(label = "Plate 3 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))

ggplot(data = raw4, aes(x=raw4$wellindex, y=raw4$nt.norm[order(raw4$nt.norm)]))+
  ggtitle(label = "Plate 4 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))


samp$nt.norm = c(raw1$nt.norm, raw2$nt.norm, raw3$nt.norm, raw4$nt.norm)

ggplot(data = samp, aes(x=samp$wellindex, y=samp$nt.norm[order(samp$nt.norm)]))+
  ggtitle(label = "Screen 1 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))


# Normalisation: z-score (x-u/o) ====

raw1$z = scale(raw1$p24.supt1)
raw2$z = scale(raw2$p24.supt1)
raw3$z = scale(raw3$p24.supt1)
raw4$z = scale(raw4$p24.supt1)

samp = rbind.fill(raw1,raw2, raw3, raw4)
summary(samp$z)
samp = na.omit(samp)
attach(samp)

# Quality metrics ====
# Z-factor

zfactor = function(sigct,sigcb,mu.cb,mu.ct){
  z1=1-(3*sigct+3*sigcb)/abs(mu.ct-mu.cb)
  return(z1=z1)
}
