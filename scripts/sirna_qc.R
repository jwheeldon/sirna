### siRNA Quality Control

# Import libraries  ====
#source("http://bioconductor.org/biocLite.R")

library(stats)
library(ggplot2)

# Import data
source("scripts/sirna_io.R")

# Raw data visualisation ====

ggplot(data = samp, aes(x=samp$wellindex, y=samp$p24.supt1))+
  ggtitle(label = "Screen 1: Raw transfer")+
  geom_point()

ggplot(data = samp2, aes(x=samp2$wellindex, y=samp2$p24.supt1))+
  ggtitle(label = "Screen 2: Raw transfer")+
  geom_point()


# Normalisation to non-target siRNA ====
raw1$nt.norm = ((raw1$p24.supt1/raw1.c[which(raw1.c$sample=="Control_NT.fcs"),"p24.supt1"])-1)*100
raw2$nt.norm = ((raw2$p24.supt1/raw2.c[which(raw2.c$sample=="Control_NT.fcs"),"p24.supt1"])-1)*100
raw3$nt.norm = ((raw3$p24.supt1/mean(raw3.c[grep("NT", raw3.c$sample),"p24.supt1"]))-1)*100
raw4$nt.norm = ((raw4$p24.supt1/mean(raw4.c[grep("NT", raw4.c$sample),"p24.supt1"]))-1)*100

raw5$nt.norm = ((raw5$p24.supt1/mean(raw5.c[grep("NT", raw5.c$sample),"p24.supt1"]))-1)*100
raw6$nt.norm = ((raw6$p24.supt1/mean(raw6.c[grep("NT", raw6.c$sample),"p24.supt1"]))-1)*100
raw7$nt.norm = ((raw7$p24.supt1/mean(raw7.c[grep("NT", raw7.c$sample),"p24.supt1"]))-1)*100
raw8$nt.norm = ((raw8$p24.supt1/mean(raw8.c[grep("NT", raw8.c$sample),"p24.supt1"]))-1)*100


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

ggplot(data = raw5, aes(x=raw5$wellindex, y=raw5$nt.norm[order(raw5$nt.norm)]))+
  ggtitle(label = "Plate 5 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))

ggplot(data = raw6, aes(x=raw6$wellindex, y=raw6$nt.norm[order(raw6$nt.norm)]))+
  ggtitle(label = "Plate 6 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))


samp$nt.norm = c(raw1$nt.norm, raw2$nt.norm, raw3$nt.norm, raw4$nt.norm)
samp2$nt.norm = c(raw5$nt.norm, raw6$nt.norm, raw7$nt.norm, raw8$nt.norm)

ggplot(data = samp, aes(x=samp$wellindex, y=samp$nt.norm[order(samp$nt.norm)]))+
  ggtitle(label = "Screen 1 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))

ggplot(data = samp2, aes(x=samp2$wellindex, y=samp2$nt.norm[order(samp2$nt.norm)]))+
  ggtitle(label = "Screen 2 Non-target normalised transfer")+
  geom_bar(stat = "identity", width=0.5)+
  geom_hline(yintercept = c(-20,50))


# Normalisation: z-score (x-u/o) ====
raw1$z = scale(raw1$p24.supt1)
raw2$z = scale(raw2$p24.supt1)
raw3$z = scale(raw3$p24.supt1)
raw4$z = scale(raw4$p24.supt1)

raw5$z = scale(raw5$p24.supt1)
raw6$z = scale(raw6$p24.supt1)
raw7$z = scale(raw7$p24.supt1)
raw8$z = scale(raw8$p24.supt1)

samp = rbind.fill(raw1,raw2, raw3, raw4)
summary(samp$z)
samp = na.omit(samp)

samp2 = rbind.fill(raw5,raw6, raw7, raw8)
samp2 = na.omit(samp2)

# Quality metrics: Z-factor ====
neg.cont = samp.c[which(samp.c$sample=="Control_F522Y.fcs"),]
pos.cont = samp.c[grep("CD4", samp.c$sample),]

zfactor = function(pos, neg){
  sig1 = sd(pos[,2])
  sig2 = sd(neg[,2])
  mu1 = mean(pos[,2])
  mu2 = mean(neg[,2])
  z1=1-(3*sig1+3*sig2)/abs(mu1-mu2)
  return(z1)
}
zfactor(pos.cont, neg.cont)
zfactor(samp, neg.cont)

# Quality metrics: SSMD ====
ssmd = function(pos,neg){
  mu1 = median(pos[,2])
  mu2 = median(neg[,2])
  sig1 = var(pos[,2])
  sig2 = var(neg[,2])
  ssmd= (mu1-mu2)/sqrt(sig1 + sig2)
  return(ssmd)
}
ssmd(pos.cont, neg.cont)
ssmd(samp, neg.cont)


