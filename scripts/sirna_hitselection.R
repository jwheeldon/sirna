##################################################
###                                            ###
###     siRNA hit selection - 15/11/16         ###
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

# Mean +/- k-SD ====
ggplot(data = samp, aes(x=samp$wellindex, y=samp$z))+
  geom_point()+
  ggtitle(label = "k-standard deviations hit identification")+
  geom_hline(yintercept = c(-1.5,1.5))+
  geom_text(aes(label=ifelse(samp$z>1.5,as.character(samp$genesymbol),'')),hjust=-0.1,vjust=-0.1)+
  geom_text(aes(label=ifelse(samp$z<(-1.5),as.character(samp$genesymbol),'')),hjust=-0.1,vjust=-0.1)


# Median +/- k-MAD Hit Identification ====
median(samp$z)
samp$zmad = abs(samp$z- median(samp$z))

ggplot(data = samp, aes(x=samp$wellindex, y=samp$zmad, label=samp$genesymbol))+
  ggtitle(label = "k-Median absolute deviation hit identification")+
  geom_point()+
  geom_hline(yintercept = 2)+
  geom_text(aes(label=ifelse(samp$zmad>2,as.character(samp$genesymbol),'')),hjust=-0.1,vjust=-0.1)


# ROC curves ====

samp$hit.cut = ifelse(samp$nt.norm>50 | samp$nt.norm<(-20), 1, 0)
samp$hit.zmad = as.numeric(ifelse(samp$zmad>1.98, 1, 0))
samp$hit.ksd = as.vector(ifelse(samp$z>1.5 | samp$z<(-1.5), 1, 0))

model1=glm(samp$hit.cut ~ samp$p24.supt1 + samp$dead.supt1 + samp$p24.mddc + samp$dead.mddc,
           data=samp,
           family="binomial")
summary(model1)
predict(model1, type="response")
roc.cut=roc(samp$hit.cut, predict(model1, type="response"))
#roc(samp$hit.cut, runif(240,0,1), plot=T) #Coin flipping ROC curve
coords(roc.cut, x="best")


model2=glm(samp$hit.zmad ~ samp$p24.supt1 + samp$dead.supt1 + samp$p24.mddc + samp$dead.mddc,
           data=samp, 
           family="binomial")
summary(model2)
predict(model2, type="response")
roc.zmad=roc(samp$hit.zmad, predict(model2, type="response"))


model3=glm(samp$hit.ksd ~ samp$p24.supt1 + samp$dead.supt1 + samp$p24.mddc + samp$dead.mddc,
           data=samp, 
           family="binomial")
summary(model3)
predict(model3, type="response")
roc.ksd=roc(samp$hit.ksd, predict(model3, type="response"))


plot(roc.cut);plot(roc.zmad, add=T, col="red"); plot(roc.ksd, add=T, col="green")
legend(0.4,0.2, c("nt.cut 0.7656", "zmad 0.8056", "ksd 0.699"), lty=1, col=c('black','red','green'), bty='n', cex=.75)

# ROC Cross-validation: 'leave-one-out' ====
# Model 1 (nt.norm cut-off)
plot(roc.cut);plot(roc.zmad, add=T, col="red"); plot(roc.ksd, add=T, col="green")

prediction.loo = as.vector(samp$hit.cut)

for(i in 1:length(samp$hit.cut)){
  train = samp[-i,]
  test = samp[i,]
  model1 = glm(hit.cut ~ p24.supt1 + dead.supt1 + p24.mddc + dead.mddc, data = train, family = "binomial")
  prediction.loo[i] = predict(model1, newdata = test, type = "response")
}

roc.cut = roc(samp$hit.cut, prediction.loo)


# Model 2 (zMAD)
prediction2.loo = as.vector(samp$hit.zmad)

for(i in 1:length(samp$hit.zmad)){
  train = samp[-i,]
  test = samp[i,]
  model2 = glm(hit.zmad ~ p24.supt1 + dead.supt1 + p24.mddc + dead.mddc, data = train, family = "binomial")
  prediction2.loo[i] = predict(model2, newdata = test, type = "response")
}
roc.zmad = roc(samp$hit.zmad, prediction2.loo)



# Model 3 (ksd)
prediction3.loo = as.vector(samp$hit.ksd)

for(i in 1:length(samp$hit.ksd)){
  train = samp[-i,]
  test = samp[i,]
  model3 = glm(hit.ksd ~ p24.supt1 + dead.supt1 + p24.mddc + dead.mddc, data = train, family = "binomial")
  prediction3.loo[i] = predict(model3, newdata = test, type = "response")
}
roc.ksd = roc(as.vector(samp$hit.ksd), prediction3.loo)

plot(roc.cut, add=T, col="blue")
plot(roc.zmad, add=T, col="blue")
plot(roc.ksd, add=T, col="blue")

plot(roc.cut, main="Cross validation: Leave-one-out");plot(roc.zmad, add=T, col="red"); plot(roc.ksd, add=T, col="green")
legend(0.4,0.2, c("nt.cut 0.7504", "zmad 0.6102", "ksd 0.6293"), lty=1, col=c('black','red','green'), bty='n', cex=.75)



# Strictly standardised mean difference (SSMD) hit identification ====
summary(z)



# PCA ====
pca = PCA(samp[,2:5])
plot.PCA(pca, choix="ind", axes=c(1,2))
pca$eig
pca$var$coord
pca$var$contrib


# Output ====
write.table(samp$accession[which(samp$nt.norm > 50)], "output/ntup.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(raw1$nt.norm <(-20))], "output/ntdown.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$z>1.5)], "output/ksdup.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$z<1.5)], "output/ksddown.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$zmad>1.98)], "output/kmad.txt", quote=F, row.names = F, col.names = F)


#cor1=cor(samp$p24.supt1[which(samp$plateno=="plate1")], samp$p24.supt1[which(samp$plateno=="plate2")])
#plot(samp$p24.supt1[which(samp$plateno=="plate1")], samp$p24.supt1[which(samp$plateno=="plate2")])
detach()
