### siRNA hit selection

# Import libraries  ====
#source("http://bioconductor.org/biocLite.R")
library(pROC)

# Import data ====
source("scripts/sirna_qc.R")

# Mean +/- k-SD ====
quantile(samp$z, c(0.05, 0.95), type = 4)

ggplot(data = samp, aes(x=samp$wellindex, y=samp$z))+
  geom_point()+
  ggtitle(label = "Screen 1: k-standard deviations hit identification")+
  geom_hline(yintercept = c(-1.321557,1.67679))+
  geom_text(aes(label=ifelse(samp$z>1.67679,as.character(samp$genesymbol),'')),hjust=-0.1,vjust=-0.1)+
  geom_text(aes(label=ifelse(samp$z<(-1.321557),as.character(samp$genesymbol),'')),hjust=-0.1,vjust=-0.1)

ggplot(data = samp2, aes(x=samp2$wellindex, y=samp2$z))+
  geom_point()+
  ggtitle(label = "Screen 2: k-standard deviations hit identification")+
  geom_hline(yintercept = c(-1.5,1.5))+
  geom_text(aes(label=ifelse(samp2$z>1.5,as.character(samp2$genesymbol),'')),hjust=-0.1,vjust=-0.1)+
  geom_text(aes(label=ifelse(samp2$z<(-1.5),as.character(samp2$genesymbol),'')),hjust=-0.1,vjust=-0.1)

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
samp$hit.ksd = as.vector(ifelse(samp$z>1.5 | samp$z<(-1.5), 1, 0))
samp$hit.zmad = as.numeric(ifelse(samp$zmad>1.98, 1, 0))


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


# Hit confirmation ====
avg = matrix(c(samp$z, samp2$z), 319, 2)
avg = cbind(avg, rowMeans(avg[,1:2]))


samp2$z.avg = avg[,3]

ggplot(data = samp2, aes(x=samp2$wellindex, y=samp2$z.avg))+
  geom_point()+
  ggtitle(label = "Screen 1&2 z.avg: k-standard deviations hit identification")+
  geom_hline(yintercept = c(-1.5,1.5))+
  geom_text(aes(label=ifelse(samp2$z.avg>1.5,as.character(samp2$genesymbol),'')),hjust=-0.1,vjust=-0.1)+
  geom_text(aes(label=ifelse(samp2$z.avg<(-1.5),as.character(samp2$genesymbol),'')),hjust=-0.1,vjust=-0.1)+
  geom_text(aes(label=ifelse(samp2$genesymbol=="IL8",as.character(samp2$genesymbol),'')),hjust=-0.1,vjust=-0.1)

# Output ====
write.table(samp$accession[which(samp$nt.norm > 50)], "output/ntup.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$nt.norm <(-20))], "output/ntdown.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$z>1.5)], "output/ksdup.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$z<(-1.5))], "output/ksddown.txt", quote=F, row.names = F, col.names = F)
write.table(samp$accession[which(samp$zmad>1.98)], "output/kmad.txt", quote=F, row.names = F, col.names = F)

#cor1=cor(samp$p24.supt1[which(samp$plateno=="plate1")], samp$p24.supt1[which(samp$plateno=="plate2")])
#plot(samp$p24.supt1[which(samp$plateno=="plate1")], samp$p24.supt1[which(samp$plateno=="plate2")])
