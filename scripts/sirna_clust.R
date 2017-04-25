# siRNA data mining

# Import libraries  ====
library(FactoMineR)
library(mixOmics)

# Import data  ====
source("scripts/siRNA_hitselection.R")

# Hit classification
#samp2$hit = ifelse(samp2$z>1.5 | samp2$z<(-1.5), T,F)
#samp2$hit = ifelse(samp2$z.avg>1.5 | samp2$z.avg<(-1.5), T,F)
#samp2$hit = ifelse(scale(samp2$dead.mddc)>1.5 | scale(samp2$dead.mddc)<(-1.5), T,F)
samp2$hit = ifelse(scale(samp2$p24.supt1)>1.5 | scale(samp2$p24.supt1)<(-1.5), 1,0)

#z = samp2$z
z= scale(samp2$p24.supt1)

type.hits=factor(c("No","Low","High"))
for(i in 1:nrow(samp2)){
  type.hits[i]="No"
  if (z[i]<(-1.5)){
    type.hits[i]="Low"
  }
  if (z[i]>1.5){
    type.hits[i]="High"
  }
}

# PCA ====
data.pca=samp2[,c(2:5)]
data.pca=cbind(data.pca,type.hits)

pca = PCA(data.pca, quali.sup = 5)
plot.PCA(pca,axes=c(1,2), habillage=5)
pca$eig
pca$var$coord
pca$var$contrib


# PCA without p24.supt1
data.pca=data.pca[,-1]

pca = PCA(data.pca, quali.sup = 4)
pca$eig
pca$var$coord
pca$var$contrib

plot.PCA(pca,axes=c(1,2), habillage=4)
clust=HCPC(pca, nb.clust=4)

# PLS-DA
plsda1=plsda(Y=data.pca$type.hits, X=samp2[,3:5])
plotIndiv(plsda1, ellipse=T)



