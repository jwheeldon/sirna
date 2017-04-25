# siRNA data mining
# 080317

# Import libraries  ====
library(FactoMineR)

# Import data  ====
source("scripts/siRNA_hitselection.R")

# Hit classification
#samp2$hit = ifelse(samp2$z>1.5 | samp2$z<(-1.5), T,F)
#samp2$hit = ifelse(samp2$z.avg>1.5 | samp2$z.avg<(-1.5), T,F)
#samp2$hit = ifelse(scale(samp2$dead.mddc)>1.5 | scale(samp2$dead.mddc)<(-1.5), T,F)
samp2$hit = ifelse(scale(samp2$p24.supt1)>1.5 | scale(samp2$p24.supt1)<(-1.5), 1,0)
z = samp2$z.avg

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
plot.PCA(pca,axes=c(1,2), habillage=4)


clust=HCPC(pca, nb.clust=4)
