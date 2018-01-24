### siRNA file i/o

# Import libraries ====
if(!require(plyr)){
  install.packages(c("plyr"))
  require(plyr)
}

if(!require(dplyr)){
  install.packages(c("dplyr"))
  require(dplyr)
}

if(!require(FactoMineR)){
  install.packages(c("FactoMineR"))
  require(FactoMineR)
}

if(!require(pROC)){
  install.packages(c("pROC"))
  require(pROC)
}

if(!require(mixOmics)){
  install.packages(c("mixOmics"))
  require(mixOmics)
}

# Import siRNA cytokine receptor library ====
cclib = read.csv("D:/My PhD/siRNA screen/Libraries/G-104000 OTP Human Cytokine Receptors Lot 060724.csv")
cclib = `colnames<-`(cclib, c("plate","well","catalog1","catalog2","genesymbol", "geneid", "accession", "ginumber"))
cclib = na.omit(cclib)

# Data import ====
import = function(path){
  x = read.csv(path, sep=",")
  x = `colnames<-`(x, c("sample","p24.supt1","dead.supt1", "p24.mddc", "dead.mddc")); x=x[,1:5]
  ind = which(with(x, x$sample == "Mean" | x$sample == "SD"))
  if("Mean" %in% x$sample) x = x[-ind,] else x=x
}

raw1 = import("data/180117 p1-1.csv")
raw2 = import("data/220117 p1-2.csv")
raw3 = import("data/290117 p1-3.csv")
raw4 = import("data/060217 p1-4.csv")


# Subset and clean ====
clean = function(x, plateno, wellstart){
  x = filter(x, !grepl("Control_*", sample))
  n = length(x[,1])-1
  x$plateno = rep(plateno, n+1)
  x$wellindex = wellstart:(wellstart+n)
  x$genesymbol = as.character(unique(cclib$genesymbol)[wellstart:(wellstart+n)])
  x$accession = as.character(unique(cclib$accession)[wellstart:(wellstart+n)])
   return(x)
}

raw1.c = filter(raw1, grepl("Control_*", raw1$sample))
raw2.c = filter(raw2, grepl("Control_*", raw2$sample))
raw3.c = filter(raw3, grepl("Control_*", raw3$sample))
raw4.c = filter(raw4, grepl("Control_*", raw4$sample))

raw1 = clean(raw1,1,1)
raw2 = clean(raw2,2,81)
raw3 = clean(raw3,3,161)
raw4 = clean(raw4,4,241)

# Collate data
samp = rbind.fill(raw1,raw2,raw3, raw4)
samp.c = rbind.fill(raw1.c, raw2.c, raw3.c, raw4.c)


# Screen 2 ====
raw5 = import("data/010317 p2-1.csv")
raw5.c = filter(raw5, grepl("Control_*", raw5$sample))
raw5 = clean(raw5,1,1)

raw6 = import("data/070317 p2-2.csv")
raw6.c = filter(raw6, grepl("Control_*", raw6$sample))
raw6 = clean(raw6,2,81)

raw7 = import("data/130317 p2-3.csv")
raw7.c = filter(raw7, grepl("Control_*", raw7$sample))
raw7 = clean(raw7,3,161)

raw8 = import("data/200317 p2-4.csv")
raw8.c = filter(raw8, grepl("Control_*", raw8$sample))
raw8 = clean(raw8,4,241)

samp2 = rbind.fill(raw5, raw6, raw7, raw8)
samp2.c = rbind.fill(raw5.c, raw6.c, raw7.c, raw8.c)

