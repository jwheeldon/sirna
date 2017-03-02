##################################################
###                                            ###
###         siRNA file i/o - 15/11/16          ###
###                                            ###
##################################################

# Import libraries ====
library(plyr)
library(dplyr)


# Import siRNA cytokine receptor library ====
cclib = read.csv("H:/My Documents/My PhD/siRNA screen/Libraries/G-104000 OTP Human Cytokine Receptors Lot 060724.csv")
cclib = `colnames<-`(cclib, c("plate","well","catalog1","catalog2","genesymbol", "geneid", "accession", "ginumber"))
cclib = na.omit(cclib)

# Data import function ====
import = function(path){
  x = read.csv(path, sep=",")
  x = `colnames<-`(x, c("sample","p24.supt1","dead.supt1", "p24.mddc", "dead.mddc")); x=x[,1:5]
  ind = which(with(x, x$sample == "Mean" | x$sample == "SD"))
  if("Mean" %in% x$sample) x = x[-ind,] else x=x
}

# Subset and clean function ====
clean = function(x, plateno, wellstart){
  x = filter(x, !grepl("Control_*", sample))
  n = length(x[,1])-1
  x$plateno = rep(plateno, n+1)
  x$wellindex = wellstart:(wellstart+n)
  x$genesymbol = as.character(unique(cclib$genesymbol)[wellstart:(wellstart+n)])
  x$accession = as.character(unique(cclib$accession)[wellstart:(wellstart+n)])
   return(x)
}

#Input data
raw1 = import("data/180117 p1-1.csv")
raw2 = import("data/220117 p1-2.csv")
raw3 = import("data/290117 p1-3.csv")
raw4 = import("data/060217 p1-4.csv")

raw1.c = filter(raw1, grepl("Control_*", raw1$sample))
raw2.c = filter(raw2, grepl("Control_*", raw2$sample))
raw3.c = filter(raw3, grepl("Control_*", raw3$sample))
raw4.c = filter(raw4, grepl("Control_*", raw4$sample))

raw1 = clean(raw1,1,1)
raw2 = clean(raw2,2,81)
raw3 = clean(raw3,3,161)
raw4 = clean(raw4,4,241)


# Collate data
samp = rbind.fill(raw1,raw2,raw3, raw4, raw5)
