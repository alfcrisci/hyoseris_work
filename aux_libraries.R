##########################################################################
# load libraries

#######################################################################################################
# To perform packages' installation before to run code please use pkg::pkg_install function available in R packages install.packages("pkg")


library(tidyverse)
library(ggplot2)
library(ordr)
library(ordr.extra)
library(openxlsx)
library(multiColl)
library(matrixTests)
library(psych)              # psychological research: descr. stats, FA, PCA etc.
library(ggpubr)             # publication ready data visualization in R
library(papeR)
library(rtables) 
library(data.table) 
library(postHoc) 
library(PMCMRplus)
library(ggord)
library(caret)
library(Rmisc)
library(readr)
library(FactoMineR)
library(folda)
library(klaR)
library(knitr)
library(kableExtra)
########################################################################

# install.packages(c("knitr","kableExtra",papeR","rtables","postHoc","PMCMRplus","readr","psych","janitor")

# local installation  postHoc R packages  postHoc1.4.tar.gz  available in root


############################################################################################
# functions

standardize <- function(x) {(x - mean(x))/sd(x)}

stat_groups=function(x,g,xlabv="",ylabv="") {
  myls=list();
  myls[[1]]=aov(x ~ g)
  myls[[2]]=shapiro.test(myls[[1]]$residuals)
  # not normal if 0.05 >
  myls[[3]]=bartlett.test(x ~ g) # not ok if 0.05 >
  # not homogeneous if 0.05 >
  myls[[4]]=kruskal.test(x ~ g)
  # differs 0.05 >
  myls[[5]]=boxplot(x ~ g, xlab=xlabv, ylab=ylabv)
  return(myls)
}

outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

#######################################################################################################
# useful functions

cal_z_score <- function(x){(x - mean(x)) / sd(x)}



############################################################################################
# references


# [1] https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

