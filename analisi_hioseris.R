################################################################################
# To achieve  reproducibility setup wdir and seeds

# setwd("")

set.seed(123)

##########################################################################################################
# Load R packages

# check if all libraries are installed
source("aux_libraries.R")
library(emmeans)
library(corrplot)
library(MKinfer) # https://github.com/stamats/MKinfer

cat("\014") 

################################################################################
# Loading data 

dati_rese=read.xlsx("truschi_work_mat_27052025.xlsx",1) # dati relativi alle rese
dati_parameters=read.xlsx("truschi_work_mat_27052025.xlsx",2) # dati relativi alle rese
dati_chemicals=read.xlsx("truschi_work_mat_27052025.xlsx",3) # dati relativi alle rese
dati_microgreens=read.xlsx("truschi_work_mat_27052025.xlsx",4) # dati relativi alle rese
dati_baby_leafs=read.xlsx("truschi_work_mat_27052025.xlsx",5) # dati relativi alle rese


dati_hioseris=list(dati_rese=dati_rese,
                   dati_parameters=dati_parameters,
                   dati_chemicals=dati_chemicals,
                   dati_microgreens=dati_microgreens,
                   dati_baby_leafs=dati_baby_leafs
                   )

cat("\014") 

################################################################################
# eliminazione dei dati outlier e salvataggio su un file di archiviazione R formato rds 


saveRDS(dati_hioseris,"dati_hioseris.rds")
dati_hioseris=readRDS("dati_hioseris.rds")


################################################################################
# Visualizzazione delle matrici di lavoro  e dei  nomi delle variabili

lapply(dati_hioseris,names)
lapply(dati_hioseris,nrow)


################################################################################
# Rese hioseris anova 2 way

dati_rese=dati_hioseris$dati_rese
result = bartlett.test(log(harvest) ~ interaction(treatment,stage), data = dati_rese)
lm_rese=lm(log(harvest)~treatment+stage,data=dati_rese)
summary(lm_rese)
rese.s <- emmeans(lm_rese,~ stage+treatment)
pairs(rese.s)
plot(rese.s, comparisons = TRUE)

dati_rese$f1f2 <- interaction(dati_rese$treatment, dati_rese$stage)
ggplot(aes(y = harvest, x = treatment, fill = stage), data = dati_rese) + geom_boxplot()

ggsave("rese_boxplot.png")

cat("\014") 

##########################################################################
# Creazione delle matrici di analisi per PCA

dati_parameters=dati_hioseris$dati_parameters

#########################################################################################
# Parameters: PCA explore parameters and outlier detection

X=dati_parameters[,7:11]
Y=paste(dati_parameters$stage,dati_parameters$treatment)

data_pca=data.frame(X,Treatments=Y)

hioseris_pca=ordinate(data_pca , cols = 1:5, model = ~ prcomp(., scale. = TRUE)) 

hioseris_pca %>%
  augment_ord() %>%
  ggbiplot(axis.type = "predictive") +
  theme_bw() +
  geom_rows_point(alpha = .3,aes(color=Treatments,shape=Treatments))+
  stat_rows_center(alpha = .8, fun.center = "mean",size = 5,aes(color = Treatments))+
  geom_origin() +
  ggtitle("PCA - Hioseris data") +
  labs(color = "Treatments")

ggsave("PCA_biplot_hioseris.png") # by using ord R package

######################################################################################
# by using factoextra R packages

hioseris_pca <- prcomp(X,scale = T) # scaled!

color_cluster_dendro <- c("#FF0000", #cluster 1
                          "#669900", #cluster 2
                          "#00CCCC", #cluster 3
                          "#9933FF") #cluster 4



PCA_final_hioseris <- factoextra::fviz_pca_biplot(hioseris_pca, 
                                                   # fill individuals by groups
                                                   geom.ind = c("point"),
                                                   col.ind = "black",
                                                   fill.ind = Y,
                                                   pointshape = 21, #i numeri definiscono  uno stile di punti!
                                                   palette = "color_cluster_dendro", 
                                                   addEllipses = T, 
                                                   ellipse.level = 0.10,
                                                   ellipse.type = "convex",
                                                   geom.var = c("arrow", "text"), 
                                                   arrowsize = 0.3,
                                                   labelsize = 2,
                                                   legend.title = "Treatment",
                                                   title = "PCA - Hioseris data",
                                                   repel = T  
) +
  ggpubr::color_palette("color_cluster_dendro") +
  theme(legend.position = "bottom")


PCA_final_hioseris

ggsave(filename="PCA_hioseris_novel.png")

cat("\014")

#########################################################################################
# data analisys ANOVA & Kruskal Wallis ( optional) 

# [1] "A_leaf_p"  "fresh_w_p" "dry_w_p"   "ps_p"      "w_tot"  

########################################################################################


#########################################################################################
# Preparazione dati chimici

X=dati_chemicals[,4:34]

X$vitb6=NULL
X

Y=paste(dati_chemicals$stage,dati_chemicals$treatment)

hioseris_pca_chem <- prcomp(X,scale = T) # scaled!

color_cluster_dendro <- c("#FF0000", #cluster 1
                          "#669900", #cluster 2
                          "#00CCCC", #cluster 3
                          "#9933FF") #cluster 4


PCA_final_hioseris_chem <- factoextra::fviz_pca_biplot(hioseris_pca_chem, 
                                                  # fill individuals by groups
                                                  geom.ind = c("point"),
                                                  col.ind = "black",
                                                  fill.ind = Y,
                                                  pointshape = 21, #i numeri definiscono  uno stile di punti!
                                                  palette = "color_cluster_dendro", 
                                                  addEllipses = T, 
                                                  ellipse.level = 0.10,
                                                  ellipse.type = "convex",
                                                  geom.var = c("arrow", "text"), 
                                                  arrowsize = 0.3,
                                                  labelsize = 2,
                                                  legend.title = "Treatment",
                                                  title = "PCA - Hioseris chem data",
                                                  repel = T  
) +
  ggpubr::color_palette("color_cluster_dendro") +
  theme(legend.position = "bottom")


PCA_final_hioseris_chem

ggsave(filename="PCA_hioseris_chem.png")

cat("\014")

###########################################################################################

X=X$vitB6=NULL

names(X)[which(col_oneway_welch(X, Y)$pvalue<0.05)] # quelli che danno diffrenze significative


# [1] "hr"       "phenols"  "prot"     "nitritis" "chl_a"    "chl_b"    "chl_tot" 
# [8] "crom"     "mn"       "copper"   "arsenic"  "cadmium"  "sb"       "vitC"    
# [15] "vitb3"    "vitb1"  

########################################################################################
# confronti paired tra le due tesi stage e terreno per veder qualis sono quelle che generano differenze

ls_res_paired_chem=list()

for ( i in 1:ncol(X)) {
  ls_res_paired_chem[[i]]=pairwise.wilcox.exact(X[,i], Y, p.adjust.method = "bonferroni") # citation("MKinfer")
}

names(ls_res_paired_chem)=names(X)

########################################################################################
# Analisi di correlazione

X_chem<-split(X, Y)

corr_mat_bl_all=cor(X_chem$`Baby-leaf all`,method="s")
png(height=1000, width=1000, file="corr_mat_bl_all.png", type = "cairo")
corrplot(corr_mat_bl_all,title="Baby-leaf all",order = 'FPC', type = 'lower', diag = FALSE)
dev.off()



corr_mat_bl_half=cor(X_chem$`Baby-leaf half`,method="s")
png(height=1000, width=1000, file="corr_mat_bl_half.png", type = "cairo")
corrplot(corr_mat_bl_half,title="Baby-leaf half",order = 'FPC', type = 'lower', diag = FALSE)
dev.off()


corr_mat_micro_all=cor(X_chem$`Microgreens all`,method="s")
png(height=1000, width=1000, file="corr_mat_micro_all.png", type = "cairo")
corrplot(corr_mat_micro_all,title="Microgreens all",order = 'FPC', type = 'lower', diag = FALSE)
dev.off()



corr_mat_micro_half=cor(X_chem$`Microgreens half`,method="s")
png(height=1000, width=1000, file="corr_mat_micro_half.png", type = "cairo")
corrplot(corr_mat_micro_half,title="Microgreens half",order = 'FPC', type = 'lower', diag = FALSE)
dev.off()




#####################################################################
# References

# https://pages.cms.hu-berlin.de/EOL/gcg_quantitative-methods/Lab11_LDA_Model-assessment.html
# https://dicook.github.io/mulgar_book/14-lda.html
# https://stats.stackexchange.com/questions/560375/how-to-interpret-the-output-of-lda-discriminant-analysis-in-r
# https://cmdlinetips.com/2020/12/canonical-correlation-analysis-in-r/
# https://www.r-bloggers.com/2021/05/linear-discriminant-analysis-in-r/
# https://vitalflux.com/pca-vs-lda-differences-plots-examples/
# https://towardsai.net/p/data-science/lda-vs-pca
# https://stats.stackexchange.com/questions/23353/pca-lda-cca-and-pls
# https://www.geeksforgeeks.org/classifying-data-using-support-vector-machinessvms-in-r/
# https://mdatools.com/docs/pca.html
# https://rdrr.io/cran/mixOmics/man/plsda.html
# http://mixomics.org/methods/spls-da/
# https://stackoverflow.com/questions/46977743/pca-analysis-remove-centroid
