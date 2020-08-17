library(openxlsx); library(data.table);library(tidyverse)
library(reshape2); #library(ClassDiscovery); library(phytools)
library(dendextend); library(circlize); library(lemon)
library(ggpubr); library(ggthemes); library(AGHmatrix)

rm(list = ls())
path <- "C:/Users/deer1/Google Drive/Alejandro/02. IBDLD/"
path1 = "C:/Users/deer1/Google Drive/Alejandro/02. IBDLD/03. Resultados/"

GIBDLD <- as.data.frame(fread(paste0(path,"GIBDLD_genome.kinship")))
colnames(GIBDLD) <- c("fid1","iid1","fid2","iid2","nSNPs",
                      paste0("D",1:9),"D0")
  
GIBDLD <- GIBDLD[,c("iid1","iid2","D0")]
a <-  GIBDLD[,c("iid2","iid1","D0")]
colnames(a) <- colnames(GIBDLD)
a$flag <- ifelse(a$iid1 == a$iid2,1,0)
a <- subset(a, a$flag == 0)
a <- a[,-4]
GIBDLD1 <- rbind(GIBDLD,a)

GIBDLD1$IBD_Diss <- 1 - GIBDLD1$D0

# GIBDLD Dendrogram

IBDLD_Den <- GIBDLD1[, c("iid1","iid2","IBD_Diss")]

IBDLD_Den <- reshape(IBDLD_Den, timevar = "iid1", idvar = "iid2", direction = "wide")
IBDLD_Den <- IBDLD_Den %>% remove_rownames() %>% column_to_rownames("iid2")
colnames(IBDLD_Den) <- gsub("IBD_Diss.","",
                            colnames(IBDLD_Den))

IBDLD_Den <- IBDLD_Den[order(rownames(IBDLD_Den)), 
                       order(colnames(IBDLD_Den))]

IBDLD_Den_dis <- as.dist(as.matrix(IBDLD_Den))
IBDLD_Den_clus <- hclust(IBDLD_Den_dis, method = "average", 
                         members = NULL)

dend <- as.dendrogram(IBDLD_Den_clus)
clucols <- c("darkmagenta",
             'darkgoldenrod3',
             "gray50"
             )
dend <- dend %>% 
  color_branches(k=3, col = clucols, lty = 3) %>% 
  color_labels(k = 3, col = clucols)%>%
  set("branches_lwd",3.2)

circlize_dendrogram(dend)



dend <- as.dendrogram(IBDLD_Den_clus)
