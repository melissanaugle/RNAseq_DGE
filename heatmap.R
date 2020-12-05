
setwd("/Users/Melissa/Desktop/GitHub/corals_RNAseq_DGE/")

library(tidyverse) 
library(reshape2)  
library(ggplot2)   
library(ggdendro)
library(gplots)
library(RColorBrewer)


#now using p = 0.05 data 
data <- read.csv("data/diffExpr.P0.05_C2.matrix.log2.centered.dat.csv", header=TRUE)
head(data)

#Check for row duplicates

data[duplicated(data) | duplicated(data, fromLast=TRUE), ]
duplicated(data)


#Re-organize data 

#include only numerical columns (expression values)
dat <- data[,2:length(data)] 
#row names are 1st column (Trinity_ID)
rownames(dat) <- data[,1]
head(dat)
nrow(dat)
#cluster rows
row.order <- hclust(dist(dat))$order 
#cluster columns by sample order 
col.order <-c(1,2,3,4,9,10,11,12,16,17,18,5,6,7,8,13,14,15,19,20,21,22) 
#re-order matrix according to clustering
dat_clust <- dat[row.order, col.order] # re-order matrix according to clustering
#reshape into data frame
df_molten_dat <- melt(as.matrix(dat_clust)) 
names(df_molten_dat)[c(1:2)] <- c("Trinity_ID", "treatment")
df_molten_dat

#Find min and max of value column in data frame to determine range of heatmap

max(df_molten_dat$value)
min(df_molten_dat$value)
mean(df_molten_dat$value)
IQR(df_molten_dat$value)


#plot heatmap
g<-ggplot(df_molten_dat, aes(x=treatment,y=Trinity_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue", mid="black", high="orange", midpoint=0,    limits=c(-5,5)) + 
  ylab("Genes") +
  xlab("Treatment") +
  ggtitle("DGE in corals under heat stress at 3 sites") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=10),
        axis.title = element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 20, size = 5))+
  labs(fill = "Log2 Fold Change")

g

#ggsave("Heatmap1.png",plot=g)

