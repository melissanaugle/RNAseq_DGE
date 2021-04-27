
setwd("/Users/Melissa/Desktop/GitHub/RNAseq_DGE/")

library(tidyverse) 
library(reshape2)  
library(ggplot2)   
library(ggdendro)
library(gplots)
library(RColorBrewer)


# using glm log cpm data
data <- read.csv("data/glmlogcpm_FDR0.5_allheatvcontrol.csv", header=TRUE)
head(data)
nrow(data) #7021

#Check for row duplicates
data[duplicated(data) | duplicated(data, fromLast=TRUE), ]
duplicated(data)

#Re-organize data 


#cluster rows
row.order <- hclust(dist(data))$order 
#cluster columns by sample order 
col.order <-c(1,3,5,7,9,11,12,14,17,19,21,2,4,6,8,10,13,15,16,18,20,22)
#re-order matrix according to clustering
dat_clust <- data[row.order, col.order] # re-order matrix according to clustering
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
  scale_fill_gradient2(low="blue", mid="black", high="orange", midpoint=3.9,    limits=c(-3,19)) + 
  ylab("Genes") +
  xlab("") +
  ggtitle("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=10),
        axis.title = element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 20, size = 0))+
  labs(fill = "Log CPM")

g

ggsave("Figures/glmlogcpm_FDR0.5_allheatvcontrol.png",plot=g)


data_coco <- read.csv("data/glmlogcpm_FDR0.5_COCOheatvcontrol.csv", header=TRUE)
data_alu <- read.csv("data/glmlogcpm_FDR0.5_ALUheatvcontrol.csv", header=TRUE)
data_tele <- read.csv("data/glmlogcpm_FDR0.5_TELEheatvcontrol.csv", header=TRUE)
head(data_tele)
data_tele <- data_tele[,c(1:8)]

row.order <- hclust(dist(data_tele))$order 
col.order <-c(1,3,5,7,2,4,6,8)
dat_clust <- data_tele[row.order, col.order]
df_molten_dat <- melt(as.matrix(dat_clust)) 
names(df_molten_dat)[c(1:2)] <- c("Trinity_ID", "treatment")
df_molten_dat

g<-ggplot(df_molten_dat, aes(x=treatment,y=Trinity_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue", mid="black", high="orange", midpoint=3.9,    limits=c(-3,19)) + 
  ylab("Genes") +
  xlab("") +
  ggtitle("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=10),
        axis.title = element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 20, size = 0))+
  labs(fill = "Log CPM")

g

#ggsave("Figures/glmlogcpm_FDR0.5_TELEheatvcontrol.png",plot=g)
