
setwd(dir = "~/Desktop/GitHub/corals_RNAseq_DGE/")
#setwd(dir = "~/Desktop/CSUMB/Thesis/Data analysis/Bioinformatics/DGE_RNAseq_Aug2020/")
rm( list = ls())
graphics.off()
library(tidyverse) 
library(reshape2) 
library(ggthemes) 
library(VennDiagram)


#Read in data files with DEGs between control and each hypoxia treatment, selecting only the contig name (ContigID) column


#only needed ones

CoPtE_CoPtH <-read.csv("data/RSEM.isoform.counts.matrix.CoPt_E_vs_CoPt_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")

'''
CoPtE_CoPtH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_CoPt_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_CoPtH <- setNames(cbind(rownames(CoPtE_CoPtH), CoPtE_CoPtH, row.names = NULL), c("ContigID", colnames(CoPtE_CoPtH)))
head(CoPtE_CoPtH)
CoPtE_CoPtH <- select(CoPtE_CoPtH, "ContigID")

CoPtE_AluE <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_Falu_E.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_AluE <- setNames(cbind(rownames(CoPtE_AluE), CoPtE_AluE, row.names = NULL), c("ContigID", colnames(CoPtE_AluE)))
head(CoPtE_AluE)

CoPtE_AluH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_Falu_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_AluH <- setNames(cbind(rownames(CoPtE_AluH), CoPtE_AluH, row.names = NULL), c("ContigID", colnames(CoPtE_AluH)))
head(CoPtE_AluH)

CoPtE_TeleE <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_Ftele_E.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_TeleE <- setNames(cbind(rownames(CoPtE_TeleE), CoPtE_TeleE, row.names = NULL), c("ContigID", colnames(CoPtE_TeleE)))
head(CoPtE_TeleE)

CoPtE_TeleH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_TeleH <- setNames(cbind(rownames(CoPtE_TeleH), CoPtE_TeleH, row.names = NULL), c("ContigID", colnames(CoPtE_TeleH)))
head(CoPtE_TeleH)

CoPtH_AluE <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_H_vs_Falu_E.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtH_AluE <- setNames(cbind(rownames(CoPtH_AluE), CoPtH_AluE, row.names = NULL), c("ContigID", colnames(CoPtH_AluE)))
head(CoPtH_AluE)

CoPtH_AluH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_H_vs_Falu_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtH_AluH <- setNames(cbind(rownames(CoPtH_AluH), CoPtH_AluH, row.names = NULL), c("ContigID", colnames(CoPtH_AluH)))
head(CoPtH_AluH)

CoPtH_TeleE <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_H_vs_Ftele_E.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtH_TeleE <- setNames(cbind(rownames(CoPtH_TeleE), CoPtH_TeleE, row.names = NULL), c("ContigID", colnames(CoPtH_TeleE)))
head(CoPtH_TeleE)

CoPtH_TeleH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_H_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtH_TeleH <- setNames(cbind(rownames(CoPtH_TeleH), CoPtH_TeleH, row.names = NULL), c("ContigID", colnames(CoPtH_TeleH)))
head(CoPtH_TeleH)

AluE_TeleE <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.Falu_E_vs_Ftele_E.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
AluE_TeleE <- setNames(cbind(rownames(AluE_TeleE), AluE_TeleE, row.names = NULL), c("ContigID", colnames(AluE_TeleE)))
head(AluE_TeleE)

AluE_TeleH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.Falu_E_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
AluE_TeleH <- setNames(cbind(rownames(AluE_TeleH), AluE_TeleH, row.names = NULL), c("ContigID", colnames(AluE_TeleH)))
head(AluE_TeleH)

AluH_TeleE <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.Falu_H_vs_Ftele_E.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
AluH_TeleE <- setNames(cbind(rownames(AluH_TeleE), AluH_TeleE, row.names = NULL), c("ContigID", colnames(AluH_TeleE)))
head(AluH_TeleE)

AluH_TeleH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.Falu_H_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
AluH_TeleH <- setNames(cbind(rownames(AluH_TeleH), AluH_TeleH, row.names = NULL), c("ContigID", colnames(AluH_TeleH)))
head(AluH_TeleH)
'''

#only needed ones

CoPtE_CoPtH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_CoPt_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_CoPtH <- setNames(cbind(rownames(CoPtE_CoPtH), CoPtE_CoPtH, row.names = NULL), c("ContigID", colnames(CoPtE_CoPtH)))
head(CoPtE_CoPtH)
CoPtE_CoPtH <- select(CoPtE_CoPtH, "ContigID")

TeleE_TeleH <-read.csv("data/RSEM.isoform.counts.matrix.Ftele_E_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
TeleE_TeleH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.Ftele_E_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
TeleE_TeleH <- setNames(cbind(rownames(TeleE_TeleH), TeleE_TeleH, row.names = NULL), c("ContigID", colnames(TeleE_TeleH)))
head(TeleE_TeleH)
TeleE_TeleH <- select(TeleE_TeleH, "ContigID")

AluE_AluH <-read.csv("data/RSEM.isoform.counts.matrix.Falu_E_vs_Falu_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
AluE_AluH <-read.csv("edgeR.22083.dir_DGE_ALL/subset_files_p.05_all/RSEM.isoform.counts.matrix.Falu_E_vs_Falu_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
AluE_AluH <- setNames(cbind(rownames(AluE_AluH), AluE_AluH, row.names = NULL), c("ContigID", colnames(AluE_AluH)))
head(AluE_AluH)
AluE_AluH <- select(AluE_AluH, "ContigID")



#Use join function to find contig names shared between comparisons. Head(joinX) will print your list of shared genes.

joinCoPtEHvAluEH <- intersect(CoPtE_CoPtH, AluE_AluH, by="ContigID")
joinCoPtEHvTeleEH <- intersect(CoPtE_CoPtH,TeleE_TeleH, by="ContigID")
joinAluEHvTeleEH <- intersect(AluE_AluH,TeleE_TeleH, by="ContigID")
join_all_2 <- intersect(joinCoPtEHvAluEH,joinCoPtEHvTeleEH,joinAluEHvTeleEH)

#head(joinCoPtEHvTeleEH)
#Plot Venn diagram with numbers generated above

#Plot Venn diagram for upregulated genes, specifying number of genes in each comparison (area1 = Cv2, area2 = Cv4, area 3 = Cv6)


grid.newpage()
venn.plot <- draw.triple.venn( 
  area1 = 3166, area2 = 2820, area3 = 1665, n12 = 1786, n13 = 1018, n23 = 1042, n123 = 830, 
  category = c("Coconut Point", "Faga'alu", "Faga'tele"), fill = c("skyblue", "mediumorchid", "salmon"), 
  lty = "blank", cex = 3);
venn.plot <- draw.triple.venn( area1 = 3166, area2 = 2820, area3 = 1665, n12 = 1786, n13 = 1018, n23 = 1042, n123 = 830, category = c("Coconut Point", "Faga'alu", "Faga'tele"), fill = c("skyblue", "mediumorchid", "salmon"), lty = "blank");


#Save the  Venn diagram


png("Venn_allsites.png")
grid.draw(venn.plot)
dev.off()

png("figs/Venn_allsites.png")
grid.draw(venn.plot)
dev.off()




