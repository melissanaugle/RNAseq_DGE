
setwd(dir = "~/Desktop/GitHub/RNAseq_DGE/")
#setwd(dir = "~/Desktop/CSUMB/Thesis/Data analysis/Bioinformatics/DGE_RNAseq_Aug2020/")
rm( list = ls())
graphics.off()
library(tidyverse) 
library(reshape2) 
library(ggthemes) 
library(VennDiagram)


#Read in data files with DEGs between control and each  treatment, selecting only the contig name (ContigID) column
#only needed ones

CoPtE_CoPtH <-read.csv("data/edgeRoutput/subset_files_p.05_all/RSEM.isoform.counts.matrix.CoPt_E_vs_CoPt_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
CoPtE_CoPtH <- setNames(cbind(rownames(CoPtE_CoPtH), CoPtE_CoPtH, row.names = NULL), c("ContigID", colnames(CoPtE_CoPtH)))
head(CoPtE_CoPtH)
CoPtE_CoPtH <- select(CoPtE_CoPtH, "ContigID")

TeleE_TeleH <-read.csv("data/edgeRoutput/subset_files_p.05_all/RSEM.isoform.counts.matrix.Ftele_E_vs_Ftele_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
TeleE_TeleH <- setNames(cbind(rownames(TeleE_TeleH), TeleE_TeleH, row.names = NULL), c("ContigID", colnames(TeleE_TeleH)))
head(TeleE_TeleH)
TeleE_TeleH <- select(TeleE_TeleH, "ContigID")

AluE_AluH <-read.csv("data/edgeRoutput/subset_files_p.05_all/RSEM.isoform.counts.matrix.Falu_E_vs_Falu_H.edgeR.DE_results.P0.05_C2.DE.subset", sep = "\t")
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

#Plot Venn diagram for genes, specifying number of genes in each comparison 
#(area1 = CoPtEvH, area2 = AluEvH, area 3 = TeleEvH, n12 = CoPtAlu, n13 = CtTele, n23 = AluTele, n123 = joinall)


grid.newpage()
venn.plot <- draw.triple.venn( 
  area1 = 3166, area2 = 2820, area3 = 1665, n12 = 1786, n13 = 1018, n23 = 1042, n123 = 830, 
  category = c("Coconut Point", "Faga'alu", "Faga'tele"), fill = c("skyblue", "mediumorchid", "salmon"), 
  lty = "blank", cex = 3);
venn.plot <- draw.triple.venn( area1 = 3166, area2 = 2820, area3 = 1665, n12 = 1786, n13 = 1018, n23 = 1042, n123 = 830, category = c("Coconut Point", "Faga'alu", "Faga'tele"), fill = c("skyblue", "mediumorchid", "salmon"), lty = "blank");
1192+ 956+822+830+188+212+435

#Save the  Venn diagram


png("Venn_allsites.png")
grid.draw(venn.plot)
dev.off()

png("figs/Venn_allsites.png")
grid.draw(venn.plot)
dev.off()


#same with p = 0.001


CoPtE_CoPtH <-read.csv("data/edgeRoutput/subset_files_p.001_all/RSEM.isoform.counts.matrix.CoPt_E_vs_CoPt_H.edgeR.DE_results.P0.001_C2.DE.subset", sep = "\t")
CoPtE_CoPtH <- setNames(cbind(rownames(CoPtE_CoPtH), CoPtE_CoPtH, row.names = NULL), c("ContigID", colnames(CoPtE_CoPtH)))
head(CoPtE_CoPtH)
CoPtE_CoPtH <- select(CoPtE_CoPtH, "ContigID")

TeleE_TeleH <-read.csv("data/edgeRoutput/subset_files_p.001_all/RSEM.isoform.counts.matrix.Ftele_E_vs_Ftele_H.edgeR.DE_results.P0.001_C2.DE.subset", sep = "\t")
TeleE_TeleH <- setNames(cbind(rownames(TeleE_TeleH), TeleE_TeleH, row.names = NULL), c("ContigID", colnames(TeleE_TeleH)))
head(TeleE_TeleH)
TeleE_TeleH <- select(TeleE_TeleH, "ContigID")

AluE_AluH <-read.csv("data/edgeRoutput/subset_files_p.001_all/RSEM.isoform.counts.matrix.Falu_E_vs_Falu_H.edgeR.DE_results.P0.001_C2.DE.subset", sep = "\t")
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

#Plot Venn diagram for genes, specifying number of genes in each comparison 
#(area1 = CoPtEvH, area2 = AluEvH, area 3 = TeleEvH, n12 = CoPtAlu, n13 = CtTele, n23 = AluTele, n123 = joinall)

grid.newpage()
venn.plot <- draw.triple.venn( 
  area1 = 1869, area2 = 1636, area3 = 696, n12 = 1019, n13 = 471, n23 = 509, n123 = 419, 
  category = c("Coconut Point", "Faga'alu", "Faga'tele"), fill = c("skyblue", "mediumorchid", "salmon"), 
  lty = "blank", cex = 3);

#Save the  Venn diagram


png("Venn_allsites_p0.001.png")
grid.draw(venn.plot)
dev.off()




