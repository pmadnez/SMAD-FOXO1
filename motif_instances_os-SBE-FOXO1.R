library(dplyr)
library(stringr)
library(tidyr)

###first, get only differentially accesible regions that are up in OS
setwd("Z:/data/")
diff_regs <- read.csv("desults_RUVr_deseq2_k3_default_os-ps.csv", sep = ",")
diff_regs <- diff_regs %>%
  mutate(sig = factor(case_when(log2FoldChange <= -0.585 & padj <= 0.05 ~ "Significantly down",
                                log2FoldChange >= 0.585 & padj <= 0.05 ~ "Significantly up",
                                TRUE ~ "Not significant")))
sig_regs_up <- diff_regs[diff_regs$sig == "Significantly up",]
sig_regs_up.bed <- as.data.frame(sig_regs_up[1])
sig_regs_up.bed <- sig_regs_up.bed %>% separate(X, c("chrom", "start", "end"))
write.table(sig_regs_up.bed, "OS_diff-regs.bed", quote = F, row.names = F, sep = "\t")

###perform motif instance calling with HOMER on linux server
# annotatePeaks.pl /project/Mechanogenomics_data/ATAC_ENCODE_Analysis_pm/Homer/os-ps/find_instances_of_motifs/Bed/OS_diff-regs.bed hg38 -m /home/mendez/bin/homer/motifs/SMAD_FOX_GATA_lib.motif > /project/Mechanogenomics_data/ATAC_ENCODE_Analysis_pm/Homer/os-ps/find_instances_of_motifs/output/SMAD_FOX_GATA_on_differential_OS.txt
#now identify regions with shared motifs
setwd("Z:/data/")
motif_instances <- read.csv("SMAD_FOX_GATA_on_differential_OS.txt", sep = "\t", header = T)

smad_foxo1_foxm1 <- motif_instances[,c(2,3,4,16,22,23,24,25,26,27,28,29,35)]

SBE_only <- smad_foxo1_foxm1[smad_foxo1_foxm1$SBE.Distance.From.Peak.sequence.strand.conservation. != "",]
foxo1_only <- smad_foxo1_foxm1[smad_foxo1_foxm1$Foxo1.Forkhead..RAW.Foxo1.ChIP.Seq.Fan.et.al...Homer.Distance.From.Peak.sequence.strand.conservation. != "",]
fox_SMAD_double <- merge(SBE_only, foxo1_only, by.x = 0, by.y = 0)

##read RNAseq data to only find motifs in DEGs
res <- read.table("deseq2_os-vs-ps.salmon.hgnc-symbol.csv", sep = ";", header = T, stringsAsFactors = F)
res <- res %>%
  mutate(sig = factor(case_when(log2FoldChange <= -0.585 & padj <= 0.05 ~ "Significantly down",
                                log2FoldChange >= 0.585 & padj <= 0.05 ~ "Significantly up",
                                TRUE ~ "Not significant")))
sig_genes <- res[res$sig == "Significantly up",]

diff_reg_and_double_pos <- merge(SBE_and_Fox, sig_genes, by.x = 4, by.y = 9)

diff_reg_and_double_pos$Start <- diff_reg_and_double_pos$Start-1 #homer adds a +1 to start
diff_reg_and_double_pos$X <- paste(diff_reg_and_double_pos$Chr, diff_reg_and_double_pos$Start, diff_reg_and_double_pos$End, sep = ".")
final_table <- merge(diff_reg_and_double_pos, sig_regs_up, by = "X")

##make the plot
library(ggplot2)
top = 1
top_genes <- bind_rows(
  final_table %>% 
    arrange(desc(log2FoldChange.y)) %>% 
    head(top))
genes <- c("EDN1", "ID3") # get BMP target genes identified before
genes_to_pot <- final_table %>% dplyr::filter(Gene.Name %in% genes)
p1 <- ggplot(final_table, aes(x = log2FoldChange.x, y =log2FoldChange.y,color="red", labels=Gene.Name))+
  ggtitle(label = "", subtitle = "SBE and FOXO1 double positive regions") +
  geom_point(size = 4) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "right") +
  xlab(expression(log[2]("FoldChange RNA-Seq"))) + 
  ylab(expression(log[2]("FoldChange ATAC-Seq"))) +
  theme(text =element_text(size=18))+
  theme(aspect.ratio = 1)+
  geom_text_repel(data=genes_to_pot, aes(label=Gene.Name), size =6,max.overlaps =25, color = "black", segment.color = "black", hjust= "outward", box.padding = 0.2, point.padding = 0.35, show.legend = FALSE, min.segment.length = 35)+
  NULL
p1

###make Venn diagram to visualize overlap
library(VennDiagram)
#make vectors of regions for each file by concatenating Chr, Start, End columns
overall_regs <- paste(smad_foxo1_foxm1$Chr, smad_foxo1_foxm1$Start, smad_foxo1_foxm1$End, sep = ".")
SBE_only_regs <- paste(SBE_only$Chr, SBE_only$Start, SBE_only$End, sep = ".")
foxo1_only_regs <- paste(foxo1_only$Chr, foxo1_only$Start, foxo1_only$End, sep = ".")
vector_list_venn <- list("SBE pos" = SBE_only_regs, "FOXO1 pos" = foxo1_only_regs)
venn_overlaps <- calculate.overlap(vector_list_venn)

dev.off()
draw.pairwise.venn(area1 = 669, area2 = 662, cross.area = 486, category = c("FOXO1 positive","SBE positive"), fill = c("firebrick1", "dodgerblue2"), cat.pos = c(-20,0))
dev.off()
draw.triple.venn(area1 = 950, area2 = 669, area3 = 662, n12 = 669, n23 = 486, n13 = 662, n123 = 486, category = c("Diff_peaks_OS","FOXO1 positive","SBE positive"))

dev.off()
