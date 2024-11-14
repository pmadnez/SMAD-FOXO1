library(dplyr)
library(stringr)
library(tidyr)

###first, get only differentially accesible regions that are up in OS
setwd("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/ATAC-seq/Encode pipeline/RUVSeq norm/os-ps/")
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
setwd("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/ATAC-seq/Encode pipeline/NarrowPeaks/")
motif_instances <- read.csv("SMAD_FOX_GATA_on_differential_OS.txt", sep = "\t", header = T)

smad_foxo1_foxm1 <- motif_instances[,c(2,3,4,16,22,23,24,25,26,27,28,29,35)]

GC_LA_SBE_and_Fox <- smad_foxo1_foxm1[smad_foxo1_foxm1$GC_LA_SBE.Distance.From.Peak.sequence.strand.conservation. != "" & smad_foxo1_foxm1$Foxo1.Forkhead..RAW.Foxo1.ChIP.Seq.Fan.et.al...Homer.Distance.From.Peak.sequence.strand.conservation. != "",]
GC_SBE_and_Fox <- smad_foxo1_foxm1[smad_foxo1_foxm1$GC_SBE.Distance.From.Peak.sequence.strand.conservation. != "" & smad_foxo1_foxm1$Foxo1.Forkhead..RAW.Foxo1.ChIP.Seq.Fan.et.al...Homer.Distance.From.Peak.sequence.strand.conservation. != "",]
SBE_and_Fox <- smad_foxo1_foxm1[smad_foxo1_foxm1$SBE.Distance.From.Peak.sequence.strand.conservation. != "" & smad_foxo1_foxm1$Foxo1.Forkhead..RAW.Foxo1.ChIP.Seq.Fan.et.al...Homer.Distance.From.Peak.sequence.strand.conservation. != "",]
SBE_only <- smad_foxo1_foxm1[smad_foxo1_foxm1$SBE.Distance.From.Peak.sequence.strand.conservation. != "",]
foxo1_only <- smad_foxo1_foxm1[smad_foxo1_foxm1$Foxo1.Forkhead..RAW.Foxo1.ChIP.Seq.Fan.et.al...Homer.Distance.From.Peak.sequence.strand.conservation. != "",]
GC_SBE_only <- smad_foxo1_foxm1[smad_foxo1_foxm1$GC_SBE.Distance.From.Peak.sequence.strand.conservation. != "",]
fox_SMAD_double <- merge(SBE_only, foxo1_only, by.x = 0, by.y = 0)

setwd("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/RNAseq/salmon_pm")
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
genes <- c("EDN1", "ID3")
genes_to_pot <- final_table %>% dplyr::filter(Gene.Name %in% genes)
p1 <- ggplot(final_table, aes(x = log2FoldChange.x, y =log2FoldChange.y,color="red", labels=Gene.Name))+
  ggtitle(label = "", subtitle = "SBE and FOXO1 double positive regions") +
  geom_point(size = 4) +
  theme_bw(base_size = 18) + 
  #xlim(-10, 10)+
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
#GC_SBE_only <- paste(GC_SBE_only$Chr, GC_SBE_only$Start, GC_SBE_only$End, sep = ".")
foxo1_only_regs <- paste(foxo1_only$Chr, foxo1_only$Start, foxo1_only$End, sep = ".")
vector_list_venn <- list("SBE pos" = SBE_only_regs, "FOXO1 pos" = foxo1_only_regs)
venn_overlaps <- calculate.overlap(vector_list_venn)

#vector_list_venn_2 <- list("GC SBE pos" = GC_SBE_only, "FOXO1 pos" = foxo1_only_regs)
#venn_overlaps_2 <- calculate.overlap(vector_list_venn_2)


dev.off()
draw.pairwise.venn(area1 = 669, area2 = 662, cross.area = 486, category = c("FOXO1 positive","SBE positive"), fill = c("firebrick1", "dodgerblue2"), cat.pos = c(-20,0))
dev.off()
draw.triple.venn(area1 = 950, area2 = 669, area3 = 662, n12 = 669, n23 = 486, n13 = 662, n123 = 486, category = c("Diff_peaks_OS","FOXO1 positive","SBE positive"))



##for supp figure, also analyze FOXM1
foxm1_only <- smad_foxo1_foxm1[smad_foxo1_foxm1$FOXM1.Distance.From.Peak.sequence.strand.conservation. != "",]
fox_SMAD_double <- merge(SBE_only, foxm1_only, by.x = 0, by.y = 0)
diff_reg_and_double_pos <- merge(fox_SMAD_double, sig_genes, by.x = 5, by.y = 9)

dev.off()
draw.triple.venn(area1 = 950, area2 = 29, area3 = 662, n12 = 29, n23 = 21, n13 = 662, n123 = 21, category = c("Diff_peaks_OS","FOXM1 positive","SBE positive"), scaled = T)
dev.off()