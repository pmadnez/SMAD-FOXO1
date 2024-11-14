library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("Z:/data/")
###read RNAseq and ATAC-Seq data
os_ps.rnaseq <- read.table("deseq2_os-vs-ps.salmon.hgnc-symbol.pval.fc.csv")
os_ps_atacseq <- read.table("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/ATAC-seq/Encode pipeline/RUVSeq norm/os-ps/results_RUVr_deseq2_k3_default_os-ps.annotated.csv", sep = " ")

#get rid of duplicates to be able to display in plot (multiple peaks per gene is difficult, keep only highest absolute peaks)
atac_wo_dups <- os_ps_atacseq[order(os_ps_atacseq$SYMBOL, -abs(os_ps_atacseq$log2FoldChange) ), ]
atac_wo_dups <- atac_wo_dups[ !duplicated(atac_wo_dups$SYMBOL), ]

##to avoid irritation, also get rid of transcript variants and only keep the highest expressed one (highest FC)
rna_wo_dups <- os_ps.rnaseq[order(os_ps.rnaseq$hgnc_symbol, -abs(os_ps.rnaseq$log2FoldChange) ), ]
rna_wo_dups <- rna_wo_dups[ !duplicated(rna_wo_dups$hgnc_symbol), ]

atac_sig_fc <- atac_wo_dups[(abs(atac_wo_dups$log2FoldChange)) > log(1.5, 2) & atac_wo_dups$padj < 0.05,]

atac_sig_fc <- atac_sig_fc[complete.cases(atac_sig_fc),]

atacseq_for_plot <- atac_sig_fc[,c(16,24)]
rnaseq_for_plot <- rna_wo_dups[,c(3,8)]

merged_for_plot <- merge(rnaseq_for_plot, atacseq_for_plot, by.x = "hgnc_symbol", by.y = "SYMBOL")

merged_for_plot <- merged_for_plot %>%
  mutate(regulation = factor(case_when(log2FoldChange.x < 0 & log2FoldChange.y <0 ~ "Concomitantly down",
                                       log2FoldChange.x > 0 & log2FoldChange.y >0 ~ "Concomitantly up",
                                       TRUE ~ "Conversely regulated")))

cor(merged_for_plot$log2FoldChange.x, merged_for_plot$log2FoldChange.y, method = "pearson")
table(merged_for_plot$regulation)

gene_names <- c("EDN1", "ID3", "CCN2", "SERPINE1", "KLF4", "KLF2", "TEK", "CYP1B1") #define gene names to plot
genes_to_plot <- merged_for_plot %>% filter(hgnc_symbol %in% gene_names)

ggplot(merged_for_plot, aes(x = log2FoldChange.x, y = log2FoldChange.y, color = regulation)) +
  ggtitle(label = "Integrative Analysis RNAseq and ATACseq", subtitle = "OS vs PS, FC > 1.5, padj < 0.05") +
  geom_point() +
  theme(legend.position = "right")+
  xlab(expression(log[2]("FoldChange RNAseq"))) +
  xlim(-7.5,7.5)+
  ylab(expression(log[2]("FoldChange ATACseq"))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_minimal(base_size = 18)+
  geom_text_repel(data= genes_to_plot, mapping = aes(log2FoldChange.x, log2FoldChange.y, label = hgnc_symbol),colour = "black",size = 6, show.legend = F, max.overlaps = 30, box.padding = 0.5, segment.color= "black") +
  scale_color_manual(values = c("Concomitantly down" = "#4169e1", "Concomitantly up" = "#ff6347", "Conversely regulated" = "black"))+
  NULL

conc_up <- merged_for_plot[merged_for_plot$regulation == "Concomitantly up",]
conc_down <- merged_for_plot[merged_for_plot$regulation == "Concomitantly down",]
convers_reg <- merged_for_plot[merged_for_plot$regulation == "Conversely regulated",]

conc_up_complete <- merge(conc_up, atac_sig_fc, by.x = "hgnc_symbol", by.y = "SYMBOL")
conc_down_complete <- merge(conc_down, atac_sig_fc, by.x = "hgnc_symbol", by.y = "SYMBOL")
convers_reg_complete <- merge(convers_reg, atac_sig_fc, by.x = "hgnc_symbol", by.y = "SYMBOL")

conc_down_bed <- conc_down_complete[,c(5:7)]
names(conc_down_bed) <- c("chr", "start", "end")

conc_up_bed <- conc_up_complete[,c(5:7)]
names(conc_up_bed) <- c("chr", "start", "end")

convers_reg_bed <- convers_reg_complete[,c(5:7)]
names(convers_reg_bed) <- c("chr", "start", "end")


write.table(conc_down_bed, "conc_down_os-vs-ps.bed", sep = "\t", row.names = F, quote = F)
write.table(conc_up_bed, "conc_up_os-vs-ps.bed", sep = "\t", row.names = F, quote = F)
write.table(conc_down_complete, "conc_down_complete_os-vs-ps.csv", sep = ",", row.names = F, quote = F)
write.table(conc_up_complete, "conc_up_complete_os-vs-ps.csv", sep = ",", row.names = F, quote = F)
write.table(conc_down_complete, "conc_down_complete_os-vs-ps.txt", sep = "\t", row.names = F, quote = F)
write.table(conc_up_complete, "conc_up_complete_os-vs-ps.txt", sep = "\t", row.names = F, quote = F)

write.table(convers_reg_bed, "convers-reg-os-ps.bed", sep = "\t", row.names = F, quote = F)
write.table(convers_reg_complete, "convers-reg_os-vs-ps-complete.txt", sep = "\t", row.names = F, quote = F)

##continue with motif enrichment etc.
