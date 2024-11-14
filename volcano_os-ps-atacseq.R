library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/ATAC-seq/Encode pipeline/RUVSeq norm/os-ps")
res <- read.table("results_RUVr_deseq2_k3_default_os-ps.annotated.csv", sep =" ", header = T, stringsAsFactors = F)
res <- res %>%
  mutate(sig = factor(case_when(log2FoldChange <= -log2(1.5) & padj <= 0.05 ~ "Significantly down",
                                log2FoldChange >= log2(1.5) & padj <= 0.05 ~ "Significantly up",
                                TRUE ~ "Not significant")))
table(res$sig)
top = 10
top_genes <- bind_rows(
  res %>% 
    filter(sig == 'Significantly down') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  res %>% 
    filter(sig == 'Significantly up') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top)
)
top_genes$SYMBOL[11] <- "THBS1"
##atheroprone genes
gene_names <- c("EDN1", "ID3", "CCN2", "SERPINE1", "KLF4", "KLF2", "LDLR", "TEK", "CYP1B1", "ID1", "ID2")
atheroprone_genes <- res %>% filter(SYMBOL %in% gene_names)
atheroprone_genes <- atheroprone_genes[atheroprone_genes$sig != "Not significant",]
filtered_data <- atheroprone_genes %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>%
  ungroup()

p1 <- ggplot(res, aes(x = log2FoldChange, y =-log10(padj),color=sig, labels=SYMBOL))+
  ggtitle(label = "", subtitle = "ATAC-Seq: OS vs. PS") +
  geom_point(alpha=0.5) +
  theme_bw(base_size = 18) + 
  #xlim(-10, 10)+
  theme(legend.position = "right") +
  xlab(expression(log[2]("FoldChange"))) + 
  ylab(expression(-log[10]("adjusted pvalue"))) +
  theme(text =element_text(size=18))+
  scale_color_manual(labels = c("Significantly down", "Significantly up", "Not significant"), values=c("grey", "#3399ff","#ff3333"))+
  geom_text_repel(data=filtered_data, aes(label=SYMBOL), size =6,max.overlaps =25, color = "black",  min.segment.length = unit(0, 'lines'))+
  NULL

p1


table(res$sig)

