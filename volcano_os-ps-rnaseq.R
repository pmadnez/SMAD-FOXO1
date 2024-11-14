library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("Z:/data")
res <- read.table("deseq2_os-vs-ps.salmon.hgnc-symbol.csv", sep = ";", header = T, stringsAsFactors = F)
res <- res %>%
  mutate(sig = factor(case_when(log2FoldChange <= -0.585 & padj <= 0.05 ~ "Significantly down",
                                    log2FoldChange >= 0.585 & padj <= 0.05 ~ "Significantly up",
                                    TRUE ~ "Not significant")))

table(res$sig)


##atheroprone/ atheroprotective genes
gene_names <- c("EDN1", "ID3", "CCN2", "SERPINE1", "KLF4", "KLF2", "LDLR", "TEK", "PTGDS", "ID1", "ID2")
atheroprone_genes <- res %>% filter(hgnc_symbol %in% gene_names)
atheroprone_genes <- atheroprone_genes[atheroprone_genes$sig != "Not significant",]


p1 <- ggplot(res, aes(x = log2FoldChange, y =-log10(pvalue),color=sig, labels=hgnc_symbol))+
  ggtitle(label = "", subtitle = "OS vs. PS") +
  geom_point() +
  theme_bw(base_size = 18) + 
  xlim(-10, 10)+
  theme(legend.position = "right") +
  xlab(expression(log[2]("FoldChange"))) + 
  ylab(expression(-log[10]("pvalue"))) +
  theme(text =element_text(size=18))+
  theme(aspect.ratio = 1)+
  scale_color_manual(labels = c("Significantly down", "Significantly up", "Not significant"), values=c("grey", "#3399ff","#ff3333"))+
  geom_text_repel(data=atheroprone_genes, aes(label=hgnc_symbol), size =6,max.overlaps =300, color = "black", segment.width =   min.segment.length = unit(0, 'lines'))+
  NULL

p1
