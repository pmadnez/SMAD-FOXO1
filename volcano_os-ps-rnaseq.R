library(ggplot2)
library(dplyr)
library(ggrepel)
###Fig2D
setwd("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/RNAseq/salmon_pm")
res <- read.table("deseq2_os-vs-ps.salmon.hgnc-symbol.csv", sep = ";", header = T, stringsAsFactors = F)
res <- res %>%
  mutate(sig = factor(case_when(log2FoldChange <= -0.585 & padj <= 0.05 ~ "Significantly down",
                                    log2FoldChange >= 0.585 & padj <= 0.05 ~ "Significantly up",
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
sox_genes <- bind_rows(
  res %>% 
    filter(grepl('SOX', hgnc_symbol)))

fox_genes <- bind_rows(
  res %>% 
    filter(grepl('FOX', hgnc_symbol)))
fox_genes <- fox_genes[fox_genes$sig != "Not significant",]
fox_genes <- fox_genes[-2,]
gata_genes <- bind_rows(
  res %>% 
    filter(grepl('GATA', hgnc_symbol)))
gata_genes <- gata_genes[gata_genes$sig != "Not significant",]
gata_genes <- gata_genes[-2,]

##atheroprone genes
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
  #geom_text_repel(data= sox_genes, mapping = aes(log2FoldChange, -log(pvalue,10), label = hgnc_symbol), show.legend = F, max.overlaps = 30, box.padding = 0.5, segment.color= "black") +
  scale_color_manual(labels = c("Significantly down", "Significantly up", "Not significant"), values=c("grey", "#3399ff","#ff3333"))+
  geom_text_repel(data=atheroprone_genes, aes(label=hgnc_symbol), size =6,max.overlaps =300, color = "black", segment.width =   min.segment.length = unit(0, 'lines'))+
  NULL

p1
color = ("black")
p1 + 
  geom_point(aes(x=1.2346929, y=-log(0.000000741,10), colour=color, alpha=1))+
  annotate("text",size = 5, x=4.5, y=7.5, label= "FoxM1", colour=color, fontface = "bold")+
  geom_point(aes(x=0.4261073, y=-log(0.045139330,10), colour=color, alpha=1))+
  annotate("text",size = 5, x=3.5, y=1.6, label= "FoxO1", colour=color, fontface = "bold")+
  NULL

p1 + 
  geom_point(aes(x=0.9702429, y=-log(0.000417805,10), colour=color, alpha=1))+
  annotate("text",size = 5, x=2, y=5, label= "GATA3", colour=color, fontface = "bold")+
  geom_point(aes(x=0.4973178, y=-log(0.027994832,10), colour=color, alpha=1))+
  annotate("text",size =5, x=2.6, y=2, label= "GATA2", colour=color, fontface = "bold")+
  NULL
-log(0.000417805,10)
-log(0.030571974,10)
-log(0.027994832,10)
