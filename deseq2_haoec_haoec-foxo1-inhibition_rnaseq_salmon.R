#BiocManager::install("tximport")
#install.packages("readr")
#BiocManager::install("ensembldb")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("tximeta")
#library(tximeta)
library(Gen)
library(tximport)
library(DESeq2)
library(readr)
library(edgeR)
#library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library('biomaRt')
library(ggpp)
library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("Z:/Shared Workspaces/Paul2Jerome/Mechanogenomics/Endothel Zellen/RNAseq/salmon_foxo-inhibition/")
samples <-read.table("samples.txt")
txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id", "gene_id")])
files <- row.names(samples)
all(file.exists(files))
txi <- tximport(files, type ="salmon", tx2gene = tx2gene, ignoreTxVersion = T)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
dds <- ddsTxi
keep <- rowSums(counts(dds)) >= 3
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "DMSO")
dds <- DESeq(dds)

counts <- counts(dds, normalized=F)
write.table(counts, "deseq2_raw-counts.csv")
counts.df <- data.frame(counts)
relevant.counts <- counts.df %>% filter_all(all_vars(.>3))

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
all_coding_genes <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), 
                          values= row.names(relevant.counts), 
                          filters = "ensembl_gene_id", 
                          mart = mart, 
                          uniqueRows = TRUE)

rel.counts.hgnc <- merge(all_coding_genes,relevant.counts, by.x="ensembl_gene_id", by.y=0)
rel.counts.hgnc <- rel.counts.hgnc[,-1]
rel.counts.hgnc2 <- rel.counts.hgnc[ !(rel.counts.hgnc$hgnc_symbol == ""),]
write.csv(rel.counts.hgnc2, "raw-counts_hgnc-symbol.csv", sep ="\t", quote = F, col.names = F )

dds
res <-  results(dds, contrast=c("condition","DMSO","Inhibitor"))

resultsNames(dds)

write.table(res, "deseq2_DMSO-Inhibitor.salmon.csv")

normalized_counts <- counts(dds, normalized=T)
write.table(normalized_counts, "deseq2_normalized-counts_salmon.csv")


res.df <- as.data.frame(res)

res.hgnc <- merge(res.df, all_coding_genes, by.x=0, by.y="ensembl_gene_id")
normalized_counts_hgnc <- merge(normalized_counts, all_coding_genes, by.x=0, by.y="ensembl_gene_id")

write.table(res.hgnc, "deseq2_DMSO-Inhibitor.salmon.hgnc-symbol.csv")
write.table(normalized_counts_hgnc, "deseq2_normalized-counts_salmon_hgnc.csv")


###knock-down sample are quite heterogeneous -> not a lot of DEGs
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "run"))

##model batch effect of biological donor
ddsTxi_2 <- DESeqDataSetFromTximport(txi,
                                     colData = samples,
                                     design = ~ condition+ donor)
samples$donor <- factor(samples$donor)
ddsTxi_2 <- DESeqDataSetFromTximport(txi,
                                     colData = samples,
                                     design = ~ condition+ donor)
dds_2 <- ddsTxi_2
keep <- rowSums(counts(dds_2)) >= 3
dds_2 <- dds_2[keep,]
dds_2$condition <- relevel(dds_2$condition, ref = "DMSO")
dds_2 <- DESeq(dds_2)
counts_2 <- counts(dds_2, normalized=T)
counts_2.df <- data.frame(counts_2)
relevant.counts_2 <- counts_2.df %>% filter_all(all_vars(.>3))
rel.counts_2.hgnc <- merge(all_coding_genes,relevant.counts_2, by.x="ensembl_gene_id", by.y=0)
rel.counts_2.hgnc <- rel.counts_2.hgnc[,-1]
rel.counts_2.hgnc2 <- rel.counts_2.hgnc[ !(rel.counts_2.hgnc$hgnc_symbol == ""),]
res.2 <-  results(dds_2, contrast=c("condition","DMSO","Inhibitor"))

res.2.df <- as.data.frame(res.2)

res.2.hgnc <- merge(res.2.df, all_coding_genes, by.x=0, by.y="ensembl_gene_id")

write.table(res.2.hgnc, "deseq2_DMSO-inhibitor.salmon.hgnc-symbol_donor-normalised.csv")
res.2.hgnc <- read.table("deseq2_DMSO-inhibitor.salmon.hgnc-symbol_donor-normalised.csv", sep = " ")
table(res.2.hgnc$padj < 0.05 & res.2.hgnc$log2FoldChange >0)


###combatseq, a bit messy beacuase it was reconstructed from log after accidenbtial deletion
library(sva)

counts <- as.data.frame(txi$counts)
#Single count matrices for C2c12 and primary myoblasts; set batch factors
names(counts) <- samples$run
batch<- c(1,1,2,2,3,3)
#adjust counts using Combat-Seq from sva package
adjusted_counts <- ComBat_seq(counts, batch = batch, full_mod=TRUE)
row.names(samples) <- samples $run
dds_c2 <- DESeqDataSetFromMatrix(countData = round(adjusted_counts), colData = samples, design = ~ condition)
dds_c2 <- DESeq(dds_c2)
dds_t_c2 <- varianceStabilizingTransformation(dds_c2)
plotPCA(dds_t_c2, intgroup = "condition", ntop = 500, returnData = FALSE)
plotPCA(dds_t_c2, intgroup = c("donor"))
res2.combat <- results(dds_c2, contrast=c("condition","DMSO","Inhibitor"))
res2.combat.df <- as.data.frame(res2.combat)
res2.combat.hgnc <- merge(res2.combat.df, all_coding_genes, by.x=0, by.y="ensembl_gene_id")
write.table(res2.combat.hgnc, "deseq2_DMSO-Inhibitor.hgnc-symbol_donor-normalised_combatSeq.csv")
table(res2.combat.hgnc$padj < 0.05)

###visualisation
res.2.hgnc <- read.csv("deseq2_DMSO-inhibitor.salmon.hgnc-symbol_donor-normalised.csv", header = T, sep = " ")
sig <- res.2.hgnc%>%
  mutate(regulation = factor(case_when(padj <0.05 & log2FoldChange > 0 ~ "significantly up",
                                       padj <0.05 & log2FoldChange < 0 ~ "significantly down",
                                       TRUE ~ "not significant")))
sig_up <- sig[sig$regulation == "significantly up",]
sig_down <- sig[sig$regulation == "significantly down",]
write.table(sig_up, "deseq2_DMSO-inhibitor.salmon.hgnc-symbol_donor-normalised_sig-up.csv")
write.table(sig_down, "deseq2_DMSO-inhibitor.salmon.hgnc-symbol_donor-normalised_sig-down.csv")

gene_names <- c("CDH2", "PDGFA", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "PDGFB", "ACTA2", "COL1A1", "COL3A1", "MMP2", "MMP9", "TAGLN")
endmt_genes <- sig %>% dplyr::filter(hgnc_symbol %in% gene_names)
endmt_genes <- endmt_genes[endmt_genes$regulation != "not significant",]
p1 <- ggplot(sig, aes(x = log2FoldChange, y =-log(padj, 10),color=regulation, label = hgnc_symbol))+
  ggtitle(label = "", subtitle = "Expression of EndMT genes upon FOXO1 inhibition") +
  geom_point() +
  theme_bw(base_size = 12) + 
  #xlim(-10, 10)+
  theme(legend.position = "right") +
  xlab(expression(log[2]("FoldChange"))) + 
  ylab(expression(-log[10]("pvalue"))) +
  theme(text =element_text(size=18))+
  theme(aspect.ratio = 1)+
  scale_color_manual(labels = c("not significant", "significantly up","significantly down"), values=c("grey","#3399ff", "#ff3333"))+
  geom_text_repel(data=endmt_genes, mapping = aes(log2FoldChange, -log10(padj), label=hgnc_symbol), size =5,max.overlaps =40, color = "black", segment.color = "black", position = position_nudge_keep(x = 0.7))+
  NULL

p1


###heatmap EndMT genes
gene_names <- c("EDN1", "ID3", "CCN2", "SERPINE1", "KLF4", "KLF2", "LDLR", "THBD", "CYP1B1", "ID1", "ID2", "VIM", "CDH2", "ACTA2")
atheroprone_genes <- rel.counts_2.hgnc %>% dplyr::filter(hgnc_symbol %in% gene_names)
row.names(atheroprone_genes) <- atheroprone_genes[,1]
X <- atheroprone_genes[,-1]
names(X) <- samples$run

# make RPKMs to Z-score for proper heatmap visualization
xt = X
tX<-t(xt)
z_before_transpose_1 = scale(tX)
z_score_1 = t(z_before_transpose_1)
heatmap_df = as.data.frame(z_score_1)


# make heatmap YAY!
library(pheatmap)
my.colors = c(colorRampPalette(colors = c("blue", "white", "red"))(50))
hm = pheatmap(heatmap_df,main = "GO: negative_regulation_of_Inflammation", color = my.colors ,border_color = "black", fontsize_col = 8, cellheight = 9,cellwidth = 9, fontsize_row = 7,show_rownames=T,legend  = T, cluster_rows = T, cluster_cols = T)
print(hm)


###Gene ontology analysis
library(stringr)
up_KEGG <- read.table(file = "GO_analysis_upregulated_KEGG_for_plot.csv", sep = ";", header = T, stringsAsFactors = F)
names(up_KEGG)[10] <- "Enrichment"
up_KEGG$Enrichment <- as.numeric(up_KEGG$Enrichment)
up_KEGG$Term = str_wrap(up_KEGG$Term, width = 30)

up_KEGG$Term <- factor(up_KEGG$Term, levels = up_KEGG$Term[order(up_KEGG$Enrichment)])

ggplot(data = up_KEGG, aes(x = Enrichment, y = Term, color = FDR, size = 10)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 18) + 
  theme(aspect.ratio = 1)+
  ylab("") + 
  xlab("Enrichment") + 
  NULL

up_GO_BP <- read.table(file = "GO_analysis_upregulated_GO_BP_for_plot.csv", sep = ";", header = T, stringsAsFactors = F)
names(up_GO_BP)[10] <- "Enrichment"
up_GO_BP <- up_GO_BP[-3,]
up_GO_BP$Enrichment <- as.numeric(up_GO_BP$Enrichment)
up_GO_BP$Term = str_wrap(up_GO_BP$Term, width = 40)

up_GO_BP$Term <- factor(up_GO_BP$Term, levels = up_GO_BP$Term[order(up_GO_BP$Enrichment)])

ggplot(data = up_GO_BP, aes(x = Enrichment, y = Term, color = FDR, size = 10)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 18) + 
  theme(aspect.ratio = 1)+
  ylab("") + 
  xlab("Enrichment") + 
  NULL


##for downregulated
GO_analysis_downregulated <- read.delim(file = "GO_analysis_downregulated.txt", sep = "\t", header = T, stringsAsFactors = F)
names(GO_analysis_downregulated)[10] <- "Enrichment"
GO_analysis_downregulated$Enrichment <- as.numeric(GO_analysis_downregulated$Enrichment)
GO_analysis_downregulated$Term = str_wrap(GO_analysis_downregulated$Term, width = 40)
GO_analysis_downregulated_bp <- GO_analysis_downregulated[GO_analysis_downregulated$Category == "GOTERM_BP_DIRECT",]
GO_analysis_downregulated_bp <- GO_analysis_downregulated_bp[1:5,]
GO_analysis_downregulated_bp$Term <- sub(".*~", "", GO_analysis_downregulated_bp$Term)
GO_analysis_downregulated_bp$Term <- factor(GO_analysis_downregulated_bp$Term, levels = GO_analysis_downregulated_bp$Term[order(GO_analysis_downregulated_bp$Enrichment)])
ggplot(data = GO_analysis_downregulated_bp, aes(x = Enrichment, y = Term, color = PValue, size = 10)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 18) + 
  theme(aspect.ratio = 1)+
  ylab("") + 
  xlab("Enrichment") + 
  NULL

GO_analysis_downregulated_KEGG <- GO_analysis_downregulated[GO_analysis_downregulated$Category == "KEGG_PATHWAY",]
GO_analysis_downregulated_KEGG <- GO_analysis_downregulated_KEGG[1:5,]
GO_analysis_downregulated_KEGG$Term <- sub(".*:", "", GO_analysis_downregulated_KEGG$Term)
GO_analysis_downregulated_KEGG$Term <- factor(GO_analysis_downregulated_KEGG$Term, levels = GO_analysis_downregulated_KEGG$Term[order(GO_analysis_downregulated_KEGG$Enrichment)])
ggplot(data = GO_analysis_downregulated_KEGG, aes(x = Enrichment, y = Term, color = PValue, size = 10)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 18) + 
  theme(aspect.ratio = 1)+
  ylab("") + 
  xlab("Enrichment") + 
  NULL
