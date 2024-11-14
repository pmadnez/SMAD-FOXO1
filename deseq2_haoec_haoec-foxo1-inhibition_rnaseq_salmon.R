library(Gen)
library(tximport)
library(DESeq2)
library(readr)
library(edgeR)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library('biomaRt')
library(ggpp)
library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("Z:/data/")
### FASTQ files have been aligned and quantified with SALMON before. A samplesheet containing metadata was created seperately. Please refer to the salmon vignette for more information
samples <-read.table("samples.txt")
txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id", "gene_id")])
files <- row.names(samples)
all(file.exists(files))
txi <- tximport(files, type ="salmon", tx2gene = tx2gene, ignoreTxVersion = T)
##This is without Donor batch correction (not used in manuscript, just for completeness)
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
res <-  results(dds, contrast=c("condition","DMSO","Inhibitor")) ##this needs to be changed according to the samples present in RNAseq experiment

resultsNames(dds)

write.table(res, "deseq2_DMSO-Inhibitor.salmon.csv")

normalized_counts <- counts(dds, normalized=T)
write.table(normalized_counts, "deseq2_normalized-counts_salmon.csv")


res.df <- as.data.frame(res)

###get hgnc symbols additionally to ENSEMBL
res.hgnc <- merge(res.df, all_coding_genes, by.x=0, by.y="ensembl_gene_id")
normalized_counts_hgnc <- merge(normalized_counts, all_coding_genes, by.x=0, by.y="ensembl_gene_id")

write.table(res.hgnc, "deseq2_DMSO-Inhibitor.salmon.hgnc-symbol.csv")
write.table(normalized_counts_hgnc, "deseq2_normalized-counts_salmon_hgnc.csv")


##This is with Donor batch correction
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
