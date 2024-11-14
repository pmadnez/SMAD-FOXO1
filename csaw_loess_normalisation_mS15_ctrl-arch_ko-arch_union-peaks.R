#example experimental design: n=2/3 mouse ATAC-seq biological replicates for two conditions: treat and control 
#This approach can be adjusted to whatever conditions you might have (change number of replicates etc.)
#peak files come from the ENCODE pipeline
library(GenomicRanges)
library(csaw)

########################################
########################################
########################################

# starting from MACS2 filtered broadpeaks
# read replicate broadPeak files
setwd("Z:/data/")
a <- read.table("samplesheet_new.csv", sep = ",", header = T)

treat1.peaks <- read.table(a$Peaks[8], sep="\t")[,1:3]
treat2.peaks <- read.table(a$Peaks[9], sep="\t")[,1:3]
treat3.peaks <- read.table(a$Peaks[10], sep="\t")[,1:3]
control1.peaks <- read.table(a$Peaks[3], sep="\t")[,1:3]
control2.peaks <- read.table(a$Peaks[4], sep="\t")[,1:3]
colnames(treat1.peaks) <- c("chrom", "start", "end")
colnames(treat2.peaks) <- c("chrom", "start", "end")
colnames(treat3.peaks) <- c("chrom", "start", "end")
colnames(control1.peaks) <- c("chrom", "start", "end")
colnames(control2.peaks) <- c("chrom", "start", "end")


# convert to GRanges objects
treat1.peaks <- GRanges(treat1.peaks)
treat2.peaks <- GRanges(treat2.peaks)
treat3.peaks <- GRanges(treat2.peaks)
control1.peaks <- GRanges(control1.peaks)
control2.peaks <- GRanges(control2.peaks)

# define consensus peakset

#one method: union of all replicate peak sets for both conditions
treat.peaks <- union(treat1.peaks, treat2.peaks)#, treat3.peaks)
treat.peaks <- union(treat.peaks, treat3.peaks)
control.peaks <- union(control1.peaks, control2.peaks)
all.peaks <- union(treat.peaks, control.peaks)




##############################
# specify paired-end BAMs
pe.bams <- c(a$bamReads[3], a$bamReads[4],
             a$bamReads[8], a$bamReads[9], a$bamReads[10])
##############################
# read mm10 blacklist
blacklist <- read.table("mm10_blacklist.bed", sep="\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)

# define read parameters
standard.chr <- paste0("chr", c(1:19, "X", "Y")) # only use standard chromosomes
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)

##############################
# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)

##############################
# MACS2 peaks only: filter low abundance peaks
library("edgeR")
peak.abundances <- aveLogCPM(asDGEList(peak.counts)) 
peak.counts.filt <- peak.counts[peak.abundances > -3, ] # only use peaks logCPM > -3
# few or no peaks should be removed; modify as desired

##############################

# get paired-end fragment size distribution
control1.pe.sizes <- getPESizes(a$bamReads[3])
control2.pe.sizes <- getPESizes(a$bamReads[4])
treat1.pe.sizes <- getPESizes(a$bamReads[8])
treat2.pe.sizes <- getPESizes(a$bamReads[9])
treat3.pe.sizes <- getPESizes(a$bamReads[10])
#gc()
# plot
hist(control1.pe.sizes$sizes) # repeat for all replicates and conditions

# for analysis with csaw de novo enriched query windows, select a window size that is greater than the majority of fragments

##############################
# count BAM reads in, e.g. 300 bp windows
counts <- windowCounts(pe.bams, width=300, param=param) # set width as desired from the fragment length distribution analyses

# filter uninteresting features (windows) by local enrichment
# local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods
# filter.stat <- filterWindows(counts, wider, type="local") # the filterWindows() function is deprecated and has been replaced by filterWindowsLocal(). This is an archived step.
filter.stat <- filterWindowsLocal(counts, wider)
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired

###############################
# count BAM background bins (for TMM normalization)
binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)

##########################################
# NORMALIZATION

# MACS2 peaks only, csaw loess-normalization
peak.counts.loess <- peak.counts.filt
peak.counts.loess <- normOffsets(peak.counts.loess, se.out=TRUE) # type="loess" is now default
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

# DIFFERENTIAL ACCESSIBILITY ANALYSIS

# set working windows for the desired analysis
working.windows <- peak.counts.loess # MACS2 peaks only, for trended biases
###########

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)
colnames(y$counts) <- c("control1", "control2", "treat1", "treat2", "treat3")
rownames(y$samples) <- c("control1", "control2", "treat1", "treat2", "treat3")
y$samples$group <- c("control", "control", "treat", "treat", "treat")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("control", "treat") # CONFIRM THAT THESE COLUMNS CORRECTLY ALIGN!!
design
# IMPORTANT: the user should manually confirm that the design matrix is correctly labeled according to sample metadata!

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))
rowData(working.windows) <- cbind(rowData(working.windows), results$table) # combine GRanges rowdata with differential statistics


# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)

# concatenating all relevant statistical data for final merged windows (no redundant columns)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks

# filter by FDR threshold
FDR.thresh <- 0.2 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows
results_df <- as.data.frame(final.merged.peaks.sig)

write.table(final.merged.peaks, "arch.ctrl_vs_arch.ko_csaw_DA-windows_all.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(final.merged.peaks.sig, "arch.ctrl_vs_arch.ko_csaw_DA-windows_significant.txt", sep="\t", quote=F, col.names=T, row.names=F)

###########################################


####################################
##HOMER
#get bed files for HOMER motif analyis
##up means up in arch, as I had it the other way before its changed here
library((tibble))


bed_up <- results_df[results_df$direction == "down", c(1,2,3,5)]
bed_down <- results_df[results_df$direction == "up", c(1,2,3,5)]

up_bed <- bed_up %>%
  tibble::rownames_to_column("id") %>%
  dplyr::mutate(score = 1000) %>%
  dplyr::select(seqnames, start, end, id, score, strand) %>%
  dplyr::mutate(strand = ".")
down_bed <- bed_down %>%
  tibble::rownames_to_column("id") %>%
  dplyr::mutate(score = 1000) %>%
  dplyr::select(seqnames, start, end, id, score, strand) %>%
  dplyr::mutate(strand = ".")

write.table(up_bed, "csaw_loess_FDR.01_arch.ctrl_vs_arch.ko_up.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(down_bed, "csaw_loess_FDR.01_arch.ctrl_vs_arch.ko_down.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
