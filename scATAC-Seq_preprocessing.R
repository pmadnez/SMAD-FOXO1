library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(magrittr)

# define a convenient function to load all the data and create a Seurat object
create_obj <- function(dir) {
  count.path <- list.files(path = dir, pattern = "filtered_peak_bc_matrix.h5", full.names = TRUE)
  fragment.path <- list.files(path = dir, pattern = "fragments.tsv.gz", full.names = TRUE)[1]
  counts <- Read10X_h5(count.path)
  md.path <- list.files(path = dir, pattern = "singlecell.csv", full.names = TRUE)
  md <- read.table(file = md.path, stringsAsFactors = FALSE, sep = ",", header = TRUE)
  scATAC <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), genome = "mm9", fragments = fragment.path, min.cells = 10, min.features = 200)
  obj <- CreateSeuratObject(counts = scATAC, assay = "peaks", project = "ATAC", meta.data = md)
  return(obj)
}

##outs directories come from cellranger atac output, please refer to cellranger/ 10x manual/ vignette
dr <- create_obj("Z:/scATACseq_mm_2dR/outs/")
dl <- create_obj("Z:/scATACseq_mm_2dL/outs/")
wr <- create_obj("Z:/scATACseq_mm_2wR/outs/")
wl <- create_obj("Z:/scATACseq_mm_2wL/outs/")



# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm9"

# add the gene information to the object
Annotation(dr) <- annotations
Annotation(dl) <- annotations
Annotation(wr) <- annotations
Annotation(wl) <- annotations

#compute gene activities

gene.activities.dr <- GeneActivity(dr)
gene.activities.dl <- GeneActivity(dl)
gene.activities.wr <- GeneActivity(wr)
gene.activities.wl <- GeneActivity(wl)

dr[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities.dr)
dl[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities.dl)
wr[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities.wr)
wl[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities.wl)

dr <- subset(dr, subset = nCount_peaks > 5000)
dl <- subset(dl, subset = nCount_peaks > 5000)
wr <- subset(wr, subset = nCount_peaks > 5000)
wl <- subset(wl, subset = nCount_peaks > 5000)

dr$orig.ident <- factor(x = dr@meta.data$orig.ident, levels = "2D-R")
dl$orig.ident <- factor(x = dl@meta.data$orig.ident, levels = "2D-L")
wr$orig.ident <- factor(x = wr@meta.data$orig.ident, levels = "2W-R")
wl$orig.ident <- factor(x = wl@meta.data$orig.ident, levels = "2W-L")

dr$active.ident <- factor(x = dr@active.ident, levels = "2D-R")
dl$active.ident <- factor(x = dl@active.ident, levels = "2D-L")
wr$active.ident <- factor(x = wr@active.ident, levels = "2W-R")
wl$active.ident <- factor(x = wl@active.ident, levels = "2W-L")

dl.2 <- dl
dr.2 <- dr
wr.2 <- wr
wl.2 <- wl
dl.2$orig.ident <- "2D-L"
dr.2$orig.ident <- "2D-R"
wr.2$orig.ident <- "2W-R"
wl.2$orig.ident <- "2W-L"

#merge datasets



pbmc.atac.2 <- merge(dr.2, y = c(dl.2, wr.2, wl.2), 
                     add.cell.ids = c("2D-R", "2D-L", "2W-R", "2W-L"), project = "atac")


pbmc.2 <- pbmc.atac.2

pbmc.2 <- RunTFIDF(pbmc.2)
pbmc.2 <- FindTopFeatures(pbmc.2, min.cutoff = 20)
pbmc.2 <- RunSVD(
  pbmc.2,
  reduction.key = 'LSI_',
  reduction.name = 'lsi', 
  irlba.work = 400
)

DepthCor(pbmc.2)

pbmc.2 <- RunUMAP(object = pbmc.2, reduction = 'lsi', dims = 2:30)
pbmc.2 <- FindNeighbors(object = pbmc.2, reduction = 'lsi', dims = 2:30)
pbmc.2 <- FindClusters(object = pbmc.2, verbose = FALSE, algorithm = 3, resolution = 0.2)
DimPlot(object = pbmc.2, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(pbmc.2)
pbmc.2[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc.2 <- NormalizeData(
  object = pbmc.2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc.2$nCount_RNA)
)

DefaultAssay(pbmc.2) <- 'RNA'

DotPlot(pbmc.2, features = c("Itk", "Mmp25", "Ccr7","C5ar1", "C1qc", "C1qb", "C1qa", "Lum", "Dpep1", "Medag", "Speg", "Myh11", "Cnn1","Tie1","Icam2","Cdh5", "Pecam1"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Create factors to have levels in orig.ident
my_levels <- c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7", "SMC1","SMC2", "Fibro", "Mo1", "Mo2", "Mo3", "Mo4", "DC", "T")
pbmc.2@active.ident <- factor(x = pbmc.2@active.ident, levels = my_levels)


ec <- subset(pbmc.2, idents = c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7"))

