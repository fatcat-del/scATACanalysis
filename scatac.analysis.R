library(GenomicRanges)
library(S4Vectors)
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(20, 70, 300), end = c(120, 200, 400)))
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(AtacAnnoR)
library(harmony)
library(EnsDb.Mmusculus.v79)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)



colorset <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72","#B17BA6","#FF7F00",
              "#FDB462","#E7298A","#E78AC3","#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02",
              "#7570B3","#BEAED4","#666666","#999999","#aa8282","#d4b7b7","#8600bf","#ba5ce3",
              "#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d",
              "#ffff00", "#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#BBBBBB",
              "#027EB6", "#746FB2", "#9651A0", "#C8008F", "#ee64a4",
              "#EE0220", "#D93F00", "#795549", "#6F7C4D", "#008A25",
              "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
              "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

peaks.control <- read.table(
  file = "/data/fengw/scATAC/control-1/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.HF <- read.table(
  file = "/data/fengw/scATAC/HF-1/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.HFAT <- read.table(
  file = "/data/fengw/scATAC/HFAT/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end")
)
gr.control <- makeGRangesFromDataFrame(peaks.control)
gr.HF <- makeGRangesFromDataFrame(peaks.HF)
gr.HFAT <- makeGRangesFromDataFrame(peaks.HFAT)
combined.peaks <- reduce(x = c(gr.control, gr.HF, gr.HFAT))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


md.control <- read.table(
  file = "/data/fengw/scATAC/control-1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.HF <- read.table(
  file = "/data/fengw/scATAC/HF-1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.HFAT <- read.table(
  file = "/data/fengw/scATAC/HFAT/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

#md.control <- md.control[md.control$passed_filters > 1000, ]
#md.HF <- md.HF[md.HF$passed_filters > 1000, ]
#md.HFAT <- md.HFAT[md.HFAT$passed_filters > 1000, ]

control.cell <- read.csv("/data/fengw/scATAC/control-1/outs/filtered_peak_bc_matrix/barcodes.tsv",header=F)
md.control <- md.control[control.cell$V1, ]

HF.cell <- read.csv("/data/fengw/scATAC/HF-1/outs/filtered_peak_bc_matrix/barcodes.tsv",header=F)
md.HF <- md.HF[HF.cell$V1, ]

HFAT.cell <- read.csv("/data/fengw/scATAC/HFAT/outs/filtered_peak_bc_matrix/barcodes.tsv",header=F)
md.HFAT <- md.HFAT[HFAT.cell$V1, ]


frags.control <- CreateFragmentObject(
  path = "/data/fengw/scATAC/control-1/outs/fragments.tsv.gz",
  cells = rownames(md.control)
)
frags.HF <- CreateFragmentObject(
  path = "/data/fengw/scATAC/HF-1/outs/fragments.tsv.gz",
  cells = rownames(md.HF)
)
frags.HFAT <- CreateFragmentObject(
  path = "/data/fengw/scATAC/HFAT/outs/fragments.tsv.gz",
  cells = rownames(md.HFAT)
)


control.counts <- FeatureMatrix(
  fragments = frags.control,
  features = combined.peaks,
  cells = rownames(md.control)
)
HF.counts <- FeatureMatrix(
  fragments = frags.HF,
  features = combined.peaks,
  cells = rownames(md.HF)
)

HFAT.counts <- FeatureMatrix(
  fragments = frags.HFAT,
  features = combined.peaks,
  cells = rownames(md.HFAT)
)
control_assay <- CreateChromatinAssay(control.counts, fragments = frags.control)
control <- CreateSeuratObject(control_assay, assay = "ATAC", meta.data=md.control)
HF_assay <- CreateChromatinAssay(HF.counts, fragments = frags.HF)
HF <- CreateSeuratObject(HF_assay, assay = "ATAC", meta.data=md.HF)
HFAT_assay <- CreateChromatinAssay(HFAT.counts, fragments = frags.HFAT)
HFAT <- CreateSeuratObject(HFAT_assay, assay = "ATAC", meta.data=md.HFAT)
control$orig.ident="control"
HF$orig.ident="HF"
HFAT$orig.ident="HFAT"
combined <- merge(
  x = control,
  y = c(HFAT, HF),
  add.cell.ids = c("control", "HFAT","HF")
)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
Annotation(combined)=annotation

combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined <- TSSEnrichment(combined, fast = FALSE)
combined <- NucleosomeSignal(combined)

combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')

combined$blacklist_fraction <- FractionCountsInRegion(
  object = combined,
  assay = 'ATAC',
  regions = blacklist_mm10
)


p1 <- VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 3
)


combined <- subset(
    x = combined,
    subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 20 &
    pct_reads_in_peaks < 60 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1 &
    TSS.enrichment < 10
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

combined <- FindNeighbors(
  object = combined,
  reduction = 'lsi',
  dims = 2:30
)
combined <- FindClusters(
  object = combined,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

p1 <- DimPlot(object = combined, cols = colorset, label = TRUE, split.by = "orig.ident") + NoLegend()
ggsave("./umap.split.png", p1)


gene.activities <- GeneActivity(combined)
sce <- CreateSeuratObject(gene.activities, min.cells=3, min.features=200, project="activte")

combined[['RNA']] <- CreateAssayObject(counts = sce@assays$RNA@counts)

p1 <- VlnPlot(combined, features=c("nCount_RNA", "nFeature_RNA"), group.by = "orig.ident")
ggsave("./qc.rna.split.png", p1)

combined <- SCTransform(combined, assay = "RNA",vars.to.regress = c("nCount_RNA"))
combined <- NormalizeData(combined, assay = "RNA")
combined <- RunPCA(combined, assay ="SCT", npcs=50, verbose=FALSE)
combined <- RunHarmony(combined, assay.use = "SCT",group.by.vars="orig.ident", max.iter.harmony = 20)
combined <- RunUMAP(combined,reduction="harmony", dims=1:50 ,seed.use = 0)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:50)
combined <- FindClusters(combined, resolution = 0.2, random.seed = 0)

p1 <- DimPlot(combined, cols = colorset, pt.size = 2, label = T, split.by = "orig.ident")
ggsave("./umap.rna.split.png", p1)



sce <- readRDS("/data/fengw/singlecell/rna_scvelo/liversc_clean.rds")
sce <- sce[,sce$orig.ident %in% c("control", "HFAT")]
pred <- RunAtacAnnoR(ref_mtx = sce[['RNA']]@counts, 
                     ref_celltype = sce$celltype, 
                     query_gene_activity = combined[['RNA']]@counts, 
                     query_peak_counts = combined[['ATAC']]@counts) 


pred <- RunAtacAnnoR_Signac(query_SeuratObj = sce,
                            query_ga_assay = 'RNA',
                            query_peak_assay = 'ATAC',
                            ref_SeuratObj = sce,
                            ref_assay = 'RNA',
                            ref_ident = 'celltype')
							



da_peaks_HFAT <- FindMarkers(
  object = hep,
  group.by="orig.ident",
  ident.1 = 'HFAT',
  ident.2 = 'control',
  only.pos = F,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)




pdf("PPARA_test0715.pdf")
CoveragePlot(
  object = combined,
  group.by = 'orig.ident',
  region = "Saa1",
  assay = "ATAC"
)
dev.off()


CoveragePlot(
  object = brain,
  region = c("Neurod6", "Gad2"),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)

CoveragePlot(
  object = hep, assay = "ATAC",
  region = "chr12-103944030-103944939",
  extend.upstream = 10000,
  extend.downstream = 10000,
  group.by = "orig.ident"
)


da_peaks_HFAT <- FindMarkers(
  object = hep,
  group.by="orig.ident",
  ident.1 = 'HFAT',
  ident.2 = 'control',
  only.pos = F,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)

peaks_anno <- ClosestFeature(hep, regions=rownames(da_peaks))

da_peaks$genes <- ifelse(rownames(da_peaks) %in% peaks_anno$query_region, peaks_anno$gene_name, "N")














dfanno <- ClosestFeature(kf, regions = rownames(da_peaks_HFAT))
rownames(dfanno) <-dfanno$query_region
temp <- da_peaks_HFAT[rownames(dfanno),]
aaa <- cbind(dfanno,temp)


query_nmf_embedding <- get_nmf_embedding(peak_counts = combined[['ATAC']]@counts)


pre_processing_mtxs <- pre_processing(ref_mtx = sce[['RNA']]@counts,
                                      query_mtx = combined[['ACTIVITY']]@counts,
                                      verbose = T)
									  
ref_mtx <- pre_processing_mtxs$ref_mtx
query_mtx <- pre_processing_mtxs$query_mtx

global_markers <- get_global_markers_sc(sc_counts_mtx = ref_mtx,
                                        labels = sce$celltype,
                                        max_marker = 200)
neighbor_celltypes <- get_neighbor_celltypes(sc_count_mtx = ref_mtx,
                                             labels = sce$celltype,
                                             global_markers,
                                             min_cor = 0.6)

neighbor_markers <- get_neighbor_markers_sc(sc_counts_mtx = ref_mtx,
                                            labels = sce$celltype,
                                            neighbor_celltypes = neighbor_celltypes,
                                            global_markers = global_markers)
																						
											
plot_ref_global_markers_heatmap(ref_mtx = ref_mtx,
                                ref_labels = sce$celltype,
                                global_markers = global_markers,
                                neighbor_celltypes = neighbor_celltypes)
								
								
								
plot_ref_neighbor_markers_heatmap(ref_mtx = ref_mtx,
                                  ref_labels = sce$celltype,
                                  neighbor_markers = neighbor_markers,
                                  neighbor_celltypes = neighbor_celltypes,
                                  celltype_to_plot = 'Hepatocytes')
								  
cor_mtx <- get_cor_mtx(sc_count_mtx = ref_mtx,
                       labels = sce$celltype,
                       query_mtx = query_mtx,
                       global_markers = global_markers,
                       query_nmf_embedding = query_nmf_embedding)

cell_meta <- get_kendall_pred(cor_mtx)

cell_meta <- test_markers(query_mtx,cell_meta,global_markers,neighbor_markers,threads = 10,verbose = T)

cell_meta <- get_seed_candidates(cell_meta)


transfer.anchors <- FindTransferAnchors(  
  reference = sce,
  query = combined,
  reduction = 'cca')
