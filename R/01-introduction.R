# Notes for 01-introduction.R
# --------------------------------------
## Copy code from https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/01-introduction.R

## Notes

## ----'quick_intro_01', message = FALSE-------------------------------------------------------------------------------
library('scRNAseq')
library('scater')
library('scran')
library('plotly')


## ----'quick_intro_02', cache = TRUE----------------------------------------------------------------------------------
sce <- scRNAseq::MacoskoRetinaData()

## How big is the data?
pryr::object_size(sce)

## How does it look?
sce

## ----'quick_intro_03', cache = TRUE----------------------------------------------------------------------------------
# Quality control.
is.mito <- grepl("^MT-", rownames(sce))
qcstats <-
    scater::perCellQCMetrics(sce, subsets = list(Mito = is.mito))
filtered <-
    scater::quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- scater::logNormCounts(sce)

# Feature selection.
dec <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(dec, prop = 0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- scater::runPCA(sce, ncomponents = 25, subset_row = hvg)
sce <- scater::runUMAP(sce, dimred = 'PCA', external_neighbors = TRUE)

# Clustering.
g <- scran::buildSNNGraph(sce, use.dimred = 'PCA')
sce$clusters <- factor(igraph::cluster_louvain(g)$membership)


## ----'quick_intro_04'------------------------------------------------------------------------------------------------
# Visualization.
scater::plotUMAP(sce, colour_by = "clusters")


scater::plotUMAP(sce, colour_by = "clusters")

## # If I'm working through the cluster connection, then I need to

pdf('plot_umap.pdf')
scater::plotUMAP(sce, colour_by = "clusters")

dev.off()

## ----'quick_intro_05', eval = FALSE----------------------------------------------------------------------------------
## # Interactive visualization

# This makes the map interactive, showing little notes
p <- scater::plotUMAP(sce, colour_by = "clusters")
plotly::ggplotly(p)
