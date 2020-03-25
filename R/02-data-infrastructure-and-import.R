# Notes for 02-data-infrastructure-and-import.R
# --------------------------------------
## Copy code from https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/02-data-infrastructure-and-import.R

## Notes

## ----setup, include=FALSE--------------------------------------------------------------------------------------------
options(htmltools.dir.version = FALSE)


## ----xaringan-themer, include=FALSE----------------------------------------------------------------------------------
library(xaringanthemer)
solarized_dark(
    code_font_family = "Fira Code",
    code_font_url    = "https://cdn.rawgit.com/tonsky/FiraCode/1.204/distr/fira_code.css"
)


## /* From https://github.com/yihui/xaringan/issues/147  */

## .scroll-output {

##   height: 80%;

##   overflow-y: scroll;

## }

##
## /* https://stackoverflow.com/questions/50919104/horizontally-scrollable-output-on-xaringan-slides */

## pre {

##   max-width: 100%;

##   overflow-x: scroll;

## }

##
## /* From https://github.com/yihui/xaringan/wiki/Font-Size */

## .tiny{

##   font-size: 40%

## }

##
## /* From https://github.com/yihui/xaringan/wiki/Title-slide */

## .title-slide {

##   background-image: url(https://raw.githubusercontent.com/Bioconductor/OrchestratingSingleCellAnalysis/master/images/Workflow.png);

##   background-size: 33%;

##   background-position: 0% 100%

## }


## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------

## # We starded here:
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")

# Load the SingleCellExperiment package
library('SingleCellExperiment')
# Extract the count matrix from the 416b dataset
counts.416b <- counts(sce.416b)
# Construct a new SCE from the counts matrix
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Inspect the object we just created
sce

## How big is it?
pryr::object_size(sce)

# Access the counts matrix from the assays slot
# WARNING: This will flood RStudio with output!

# 1. The general method
assay(sce, "counts")[1:6, 1:3]
# 2. The special method for the assay named "counts"
counts(sce)[1:6, 1:3]

sce <- scater::logNormCounts(sce)
# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

# 1. The general method
assay(sce, "logcounts")[1:6, 1:3]
# 2. The special method for the assay named "logcounts"
logcounts(sce)[1:6, 1:3]

# assign a new entry to assays slot
assay(sce, "counts_100") <- assay(sce, "counts") + 100
# List the assays in the object
assays(sce)
assayNames(sce)

## How big is it?
pryr::object_size(sce)

# Extract the sample metadata from the 416b dataset
colData.416b <- colData(sce.416b)
# Add some of the sample metadata to our SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Inspect the object we just updated
sce
# Access the sample metadata from our SCE
colData(sce)
# Access a specific column of sample metadata from our SCE
table(sce$block)

# Example of function that adds extra fields to colData
sce <- scater::addPerCellQC(sce.416b)
# Access the sample metadata from our updated SCE
colData(sce)

# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

## Add the lognorm counts again
sce <- scater::logNormCounts(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just wild type cells
# Remember, cells are columns of the SCE
sce[, sce$phenotype == "wild type phenotype"]

# Access the feature metadata from our SCE
# It's currently empty!
rowData(sce)

# Example of function that adds extra fields to rowData
sce <- scater::addPerFeatureQC(sce)
# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)


# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Annotate each gene with its chromosome location
ensdb <- ah[["AH73905"]]

chromosome <- mapIds(ensdb,
                     keys = rownames(sce),
                     keytype = "GENEID",
                     column = "SEQNAME")

class(chromosome)
rowData(sce)$chromosome <- chromosome

# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just genes on chromosome 3
# NOTE: which() needed to cope with NA chromosome names
x <- rowData(sce)$chromosome
table(is.na(x))
y <- x == '3'
table(is.na(y))
table(is.na(which(y)))
head(which(y))
sce[which(rowData(sce)$chromosome == "3"), ]



# Access the metadata from our SCE
# It's currently empty!
metadata(sce)

# The metadata slot is Vegas - anything goes
metadata(sce) <- list(favourite_genes = c("Shh", "Nck1", "Diablo"),
                      analyst = c("Pete"))

# Access the metadata from our updated SCE
metadata(sce)

# E.g., add the PCA of logcounts
# NOTE: We'll learn more about PCA later
sce <- scater::runPCA(sce)
# Inspect the object we just updated
sce
# Access the PCA matrix from the reducedDims slot
reducedDim(sce, "PCA")[1:6, 1:3]

# E.g., add a t-SNE representation of logcounts
# NOTE: We'll learn more about t-SNE later
sce <- scater::runTSNE(sce)
# Inspect the object we just updated
sce
# Access the t-SNE matrix from the reducedDims slot
head(reducedDim(sce, "TSNE"))

# E.g., add a 'manual' UMAP representation of logcounts
# NOTE: We'll learn more about UMAP later and a
# 		  simpler way to compute it.
u <- uwot::umap(t(logcounts(sce)), n_components = 2)
# Add the UMAP matrix to the reducedDims slot
# Access the UMAP matrix from the reducedDims slot
reducedDim(sce, "UMAP") <- u

# List the dimensionality reduction results stored in # the object
reducedDims(sce)

# Extract the ERCC SCE from the 416b dataset
ercc.sce.416b <- altExp(sce.416b, "ERCC")
# Inspect the ERCC SCE
ercc.sce.416b

# Add the ERCC SCE as an alternative experiment to our SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

# List the alternative experiments stored in the object
altExps(sce)

# Subsetting the SCE by sample also subsets the
# alternative experiments
## I'm taking only the first TEN cells
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

## How big is it?
pryr::object_size(sce.subset)

# Extract existing size factors (these were added
# when we ran scater::logNormCounts(sce))
head(sizeFactors(sce))

# 'Automatically' replace size factors
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

# 'Manually' replace size factors
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))


## ----------------- sce object exercises ----------------------------------------------------------------------------
# Which function defines the sce class?
## class(sce)
## "SingleCellExperiment" from package "SingleCellExperiment"

# What are te minimum type of tables an sce object contains?
## Gene info, which can be foun on 'rowData()'
## Number of reads that overlap each gene for EACH cell, this can be accessed through 'assays'
## Info about cells, this can be accesed through 'colData()'
## Optinal tables: PCA, TSNE (reducedDims), alternative experiments (altExp), random info (stored in 'metadata')

# Where are the colnames(sce) used?
head(colnames(sce))
## Also, column names of the assays + rownames of the colData
# identical(colnames(assays(sce, 'counts')), rownames(colData(sce)))
## colnames(assays(sce, 'counts')) is NULL in this exmaple
identical(rownames(reducedDim(sce, 'PCA')), rownames(colData(sce)))

# Simlarly, where are the rownames(sce) used?
head(rownames(sce))
## rownames(assays(sce, 'counts'))
head(rownames(rowData(sce)))

# How many principal components did we compute?
## For this, we checked the dimension of the original object 'sce' and then compare with the 'PCA' object
dim(reducedDim(sce, 'PCA'))
dim(reducedDim(sce, 'PCA'))
## We have 50 PCA

# Which three chromosomes hace the highest mean gene expression?
rowData(sce)
## For this, we use 'sort' to arrange the data accordingly
sort(with(rowData(sce), tapply(mean, chromosome, base::mean)), decreasing = TRUE)
# Another option is
sort(tapply(rowData(sce)$mean, rowData(sce)$chromosome, base::mean), decreasing = TRUE)



## ---------Excersice Day 2----------------------------------------------------------------------------------------------------

# To read the table from the internet
## I'm using 'stringAsFactors = FALSE' to keep strings as strings, or I could use 'as.is = TRUE'

data <- read.delim("https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt",
                   header=TRUE,
                   sep="\t",
                   stringsAsFactors = FALSE)
head(data)

dim(data)
dim(sce) # 192 cells
dim(altExp(sce, 'ERCC')) # 92 ERCC sequences

head(rownames(data))
head(rownames(altExp(sce, 'ERCC')))

# Before matching, assing x and y
x <- rownames(altExp(sce, 'ERCC'))

y <- data$ERCC.ID

# Match x and y
object <- match(x,y)

## Check if it worked - all should be FALSE
table(is.na(object))
### To make check it automatically
if(!all(!is.na(object))) {
    stop("Hay al menos un NA")
    }

# Re-order the data
data <- data[object, ]

# To check it worked
identical(rownames(altExp(sce, "ERCC")), data$ERCC.ID)

# Normalize the ERCC counts
# altExp(sce, "ERCC") <- scater::logNormCounts(altExp(sce, "ERCC"))

i <- 1
plot(
        data[, "concentration in Mix 1 (attomoles/ul)"] ~
            counts(altExp(sce, "ERCC"))[, i]
        )
plot(data[, "concentration in Mix 1 (attomoles/ul)"] ~
         counts(altExp(sce, "ERCC"))[, i]
         xlab


)

## ----DAY 2. ercc_exercise, cache = TRUE, dependson='all_code'---------------------------------------------------------------
## Read the data from the web
ercc_info <-
    read.delim(
        'https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt',
        as.is = TRUE,
        row.names = 2,
        check.names = FALSE
    )

## Match the ERCC data
m <- match(rownames(altExp(sce, "ERCC")), rownames(ercc_info))
ercc_info <- ercc_info[m, ]

## Normalize the ERCC counts
altExp(sce, "ERCC") <- scater::logNormCounts(altExp(sce, "ERCC"))


## ----ercc_solution_plots, cache = TRUE, dependson='ercc_exercise'----------------------------------------------------
for (i in seq_len(2)) {
    plot(
        log2(10 * ercc_info[, "concentration in Mix 1 (attomoles/ul)"] + 1) ~
            log2(counts(altExp(sce, "ERCC"))[, i] +
                     1),
        xlab = "log norm counts",
        ylab = "Mix 1: log2(10 * Concentration + 1)",
        main = colnames(altExp(sce, "ERCC"))[i],
        xlim = c(min(logcounts(
            altExp(sce, "ERCC")
        )), max(logcounts(
            altExp(sce, "ERCC")
        )))
    )
    abline(0, 1, lty = 2, col = 'red')
}



## ----all_code_part2, cache=TRUE--------------------------------------------------------------------------------------
# Download example data processed with CellRanger
# Aside: Using BiocFileCache means we only download the
#        data once
library('BiocFileCache')
bfc <- BiocFileCache()
pbmc.url <-
    paste0(
        "http://cf.10xgenomics.com/samples/cell-vdj/",
        "3.1.0/vdj_v1_hs_pbmc3/",
        "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
    )
pbmc.data <- bfcrpath(bfc, pbmc.url)

# Extract the files to a temporary location
untar(pbmc.data, exdir = tempdir())

# List the files we downloaded and extracted
# These files are typically CellRanger outputs
pbmc.dir <- file.path(tempdir(),
                      "filtered_feature_bc_matrix")
list.files(pbmc.dir)

# Import the data as a SingleCellExperiment
library('DropletUtils')
sce.pbmc <- read10xCounts(pbmc.dir)
# Inspect the object we just constructed
sce.pbmc

## How big is it?
pryr::object_size(sce.pbmc)

# Store the CITE-seq data in an alternative experiment
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Inspect the object we just updated
sce.pbmc

## How big is it?
pryr::object_size(sce.pbmc)

# Download example data processed with scPipe
library('BiocFileCache')
bfc <- BiocFileCache()
sis_seq.url <-
    "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

# Extract the files to a temporary location
unzip(sis_seq.data, exdir = tempdir())

# List (some of) the files we downloaded and extracted
# These files are typical scPipe outputs
sis_seq.dir <- file.path(tempdir(),
                         "SIS-seq_script-master",
                         "data",
                         "BcorKO_scRNAseq",
                         "RPI10")
list.files(sis_seq.dir)

# Import the data as a SingleCellExperiment
library('scPipe')
sce.sis_seq <- create_sce_by_dir(sis_seq.dir)
# Inspect the object we just constructed
sce.sis_seq

## How big is it?
pryr::object_size(sce.sis_seq)

# Download example bunch o' files dataset
library('BiocFileCache')
bfc <- BiocFileCache()
lun_counts.url <-
    paste0(
        "https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
    )
lun_counts.data <- bfcrpath(bfc, lun_counts.url)
lun_coldata.url <-
    paste0("https://www.ebi.ac.uk/arrayexpress/files/",
           "E-MTAB-5522/E-MTAB-5522.sdrf.txt")
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

# Extract the counts files to a temporary location
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir = lun_counts.dir)

# List the files we downloaded and extracted
list.files(lun_counts.dir)

# Import the count matrix (for 1 plate)
lun.counts <- read.delim(
    file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
# Store the gene lengths for later
gene.lengths <- lun.counts$Length
# Convert the gene counts to a matrix
lun.counts <- as.matrix(lun.counts[, -1])

# Import the sample metadata
lun.coldata <- read.delim(lun_coldata.data,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)
library('S4Vectors')
lun.coldata <- as(lun.coldata, "DataFrame")

# Match up the sample metadata to the counts matrix
m <- match(colnames(lun.counts),
           lun.coldata$`Source Name`)
lun.coldata <- lun.coldata[m,]

# Construct the feature metadata
lun.rowdata <- DataFrame(Length = gene.lengths)

# Construct the SingleCellExperiment
lun.sce <- SingleCellExperiment(
    assays = list(assays = lun.counts),
    colData = lun.coldata,
    rowData = lun.rowdata
)
# Inspect the object we just constructed
lun.sce

## How big is it?
pryr::object_size(lun.sce)


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()
