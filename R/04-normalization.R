# Notes for 04-normalization.R
# --------------------------------------
## Copy code from https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/04-normalization.R

## Notes

## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
library('scRNAseq')
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)

# Quality control
library('scater')
is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
    quickPerCellQC(stats,
                   percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[, !qc$discard]


## ----all_code2, cache=TRUE, dependson='all_code'---------------------------------------------------------------------
# Library size factors
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)

class(lib.sf.zeisel)
length(lib.sf.zeisel)
head(lib.sf.zeisel)
mean(lib.sf.zeisel)

# Examine distribution of size factors
summary(lib.sf.zeisel)
hist(log10(lib.sf.zeisel), xlab = "Log10[Size factor]", col = "grey80")
ls.zeisel <- colSums(counts(sce.zeisel))
plot(
    ls.zeisel,
    lib.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)

### ----- Exercise Day 3 -----------------------------------------------------------------------------------------------

# Are ls.zeisel and lib.sf.zeisel identical?
## No, they are not
identical(ls.zeisel, lib.sf.zeisel) # FALSE

# Are they proportional?
# a * b = c
# a * b / a = c / a
# a / a * b = c / a
# 1 * b = c / a
# b = c / a

prop1 <- lib.sf.zeisel/ls.zeisel
prop2 <- ls.zeisel/lib.sf.zeisel
table(prop1)
table(prop2)

# Compute lib.sf.zeisel manually
## First we need to compute the sums
zeisel_sums <- colSums(counts(sce.zeisel))
identical(zeisel_sums, ls.zeisel)

## Next, make them have unity mean
mean(zeisel_sums)
zeisel_size_factors <- zeisel_sums/mean(zeisel_sums)
identical(zeisel_size_factors, lib.sf.zeisel)



## ----all_code3, cache=TRUE, dependson='all_code2'--------------------------------------------------------------------
# Normalization by convolution

library('scran')
# Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)
# Compute deconvolution size factors
deconv.sf.zeisel <-
    calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)

table(clust.zeisel)
class(deconv.sf.zeisel)
length(deconv.sf.zeisel)


# Examine distribution of size factors
summary(deconv.sf.zeisel)

# Now we plot the distribution
pdf("deconv.sf.zeisel.pdf") # to be able to visualize the plot working from the cluster
hist(log10(deconv.sf.zeisel), xlab = "Log10[Size factor]",
     col = "grey80")
plot(
    ls.zeisel,
    deconv.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)
dev.off()


## ---------------- Exercise Day 4 ----------------------------------------------
# How many quick clusters did we get?
## We got 12
table(clust.zeisel)

# How many cells per quick cluster did we get?
## From 113 to 325 cells
summary(clust.zeisel)
## We can also get more information if we sort and make a vector
data.clust.zeisel <- table(clust.zeisel)
sort(data.clust.zeisel)
summary(as.vector(data.clust.zeisel))

# How many quick clusters will we get if we set the minimum size to 200? Use 100 as the seed.
## We got 10

set.seed(100)
clust.zeisel.2 <- quickCluster(sce.zeisel, min.size = 200)

deconv.sf.zeisel.2 <-
    calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)

table(clust.zeisel.2)

# How many lines do you see?
## Something like 2 lines


## ----all_code4, cache=TRUE, dependson='all_code3'--------------------------------------------------------------------
# Library size factors vs. convolution size factors

# Colouring points using the supplied cell-types
pdf("color.zeisel.pdf")
plot(
    lib.sf.zeisel,
    deconv.sf.zeisel,
    xlab = "Library size factor",
    ylab = "Deconvolution size factor",
    log = 'xy',
    pch = 16,
    col = as.integer(factor(sce.zeisel$level1class))
)
abline(a = 0, b = 1, col = "red")
dev.off()

# Cells are separating in groups by cell type-specific deviations
table(sce.zeisel$level1class)

# Scaling and log-transforming
# Using the convolution size factor we computed earlier

sce.zeisel <- logNormCounts(sce.zeisel)

# The above added a "Logcounts" assay to the SCE
assayNames(sce.zeisel)

### This divides the count for each feature by the appropiate size factor for that cell (normalized values)
### From this point on, we will use lognormcounts for the downstream analysis
## the log-transformation results in differences in the log-values representing the log-fold changes in expression


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()
