#List files in directory
sorghumFiles <- list.files("read_counts")
#Create table from files names
countsTable <- data.frame(read_counts = sub(".tsv","",sorghumFiles),
                          fileName = sorghumFiles,
                          treatment = sub("-..tsv","",sorghumFiles))
#set the treatment column to factor
countsTable$treatment <- factor(countsTable$treatment)
#Load DESeq2
library(DESeq2)
#Create the DESeq2 object from the HTSeqCount tables
dds <- DESeqDataSetFromHTSeqCount(sampleTable = countsTable,
                                  directory = "read_counts",
                                  design = ~ treatment)
#Create factor levels so well_watered is considered the reference treatment
dds$treatment <- relevel(dds$treatment, ref = "well_watered")
#show the dds object
dds
#Prefilter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
#estimate size factor - library size
dds <- estimateSizeFactors(dds)
#show library size asjustment
dds$sizeFactor
#estimate dispersion
dds <- estimateDispersions(dds)
#plot the dispersions
plotDispEsts(dds)
#variance stabilization
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
#calculate distance between samples
sampleDistances <- dist(t(assay(vsd)))
#plot pheatmap
install.packages("pheatmap")
install.packages("RColorBrewer")
library(pheatmap)
library(RColorBrewer)
sampledistancematrix <- as.matrix(sampleDistances)
rownames(sampledistancematrix) <- paste(colnames(vsd), vsd$type)
colnames(sampledistancematrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues"))) (255)
pheatmap(sampledistancematrix,
         clustering_distance_rows = sampleDistances,
         clustering_distance_cols = sampleDistances,
         col=colors)
#PCA
install.packages("digest")
plotPCA(vsd, intgroup = c("treatment"))
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds)
#call differentially expressed genes with Wald Test
dds <- nbinomWaldTest(dds)
#extract results for pair-wise comparison
resultsTable <- as.data.frame(results(dds, contrast = c("treatment", "well_watered","drought")))
#show how many genes padj < 0.05
nrow(na.omit(resultsTable[resultsTable$padj <0.05,]))
#volcano plot
plotMA(results(dds, contrast = c("treatment","well_watered","drought")), ylim=c(-2,2))
