library(GO.db)
library(topGO)

#Load DESeq2 results table
res <- read.table("DESeq2_results.tsv",header=TRUE,sep="\t")
res$id <- sub('.v3.2','',res$id)

#Read in topGO terms
goTerms <- readMappings(file="Sorghum_topGO.tsv")
table(res$sig)

#format higher expressed genes
upGenes <- factor(as.integer(res$id %in% res[res$padj < 0.05 & res$log2FC >= 1,]$id))
names(upGenes) <- res$id
table(upGenes)

#GO term enrichment for BP
upBP <- new("topGOdata",description="Biological Process",ontology="BP",
          allGenes=upGenes,annot=annFUN.gene2GO,nodeSize=3,gene2GO=goTerms)
upFisherBP <- runTest(upBP,algorithm = "parentchild",statistic = "fisher")
upBPgenTable <- GenTable(upBP,Fisher=upFisherBP,ranksOf="Fisher",
                         topNodes=length(score(upFisherBP)))
head(upBPgenTable)
write.table(upBPgenTable, file='higher_expressed_BP.tsv', quote=FALSE, sep='\t')

#GO term enrichment for MF
upMF <- new("topGOdata",description="Molecular Function",ontology="MF",
            allGenes=upGenes,annot=annFUN.gene2GO,nodeSize=3,gene2GO=goTerms)
upFisherMF <- runTest(upMF,algorithm = "parentchild",statistic = "fisher")
upMFgenTable <- GenTable(upMF,Fisher=upFisherMF,ranksOf="Fisher",
                         topNodes=length(score(upFisherMF)))
head(upMFgenTable)
write.table(upMFgenTable, file='higher_expressed_MF.tsv', quote=FALSE, sep='\t')

#GO term enrichment for CC
upCC <- new("topGOdata",description="Cellular Component",ontology="CC",
            allGenes=upGenes,annot=annFUN.gene2GO,nodeSize=3,gene2GO=goTerms)
upFisherCC <- runTest(upCC,algorithm = "parentchild",statistic = "fisher")
upCCgenTable <- GenTable(upCC,Fisher=upFisherCC,ranksOf="Fisher",
                         topNodes=length(score(upFisherCC)))
head(upCCgenTable)
write.table(upCCgenTable, file='higher_expressed_CC.tsv', quote=FALSE, sep='\t')

#format lower expressed genes
downGenes <- factor(as.integer(res$id %in% res[res$padj > 0.05 & res$log2FC <= -1,]$id))
names(downGenes) <- res$id
table(downGenes)

#GO term enrichment for BP
downBP <- new("topGOdata",description="Biological Process",ontology="BP",
            allGenes=downGenes,annot=annFUN.gene2GO,nodeSize=3,gene2GO=goTerms)
downFisherBP <- runTest(downBP,algorithm = "parentchild",statistic = "fisher")
downBPgenTable <- GenTable(downBP,Fisher=downFisherBP,ranksOf="Fisher",
                         topNodes=length(score(downFisherBP)))
head(downBPgenTable)
write.table(downBPgenTable, file='lower_expressed_BP.tsv', quote=FALSE, sep='\t')

#GO term enrichment for MF
downMF <- new("topGOdata",description="Molecular Function",ontology="MF",
            allGenes=downGenes,annot=annFUN.gene2GO,nodeSize=3,gene2GO=goTerms)
downFisherMF <- runTest(downMF,algorithm = "parentchild",statistic = "fisher")
downMFgenTable <- GenTable(downMF,Fisher=downFisherMF,ranksOf="Fisher",
                         topNodes=length(score(downFisherMF)))
head(downMFgenTable)
write.table(downMFgenTable, file='lower_expressed_MF.tsv', quote=FALSE, sep='\t')

#GO term enrichment for CC
downCC <- new("topGOdata",description="Cellular Component",ontology="CC",
            allGenes=downGenes,annot=annFUN.gene2GO,nodeSize=3,gene2GO=goTerms)
downFisherCC <- runTest(downCC,algorithm = "parentchild",statistic = "fisher")
downCCgenTable <- GenTable(downCC,Fisher=downFisherCC,ranksOf="Fisher",
                         topNodes=length(score(downFisherCC)))
head(downCCgenTable)
write.table(downCCgenTable, file='lower_expressed_CC.tsv', quote=FALSE, sep='\t')





