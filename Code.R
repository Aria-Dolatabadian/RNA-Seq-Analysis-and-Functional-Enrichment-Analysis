library(pheatmap)
library(ggfortify)
library(stats)
library(ggplot2)
library(corrplot)
library(DESeq2)
library(EDASeq)
library(gProfileR)
library(knitr)
library(gage)

#place the tsv files in the WD

counts_file <- "SRP029880.raw_counts.tsv"
coldata_file <- "SRP029880.colData.tsv"

counts <- as.matrix(read.table(counts_file, header = TRUE, sep = "\t"))
colData <- read.table(coldata_file, header = TRUE, sep = "\t")


summary(counts[,1:3])

cpm <- apply(subset(counts, select = c(-width)), 2, 
             function(x) x/sum(as.numeric(x)) * 10^6)

colSums(cpm)

# create a vector of gene lengths 
geneLengths <- as.vector(subset(counts, select = c(width)))

# compute rpkm 
rpkm <- apply(X = subset(counts, select = c(-width)),
              MARGIN = 2, 
              FUN = function(x) {
                10^9 * x / geneLengths / sum(as.numeric(x))
               })

colSums(rpkm)

#find gene length normalized values 
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

colSums(tpm)

#compute the variance of each gene across samples
V <- apply(tpm, 1, var)
#sort the results by variance in decreasing order 
#and select the top 100 genes 
selectedGenes <- names(V[order(V, decreasing = T)][1:100])

library(pheatmap)
pheatmap(tpm[selectedGenes,], scale = 'row', show_rownames = FALSE)


colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
pheatmap(tpm[selectedGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData)


library(ggfortify)

library(stats)
library(ggplot2)
#transpose the matrix 
M <- t(tpm[selectedGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)

#plot PCA results making use of ggplot2's autoplot function
#ggfortify is needed to let ggplot2 know about PCA data structure. 
autoplot(pcaResults, data = colData, colour = 'group')



summary(pcaResults)

library(stats)
correlationMatrix <- cor(tpm)

library(corrplot)
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.7) 

library(pheatmap)
# split the clusters into two based on the clustering similarity 
pheatmap(correlationMatrix,  
         annotation_col = colData, 
         cutree_cols = 2)

#remove the 'width' column
countData <- as.matrix(subset(counts, select = c(-width)))
#define the experimental setup 
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
#define the design formula
designFormula <- "~ group"

library(DESeq2)
library(stats)
#create a DESeq dataset object from the count matrix and the colData 
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = as.formula(designFormula))
#print dds object to see the contents
print(dds)

#For each gene, we count the total number of reads for that gene in all samples 
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

dds <- DESeq(dds)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]

#shows a summary of the results
print(DEresults)

library(DESeq2)
DESeq2::plotMA(object = dds, ylim = c(-5, 5))

library(ggplot2)
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
  geom_histogram(bins = 100)

library(DESeq2)
# extract normalized counts from the DESeqDataSet object
countsNormalized <- DESeq2::counts(dds, normalized = TRUE)

# select top 500 most variable genes
selectedGenes <- names(sort(apply(countsNormalized, 1, var), 
                            decreasing = TRUE)[1:500])

plotPCA(countsNormalized[selectedGenes,], 
        col = as.numeric(colData$group), adj = 0.5, 
        xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.6))

rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
  ylim(-50, 50) + theme_bw()

library(EDASeq)
par(mfrow = c(1, 2))
plotRLE(countData, outline=FALSE, ylim=c(-4, 4), 
        col=as.numeric(colData$group), 
        main = 'Raw Counts')
plotRLE(DESeq2::counts(dds, normalized = TRUE), 
        outline=FALSE, ylim=c(-4, 4), 
        col = as.numeric(colData$group), 
        main = 'Normalized Counts')

library(DESeq2)
library(gProfileR)
library(knitr)
# extract differential expression results
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))

#remove genes with NA values 
DE <- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

#get the list of genes of interest
genesOfInterest <- rownames(DE)

#calculate enriched GO terms
goResults <- gprofiler(query = genesOfInterest, 
                     organism = 'hsapiens', 
                     src_filter = 'GO', 
                     hier_filtering = 'moderate')



write.csv(goResults, "GO_results.csv", row.names = TRUE)




#Let's define the first gene set as the list of genes from one of the
#significant GO terms found in the GO analysis. order go results by pvalue
goResults <- goResults[order(goResults$p.value),]
#restrict the terms that have at most 100 genes overlapping with the query
go <- goResults[goResults$overlap.size < 100,]
# use the top term from this table to create a gene set 
geneSet1 <- unlist(strsplit(go[1,]$intersection, ','))

#Define another gene set by just randomly selecting 25 genes from the counts
#table get normalized counts from DESeq2 results
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)
geneSet2 <- sample(rownames(normalizedCounts), 25)

geneSets <- list('top_GO_term' = geneSet1,
                 'random_set' = geneSet2)

library(gage)
#use the normalized counts to carry out a GSEA. 
gseaResults <- gage(exprs = log2(normalizedCounts+1), 
           ref = match(rownames(colData[colData$group == 'CTRL',]), 
                       colnames(normalizedCounts)), 
           samp = match(rownames(colData[colData$group == 'CASE',]), 
                        colnames(normalizedCounts)),
           gsets = geneSets, compare = 'as.group')


write.csv(gseaResults$greater, "greater_results.csv", row.names = TRUE)
write.csv(gseaResults$greater, "less_results.csv", row.names = TRUE)

library(pheatmap)
# get the expression data for the gene set of interest
M <- normalizedCounts[rownames(normalizedCounts) %in% geneSet1, ]
# log transform the counts for visualization scaling by row helps visualizing
# relative change of expression of a gene in multiple conditions
pheatmap(log2(M+1), 
         annotation_col = colData, 
         show_rownames = TRUE, 
         fontsize_row = 8,
         scale = 'row', 
         cutree_cols = 2, 
         cutree_rows = 2)

#https://compgenomr.github.io/book/gene-expression-analysis-using-high-throughput-sequencing-technologies.html#tab:GOanalysistable
