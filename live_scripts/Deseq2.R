library(tidyverse)
library(DESeq2)

load("Robjects/preprocessing.RData")

dim(countdata)
sampleinfo

design <- as.formula(~ Status)
modelMatrix <- model.matrix(design, data = sampleinfo)
modelMatrix

sampleinfo$Status <- factor(sampleinfo$Status, 
                            levels = c("virgin", "pregnant", "lactate"))
modelMatrix <- model.matrix(design, data = sampleinfo)
modelMatrix

# the deseq object

ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = sampleinfo,
                                     design = design)

# DESeq2 workflow

#size factors
ddsObj <- estimateSizeFactors(ddsObj.raw)

colData(ddsObj.raw)
colData(ddsObj)

logcounts <- log2(countdata + 1)
limma::plotMA(logcounts)
abline(h=0, col="red")

normalizedCounts <- counts(ddsObj, normalized = TRUE)
logNormalizedCounts <- log2(normalizedCounts + 1)

limma::plotMA(logNormalizedCounts)
abline(h=0, col="red")

# dispersion

ddsObj <- estimateDispersions(ddsObj)

plotDispEsts(ddsObj)

# GLM and Wald

ddsObj <- nbinomWaldTest(ddsObj)

# DESeq

ddsObj <- DESeq(ddsObj.raw)

# results

res <- results(ddsObj, alpha = 0.05)
res

modelMatrix
resultsNames(ddsObj)

resLvV_status <- res
rm(res)

resPvV_status <- results(ddsObj,
                         name = "Status_pregnant_vs_virgin",
                         alpha = 0.05)
resPvV_status

sum(resPvV_status$padj < 0.05)
sum(resPvV_status$padj < 0.05, na.rm = TRUE)

topGenesPvV <- as.data.frame(resPvV_status) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj) %>% 
  head(100)
topGenesPvV

# Exercise 1
#Challenge 1

design <- as.formula(~ CellType + Status)
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = sampleinfo,
                                     design = design)
ddsObj <- DESeq(ddsObj.raw)
resLvV <- results(ddsObj, alpha = 0.05)

resultsNames(ddsObj)

resLvB <- results(ddsObj, name = "CellType_luminal_vs_basal", alpha = 0.05)
resLvB

sum(resLvB$padj < 0.05 & resLvB$log2FoldChange > 0, na.rm = TRUE)

# Challenge 2

resPvL <- results(ddsObj,
                  alpha = 0.05,
                  contrast = c("Status", "pregnant", "lactate"))
resPvL

# Interaction models and the LRT test

vstcounts <- vst(ddsObj.raw, blind=TRUE)
plotPCA(vstcounts, intgroup=c("Status", "CellType"))

design.reduced <- as.formula(~ CellType)
ddsObj_C <- DESeq(ddsObj, test = "LRT", reduced = design.reduced)

resCvsCS <- results(ddsObj_C, alpha = 0.05)
resCvsCS
sum(resCvsCS$padj < 0.05, na.rm = TRUE)

# Exercise 2

design <- as.formula(~ CellType * Status)
ddsObj2.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                      colData = sampleinfo,
                                      design = design)
design.reduced <- as.formula(~ CellType + Status)

ddsObj2 <- DESeq(ddsObj2.raw, test = "LRT", reduced = design.reduced)

resCSvsCxS <- results(ddsObj2, alpha = 0.05)
table(resCSvsCxS$padj < 0.05)
