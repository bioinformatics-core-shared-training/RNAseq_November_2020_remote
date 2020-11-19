# data set:
# mammary epithelium in mouse
# two types of cells: basal and luminal
# three types of mice
# 6 sample groups
# 2 replicates per sample group
# 12 samples in total

# load packages
#--------------------

library(DESeq2)
library(tidyverse)

# load sample meta data
#--------------------

sampleinfo <- read_tsv("data/SampleInfo.txt")
sampleinfo

# load count matrix
#--------------------

seqdata <- read_tsv("data/GSE60450_Lactation.featureCounts", comment='#')
seqdata

# format the data
#--------------------

countdata <- seqdata %>%
  column_to_rownames('Geneid') %>%
  rename_all(str_remove, '.bam') %>%
  select(sampleinfo$Sample) %>%
  as.matrix()
class(countdata)
dim(countdata)  
head(countdata)  

# filter lowly expressed genes
#--------------------

keep <- rowSums(countdata) > 5
head(keep)
table(keep)

countdata <- countdata[keep,]
dim(countdata)
head(countdata)

#--------------------
# QC - count distributions and visualisation
#--------------------

# raw counts
#--------------------

hist(countdata[,1])
summary(countdata)

boxplot(countdata, las=2)

# variance-mean relationship

plot(x=rowMeans(countdata), y=rowSds(countdata),
     main="sd vs mean",
     xlim=c(0,10000),
     ylim=c(0,5000))
abline(lm(rowSds(countdata) ~ rowMeans(countdata)),
       col="green")

# transformations
#--------------------

# log2 
# Variance Stabilizing Transformation VST
# regularized log rlog

# log2 transformation
#--------------------

logcounts  <- log2(countdata + 1)
head(logcounts)

# distribution

statusCol <- match(sampleinfo$Status,
                   c("virgin", "pregnant", "lactate")) + 1

boxplot(logcounts,
        las=2,
        ylab="log2 counts",
        col=statusCol,
        main="log2 counts"
        )
abline(h=median(logcounts), col="blue")

# log2 sd vs mean

log2Means <- rowMeans(logcounts)
log2Sds <- rowSds(logcounts)
plot(x=log2Means, y=log2Sds,
     main="log2 sd vs mean")
abline(lm(log2Sds ~ log2Means), col="green")

# VST
#--------------------

vst_counts <- vst(countdata)

# distribution

statusCol <- match(sampleinfo$Status,
                   c("virgin", "pregnant", "lactate")) + 1

boxplot(vst_counts,
        las=2,
        ylab="vst counts",
        col=statusCol,
        main="vst counts"
)
abline(h=median(vst_counts), col="blue")

# sd vs mean

vstMeans <- rowMeans(vst_counts)
vstSds <- rowSds(vst_counts)
plot(x=vstMeans, y=vstSds,
     main="vst sd vs mean")
abline(lm(vstSds ~ vstMeans), col="green")

# Challenge 1:
# use DESeq2 rlog function to transform
# plot distributions
# plot sd vs mean

rlogcounts <- rlog(countdata)
head(rlogcounts)

# distribution
boxplot(rlogcounts,
        las=2,
        ylab="rlog counts",
        col=statusCol,
        main="rlog counts"
)
abline(h=median(rlogcounts), col="blue")

# sd vs mean

rlogMeans <- rowMeans(rlogcounts)
rlogSds <- rowSds(rlogcounts)
plot(x=rlogMeans, y=rlogSds,
     main="vst sd vs mean")
abline(lm(rlogSds ~ rlogMeans), col="green")

# gene expression profiles - PCA
#--------------------

# prcomp()

# rlog'ed counts
rlogcounts <- rlog(countdata)

# pca: compute
pcDat <- prcomp(t(rlogcounts))
# ?prcomp

# pca: plot
library(ggfortify)

# basic
autoplot(pcDat)

# with samples color- and shape-coded by cell type and status
autoplot(pcDat,
         data=sampleinfo,
         colour="CellType",
         shape="Status",
         size=5)

# add sample name as label
library(ggrepel)

autoplot(pcDat,
         data=sampleinfo,
         colour="CellType",
         shape="Status",
         size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample),
                  box.padding = 0.8)

# fix sample name swap
sampleinfo.orig <- sampleinfo
sampleinfo <- sampleinfo %>%
  mutate(CellType = ifelse(Sample=="MCL1.DG", "basal", CellType)) %>%
  mutate(CellType = ifelse(Sample=="MCL1.LA", "luminal", CellType))
sampleinfo

# write updated table to file
dir.create("results")
write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")

# plot again to check sample name swap is fixed
autoplot(pcDat,
         data=sampleinfo,
         colour="CellType",
         shape="Status",
         size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample),
                  box.padding = 0.8)

# ?save.image
# save R objects
save(countdata, sampleinfo, file="results/preprocessing.RData")

