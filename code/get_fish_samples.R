library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(CAGEr)
library(GenomicRanges)
library(plyr)
library(dplyr)
library(entropy)
library(mixtools)
library(ggplot2)
library(ggpubr)
library(SDMTools)

location = "/Users/sigalova/Desktop/SMTB/broad_and_narrow/data/ZF/cage_RData/ZebrafishCAGE.RData"
load(location)

temp = as.data.frame(ZebrafishCAGE[["development"]])
ce = as(temp, "CAGEset")

librarySizes(ce)
print(ce)

# combining all samples
mergeSamples(ce, rep(1, 12), "merged_samples")

# parameters for normalisation
plotReverseCumulatives(ce, fitInRange = c(10,10000), values = "raw", onePlot = TRUE)


# load genome
location = "/Users/sigalova/Desktop/SMTB/broad_and_narrow/data/ZF/model/danRer7_ensembl.dms"
txdb = makeTxDbFromGFF(file = location, format = "gtf")
