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

# subset samples
temp = temp[ , c(1:3, 10:11, 13:15)]

ce = as(temp, "CAGEset")
labels = sampleLabels(ce)
mergeSamples(ce, 1:5, labels)

librarySizes(ce)
print(ce)

# # combining all samples
# mergeSamples(ce, rep(1, 12), "merged_samples")

# parameters for normalisation

plotReverseCumulatives(ce, fitInRange = c(10,10000), values = "raw", onePlot = TRUE)

normalizeTagCount(ce, method = "powerLaw", fitInRange = c(10, 10000), 
                  alpha = 1.19, T = 10 ^ 6)


threshold <- 10
maxDist <- 20

clusterCTSS(object = ce, 
            threshold = threshold, 
            thresholdIsTpm = FALSE, 
            nrPassThreshold = 1, 
            removeSingletons = FALSE, 
            method = "distclu", 
            maxDist = maxDist, useMulticore = TRUE)

cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = FALSE)
qup <- 0.9
qlow <- 0.1
quantilePositions(ce, clusters = "tagClusters", qLow = qlow, qUp = qup)


# load genome
location = "/Users/sigalova/Desktop/SMTB/broad_and_narrow/data/ZF/model/danRer7_ensembl.dms"
txdb = makeTxDbFromGFF(file = location, format = "gtf")

prom <- promoters(txdb, upstream = 250, downstream = 250, columns = c("TXID", "GENEID"))
min <- 1

# begin work with one sample

myfun = function(prom, sampl, ce, qup = 0.9, qlow = 0.1, min = 1) {
  clusters <- tagClustersGR(ce, sample = sampl, 
                            returnInterquantileWidth = TRUE, 
                            qLow = qlow, qUp = qup)
  gr = GRanges(clusters)
  
  
  res = findOverlaps(gr, prom, minoverlap = min)
  
  df1 = data.frame(gr[queryHits(res)])
  df2 = data.frame(prom[subjectHits(res)])
  df_res = cbind(df1, df2)
  
  
  names(df_res)[15:19] = paste(names(df_res)[15:19], "gene", sep = "_")
  
  df_res %>%
    filter(tpm > 10, nchar(GENEID) > 0) %>%
    select(interquantile_width, cluster, seqnames, start, end, GENEID, tpm) %>%
    unique()
  
}


endtable <- data.frame(matrix(nrow = 0, ncol = 8))
names(endtable) = c("interquantile_width", "cluster", "seqnames", 
                    "start", "end", "GENEID", "tpm", "sample")

for(sample in labels) {
  a <- myfun(prom, sampl = sample, ce)
  a$sample <- sample
  endtable = rbind(endtable, a)
}

write.csv(endtable, file = "/Users/sigalova/Desktop/SMTB/broad_and_narrow/promoter_widths_fish_embryo_stages.csv", row.names = FALSE)


