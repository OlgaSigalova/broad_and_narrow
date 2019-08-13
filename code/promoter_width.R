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

data("FANTOM5humanSamples")
data1 <- FANTOM5humanSamples
as.data.frame(data1)
tab = as.data.frame(data1)
tab

sample1 = "blood__adult__pool1" 
ce <- importPublicData(source = "FANTOM5", dataset = "human", sample = sample1)
ce
normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.25, T = 5 * 10 ^ 4)
ce
threshold <- 0.03
maxDist <- 20

clusterCTSS(object = ce
            , threshold = threshold
            , thresholdIsTpm = TRUE
            , nrPassThreshold = 1
            , removeSingletons = FALSE
            , method = "distclu"
            , maxDist = maxDist
            , useMulticore = TRUE)

cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = TRUE)
qup <- 0.9
qlow <- 0.1
quantilePositions(ce, clusters = "tagClusters", qLow = qlow, qUp = qup)

claster = tagClusters(ce, sample = sample1, returnInterquantileWidth = TRUE, qLow = qlow,
            qUp = qup)
claster

gr = GRanges(claster)

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts(txdb)

pro = promoters(txdb, upstream = 250, downstream = 250, 
                columns = c("TXID", "GENEID"))

res = findOverlaps(gr, pro, minoverlap = 1)

df1 = data.frame(gr[queryHits(res)])
df2 = data.frame(pro[subjectHits(res)])
df = cbind(df1, df2)

df %<>%
  filter(tpm>0.3, !is.na(as.numeric(df$GENEID))) 
 df1 = df %>%
  select("seqnames", "start", "end",                
         "strand","interquantile_width", "tpm" ) %>%
   unique()

 ggplot(df1, aes(x =interquantile_width)) +
  geom_histogram(bins = 300, alpha = 0.5)
 












