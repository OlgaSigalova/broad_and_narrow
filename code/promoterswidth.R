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
r=as.data.frame(data1)

sample1 <- "heart__adult__pool1"
#sample2 <- "blood__adult__pool1"
#sample3 <- "brain__adult__pool1"

#samples <- c(sample1, sample2, sample3)
ce <- importPublicData(source = "FANTOM5", dataset = "human", sample = sample1)

normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), 
                  alpha = 1.25, T = 5 * 10 ^ 4)


threshold <- 0.03
maxDist <- 20

clusterCTSS(object = ce, threshold = threshold, thresholdIsTpm = TRUE, 
            nrPassThreshold = 1, removeSingletons = FALSE, method = "distclu", 
            maxDist = maxDist, useMulticore = TRUE)

cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = TRUE)
qup <- 0.9
qlow <- 0.1
quantilePositions(ce, clusters = "tagClusters", qLow = qlow, qUp = qup)
clusters <- tagClustersGR(ce, sample = sample1, 
              returnInterquantileWidth = TRUE, 
              qLow = qlow, qUp = qup)
gr = GRanges(clusters)


# human genes
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts(txdb)
prom <- promoters(txdb, upstream = 250, downstream = 250, columns = c("TXID", "GENEID"))
min <- 1
res = findOverlaps(gr, prom, minoverlap = min)

df1 = data.frame(gr[queryHits(res)])
df2 = data.frame(prom[subjectHits(res)])
df_heart = cbind(df1, df2)

names(df_heart)[15:19] = paste(names(df_heart)[15:19], "gene", sep = "_")

df_heart %<>%
  filter(tpm>0.3, !is.na(as.numeric(GENEID)))

tab_heart = df_heart %>% 
  select(interquantile_width, cluster, seqnames, start, end, GENEID, tpm) %>%
  unique()


ggplot(tab, aes(x = interquantile_width)) +
  geom_histogram(bins = 300)







