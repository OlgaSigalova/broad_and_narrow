library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(CAGEr)
library(GenomicRanges)
library(plyr)
library(dplyr)
library(entropy)
library(mixtools)
library(ggplot2)
library(ggpubr)
library(SDMTools)


data("FANTOM5mouseSamples")
data1 <- FANTOM5mouseSamples
r=as.data.frame(data1)

r1 = r %>% filter(type == "tissue", grepl("embryo", sample), grepl("whole body", description), 
                  !grepl("internal control", description))

samples = r1 %>%
  select(sample) %>% unlist(use.names = FALSE) 
 
samples = samples[c(1, 3, 4, 6, 9)]

ce <- importPublicData(source = "FANTOM5", dataset = "mouse", sample = samples)


plotReverseCumulatives(ce, fitInRange = c(10,10000), values = "raw", onePlot = TRUE)

normalizeTagCount(ce, method = "powerLaw", fitInRange = c(10, 10000), 
                  alpha = 1.16, T = 10 ^ 6)


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

# mouse genes

txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
columns(txdb)

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
      filter(tpm > 10, !is.na(as.numeric(GENEID))) %>%
      select(interquantile_width, cluster, seqnames, start, end, GENEID, tpm) %>%
      unique()
    
}


endtable <- data.frame(matrix(nrow = 0, ncol = 8))
names(endtable) = c("interquantile_width", "cluster", "seqnames", 
                    "start", "end", "GENEID", "tpm", "sample")

for(sample in samples) {
  a <- myfun(prom, sampl = sample, ce)
  a$sample <- sample
  endtable = rbind(endtable, a)
}

write.csv(endtable, file = "/Users/sigalova/Desktop/SMTB/broad_and_narrow/promoter_widths_mouse_embryo_embryo_stages.csv", row.names = FALSE)


