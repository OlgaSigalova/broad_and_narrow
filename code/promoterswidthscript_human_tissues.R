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

samples = r %>% filter(type == "tissue" & grepl("pool1", sample), grepl("adult", sample)) %>% select(sample) %>% unlist(use.names = FALSE)
samples = samples[c(1, 3, 6, 8, 10)]


ce <- importPublicData(source = "FANTOM5", dataset = "human", sample = samples)

normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), 
                  alpha = 1.18, T = 10 ^ 7)


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

# human genes

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts(txdb)
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

    df_res %<>%
      filter(tpm > 10, !is.na(as.numeric(GENEID))) %>%
      select(interquantile_width, cluster, seqnames, start, end, GENEID, tpm) %>%
      unique()
    
}





# ggplot(tab, aes(x = interquantile_width)) +
#   geom_histogram(bins = 300)

endtable <- data.frame(matrix(nrow = 0, ncol = 8))
names(endtable) = c("interquantile_width", "cluster", "seqnames", 
                    "start", "end", "GENEID", "tpm", "sample")

for(sample in samples) {
  a <- myfun(prom, sampl = sample, ce)
  a$sample <- sample
  endtable = rbind(endtable, a)
}

write.csv(endtable, file = "/Users/sigalova/Desktop/SMTB/broad_and_narrow/promoter_widths_human_tissues.csv", row.names = FALSE)


