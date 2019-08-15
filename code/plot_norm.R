
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

tab2 = tab %>%
  filter(type == "tissue", grepl("adult.*pool1", description))
tab2
names(tab2)

samples = tab2$sample
ce <- importPublicData(source = "FANTOM5", dataset = "human", sample = samples[c(1, 3, 6, 8, 10)])
ce

plotReverseCumulatives(ce, fitInRange = c(10,10000), values = "raw", onePlot = TRUE)


data("FANTOM5mouseSamples")
data1 <- FANTOM5mouseSamples
r=as.data.frame(data1)

r1 = r %>% filter(type == "tissue")


samples = r %>% filter(type == "tissue") %>% select(sample) %>% unlist(use.names = FALSE)
samples = samples[c(10,16,23,132)]


ce <- importPublicData(source = "FANTOM5", dataset = "mouse", sample = samples)

plotReverseCumulatives(ce, fitInRange = c(10,10000), values = "raw", onePlot = TRUE)



normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), 
                  alpha = 1.18, T = 10 ^ 7)



  