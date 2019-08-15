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
library(tidyverse)
library(magrittr)

table <- read.csv("Desktop/broad_and_narrow/smtb19/code/endtable.csv", header = TRUE)

df=table %>%
  group_by(GENEID, sample) %>%
  mutate(tmpmax = max(tpm)) %>%
  filter(tmpmax == tpm)

df1=df %>%
  mutate(iwmax = max(interquantile_width)) %>%
  filter(iwmax == interquantile_width)

newtable=df1 %>%
  select(GENEID, sample, interquantile_width) %>%
  unique() %>%
  spread(key = sample, value = interquantile_width)

tablewhoutNA=newtable %>%
  na.omit()


ggplot(df1 %>% filter(interquantile_width<100), aes(x = interquantile_width)) +
  geom_histogram(bins = 150) + 
  facet_wrap(~sample, ncol = 3)

ggplot(tablewhoutNA, aes(x = aorta__adult__pool1, y = brain__adult__pool1)) +
  geom_point() +
  xlim(0, 100) + ylim(0, 100) +
  geom_smooth(method = "lm")

newtable %>%
  group_by(!is.na(aorta__adult__pool1), !is.na(brain__adult__pool1)) %>%
  summarize(n(), mean(aorta__adult__pool1), mean(brain__adult__pool1))

round(cor(tablewhoutNA[, 2:6]), 2)
