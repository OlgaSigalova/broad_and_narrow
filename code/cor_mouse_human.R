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

#cor_mouse
table <- read.csv("/home/tmpam10/Desktop/broad_and_narrow/promoter_widths_mouse_tissues.csv", header = TRUE)

df=table %>%
  group_by(GENEID, sample) %>%
  mutate(tmpmax = max(tpm)) %>%
  filter(tmpmax == tpm)

df1=df %>%
  mutate(iwmax = max(interquantile_width)) %>%
  filter(iwmax == interquantile_width)

df2 = df1 %>%
  group_by(GENEID) %>%
  summarize(n = n(),
            mean = mean(interquantile_width))

ggplot(df2, aes(x = factor(n), y = mean)) +
  geom_violin() +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  xlab("number of tissues") +
  ylab("mean promoter width")

#cor_human
table <- read.csv("/home/tmpam10/Desktop/broad_and_narrow/promoter_widths_human_tissues.csv", header = TRUE)

df=table %>%
  group_by(GENEID, sample) %>%
  mutate(tmpmax = max(tpm)) %>%
  filter(tmpmax == tpm)

df1=df %>%
  mutate(iwmax = max(interquantile_width)) %>%
  filter(iwmax == interquantile_width)

df2 = df1 %>%
  group_by(GENEID) %>%
  summarize(n = n(),
            mean = mean(interquantile_width))

ggplot(df2 %>% filter(n<6), aes(x = factor(n), y = mean)) +
  geom_violin() +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  xlab("number of tissues") +
  ylab("mean promoter width")




  
  
  
  