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

table_compr <- read.csv("/home/tmpam10/Desktop/broad_and_narrow/smtb19/databank/DRSC_human_mouse.csv", header = TRUE)

compres = table_compr %>%
  filter(Rank == "high") %>%
  select(Human.GeneID, Mouse.GeneID)

result_table_human = merge(compres, newtable, by.x = "Human.GeneID", by.y = "GENEID")
result_table = merge(result_table_human, newtable_mouse, by.x = "Mouse.GeneID", by.y = "GENEID")

round(cor(result_table[, 3:11], use = "complete.obs"), 2)


mouse_brain = ifelse(result_table$cerebellum__adult < 11, "narrow", "broad")

human_brain = ifelse(result_table$brain__adult__pool1 < 12, "narrow", "broad")



brain = result_table %>%
  cbind.data.frame(mouse_brain) %>%
  cbind.data.frame(human_brain) 

brain_t = data.frame(table(mouse_brain, human_brain))


mouse_aorta = ifelse(result_table$aorta__adult < 11, "narrow", "broad")

human_aorta = ifelse(result_table$aorta__adult__pool1 < 12, "narrow", "broad")



aorta = result_table %>%
  cbind.data.frame(mouse_aorta) %>%
  cbind.data.frame(human_aorta) 

aorta_t = data.frame(table(mouse_aorta, human_aorta))

