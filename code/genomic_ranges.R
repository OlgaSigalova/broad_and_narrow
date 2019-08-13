library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# create object: genomic ranges
genes = GRanges( 
  seqnames = c ("chr1", "chr34", "chr55"),
  IRanges (start = c(10, 15, 20), end =c (40, 45, 50)),
  strand = c( "-", "+", "*")
  )

genes

# resize and get promoters
promotres = resize (genes, 10)
promotres
promoters(genes, upstream = 10, downstream = 10)

# create second object with genomic ranges
cage= GRanges( 
  seqnames = c ("chr1", "chr34", "chr55"),
  IRanges (start = c(10, 15, 120), end =c (25, 35, 150)),
  strand = c( "-", "+", "*")
)

# overlap two sets of intervals
res = findOverlaps(genes, cage, minoverlap = 3)
sumres = c(queryHits(res), subjectHits(res))
sumres
subset.genes[2]

df1 = as.data.frame(genes[queryHits(res)])
df2 = as.data.frame(cage[subjectHits(res)])
df = cbind(df1, df2)

# format names
names(df)[1:5] = paste(names(df)[1:5], "gene", sep = "_")
names(df)[6:10] = paste(names(df)[6:10], "cage", sep = "_")

# human genes
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts(txdb)

# promoters for transcripts
promoters(transcripts(txdb), upstream = 250, downstream = 250)
# ..or simply
promoters(txdb, upstream = 250, downstream = 250)
