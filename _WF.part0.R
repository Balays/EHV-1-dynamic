


#### Settings and annotation


##
##### Import libraries and functions
library(rtracklayer)
library(Rsamtools)
library(ggsci)
library(seqinr)
library(Rsamtools)
library(data.table)
library(stringi)
library(ggplot2)
library(tidygenomics)
library(tidyr)
library(dplyr)
library(fuzzyjoin)
#library(future.apply)
library(purrr)
library(BiocParallel)

## Own functions -->> Needs to be present
misc.dir    <- 'functions'

if (.Platform$OS.type!="windows") {

  misc.dir  <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', misc.dir))),
                      stri_replace_first_regex(misc.dir, '.*:\\/', '')
  )

}

for(f in list.files(misc.dir,    '*.R', full.names = T)) { try({source(f)}) }

bam.flags           <- fread(paste0(misc.dir, '/bam.flags.tsv'))
gff_compare_classes <- fread(paste0(misc.dir, '/gff_compare.txt'))
nproc <- 48

#### ####


##
##### Metadata ####
metadata     <- read.delim('EHV-1_metadata.tsv')
metafilt     <- metadata[,c("sample_id", "rep", "time", "time_num", "cell_line")]
colnames(metafilt) <- c('sample', 'rep', 'hpi', 'Time', 'cell_line')
metafilt$hpi <- factor(metafilt$hpi, levels=c('1h', '2h' ,'4h', '6h', '8h', '12h', '18h', '24h', '48h', 'mix'))
metafilt$group <- paste0(metafilt$cell_line, '_', metafilt$hpi)

metafilt$cell_line <- factor(metafilt$cell_line, levels=c("C6", "PK-15", "PC-12", 'RK-13'))
metafilt <- metafilt[order(metafilt$cell_line, metafilt$hpi), ]

## which columns of the metadata table to be included in the downstream analysis?
## The first element of this vector should be sample ID, which should be the same as the .bam file's names.
metacols <- c('sample', 'hpi', 'Time', 'cell_line', 'group')
#### ####
##


##
##### Genome and annotation ####

### Reference genome
virus  <- 'EHV-1'
genome <- 'NC_001491.2'
#fasta.ref  <- if (.Platform$OS.type!="windows") {'/mnt/e/data/genomes/Rnor_6.0.104_and_LT934125.1.fasta'} else {'E:/data/genomes/Rnor_6.0.104_and_LT934125.1.fasta'}
fasta.ref <- paste0(genome, '.fasta')
fasta  <- seqinr::read.fasta(fasta.ref)
l_genome <- unlist(lapply(fasta, length))

### Annotation
## TSS clusters (manual annotation)
#tss.clusters <- read.delim('PRV.clusters.txt') ## not used

## CDSs
create.ann.from.gff <- T
gff.file <- paste0(genome, '.gff3')

if (create.ann.from.gff) {
  gff           <- data.table(as.data.frame(rtracklayer::import.gff(gff.file)))
  CDS.df        <- gff[gff$type == 'CDS', c("seqnames", "start", "end", "strand", "type", "product", 'ID', 'Name', 'gene')]
  CDS.df$Name   <- CDS.df$gene

  ## differentiate multicopy genes
  CDS.df[,copy_number := 1:n_distinct(strand), by=.(gene)]
  CDS.df[copy_number > 1, gene := paste0(gene, '_', copy_number)]

  CDS.df[, part := 1:.N, by=.(gene)]
  CDS.df[, part := 1:.N, by=.(gene)]

  CDS.df[part != 1, ID   := paste(gene, part, sep='_')]

  feature.df <- CDS.df[,c("seqnames", "start", "end", "strand", "type", "gene", 'part', 'ID')]

  write.table(feature.df, 'feature.df.tsv', sep = '\t', row.names = F, quote = F)
} else {
  feature.df <- read.delim('feature.df.tsv')
  feature.colname <- 'ID'
  stopifnot(nrow(feature.df) == luniq(feature.df[,feature.colname]))

}

feature.colname <- 'gene'
by <- c('seqnames', 'start', 'end')

## Import gene clusters

gene.clusters.all <- fread('EHV-1.genes.TSS.TES.txt')
#gene.clusters.all[,TES.canonic := fifelse(TES.canonic == cluster_TES, )]
gene.clusters.all <- gene.clusters.all[order(CDS.start),]
ORFs  <- unique(gene.clusters.all$gene)
genes <- unique(gene.clusters.all$gene)

## add unknown to NA kinetic class
gene.clusters.all[,Kinetic_class := fifelse(is.na(Kinetic_class), 'unknown', Kinetic_class)]

## add non-coding genes
#
add.nc.genes <- T
if(add.nc.genes) {
  feature.nc <- as.data.frame(gene.clusters.all[type=='NC-gene',.(seqnames, start=CDS.start, end=CDS.end, strand, type, gene, part=1, ID=gene)])
  feature.nc$gene_name <- feature.nc$gene
  #feature.nc$gene_name[feature.nc$gene == 'NOIR']  <- 'NOIR1'
  #feature.nc$gene_name[feature.nc$gene == 'NOIR2'] <- 'NOIR1'
  #feature.nc$gene_name[feature.nc$gene == 'NOIR2-2'] <- 'NOIR2'
  #feature.nc$gene_name[feature.nc$gene == 'NOIR-2'] <- 'NOIR2'

  feature.df$gene_name <- feature.df$gene
  feature.df$gene_name <- gsub('_2', '', feature.df$gene_name)
  feature.df <- rbind(feature.df, feature.nc)

} else {

  feature.df$gene_name <- feature.df$gene
  feature.df$gene_name <- gsub('_2', '', feature.df$gene_name)
}

feature.dt <- data.table(feature.df)


## Gene Regions (combination of genes and gene_clusters)
# Function to append cluster name if different from last gene
append_cluster <- function(genes, cluster) {
  if (tail(genes, 1) != cluster) {
    return(c(genes, cluster))
  } else {
    return(genes)
  }
}

# Apply the function by each gene_cluster
genes_and_clusters <- gene.clusters.all[, .(CompleteList = append_cluster(gene, gene_cluster)), by = gene_cluster]
genes_and_clusters[,gene_region := unlist(CompleteList)]

genes_and_clusters <- merge(genes_and_clusters,
                            unique(gene.clusters.all[,.(gene_cluster, Kinetic_class)]),
                            by='gene_cluster')

# Combine all results into one vector and ensure uniqueness
gene_regions  <- unique(genes_and_clusters$gene_region)


#### ####
##

### Analysis setting ####

## coverage
window_size <- 50
window_step <- 50


#### ####
##
