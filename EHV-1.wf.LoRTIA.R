

#### Settings and annotation

## everything was moved to here for simplicity and reproducibility

source('_WF.part0.R')



##

#### ####
##


#### START WORKFLOW !


##
######## WF Part 1. Importing and preprocessing alignments

##
#### settings ####

outdir  <- 'LoRTIA_virus'; try({ dir.create(outdir) })
project_config <- data.table(outdir = outdir)
fwrite(project_config, 'project_config.txt')

##
write.all <- T
rename_host_contigs <- F
fix.viral.contigs <- T
make.plots <- T
save.images <- F

#### ####
##




##
#### CAGE ####

include.cage <- T
if (include.cage) {
  outdir  <- 'CAGE'; try({ dir.create(outdir) })
  project_config <- data.table(outdir = outdir)
  #fwrite(project_config, 'project_config.txt')

  bamdir   <- "./CAGE/bam"
  pattern  <- '.bam$'
  bamfiles <- list.files(bamdir, pattern, recursive = T, full.names = T)

  is.lortia <- F
  flag  <- scanBamFlag(isSupplementaryAlignment=FALSE)
  param <- ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isSupplementaryAlignment=FALSE))
  rm.gaps.in.aln <- T

  meta.cage <- data.frame(sample = gsub('.*\\/', '', gsub(pattern, '', bamfiles)))

  meta.cage$cell_line <- 'RK-13'
  meta.cage$rep       <- c(1,1,1,2,2,2,3,3,3)
  meta.cage$group     <- gsub('_2023.*', '', meta.cage$sample)

  metafilt <- plyr::rbind.fill(meta.cage, metafilt)

  ### CECKPOINT -> stop
  source('_WF.part1.R')
  ##

  ## overwriting 'sample' to coerce the same sample from different runs
  bam.all <- merge(bam.all, meta.cage, by='sample')
  bam.all[,sample := group]
  bam.all <- bam.all[,c(1:15)] ## omit the rest of metadata for now

  bam.all[,correct_tes := T] # fifelse(grepl('correct', tag.l3) | grepl('correct', tag.r3), T, F)]
  bam.all[,correct_tss := T] # fifelse(grepl('correct', tag.l5) | grepl('correct', tag.r5), T, F)]

  bam.filt <- bam.all [!is.na(seqnames)]
  bam.filt <- bam.filt[flag %in% c(0, 16), ]

  ## corrrect the metadata accordingly
  metafilt[grepl('cage', metafilt$group), 'sample'] <- metafilt[grepl('cage', metafilt$group), 'group']
  metafilt <- unique.data.frame(metafilt)

  bam.filt.cage <- bam.filt

  fwrite(bam.filt.cage, paste0(outdir, '/CAGE.bam.filt.tsv'), sep = '\t')

  ### -> continue
  #mapped.cov       <- fread(paste0(outdir, '/mapped.cov.tsv'), na.strings = '')
  #merged_cov       <- fread(paste0(outdir, '/merged.cov.tsv'), na.strings = '')
  norm_cov_summary <- fread(paste0(outdir, '/norm.cov.summary.tsv'), na.strings = '')
  readcounts       <- fread(paste0(outdir, '/readcounts.tsv'), na.strings = '')

  readcounts.CAGE  <- readcounts
}

#### ####
##


##
#### dcDNA ####

flag  <- scanBamFlag(isSupplementaryAlignment=FALSE)
param <- ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isSupplementaryAlignment=FALSE))
rm.gaps.in.aln <- T

is.lortia <- T
bamdir  <- "../LoRTIA_stranded_only"
pattern <- '_stranded_only.bam'
bamfiles <- grep('.bai',
                 list.files(bamdir, pattern, recursive = T, full.names = T),
                 invert = T, value = T)


### import coverages
source('_WF.part1.R')

try({
  mapped.cov       <- fread(paste0(outdir, '/mapped.cov.tsv'), na.strings = '')
  merged_cov       <- fread(paste0(outdir, '/merged_cov.tsv'), na.strings = '')

  win.cov.sum      <- fread(paste0(outdir, '/win.cov.sum.tsv'), na.strings = '')
  win.cov.hpi.sum  <- fread(paste0(outdir, '/win.cov.hpi.sum.tsv'), na.strings = '')

  norm_cov_summary <- fread(paste0(outdir, '/norm.cov.summary.tsv'), na.strings = '')

  readcounts       <- fread(paste0(outdir, '/readcounts.tsv'), na.strings = '', header = T)
})

### remove unnecessary
rm(list = c('bam.all.list', 'bam.cov.list', 'mapped.cov', 'merged_cov', 'norm.cov'))

#### ####
##



##
#### dRNA ####

is.lortia <- F
bamdir  <- "C:/data/EHV-1/rebasecall/mapped"
pattern <- '.bam'
bamfiles <- grep('.bai',
                 list.files(bamdir, pattern, recursive = T, full.names = T),
                 invert = T, value = T)

### import coverages
source('_WF.part1.R')

try({
  mapped.cov       <- fread(paste0(outdir, '/mapped.cov.tsv'), na.strings = '')
  merged_cov       <- fread(paste0(outdir, '/merged_cov.tsv'), na.strings = '')

  win.cov.sum      <- fread(paste0(outdir, '/win.cov.sum.tsv'), na.strings = '')
  win.cov.hpi.sum  <- fread(paste0(outdir, '/win.cov.hpi.sum.tsv'), na.strings = '')

  norm_cov_summary <- fread(paste0(outdir, '/norm.cov.summary.tsv'), na.strings = '')

  readcounts       <- fread(paste0(outdir, '/readcounts.tsv'), na.strings = '', header = T)
})

### remove unnecessary
rm(list = c('bam.all.list', 'bam.cov.list', 'mapped.cov', 'merged_cov', 'norm.cov'))

#### ####
##




##
#### ####

### Combine w CAGE
if(include.cage) {
  readcounts       <- rbind(readcounts, readcounts.CAGE)
  #bam.filt         <- rbind(bam.filt,   bam.filt.cage)
}


### SAVE / LOAD IMAGE
gc()
if(save.images) {
  save.image(paste0(outdir, '.P1.RData'))
  load(paste0(outdir, '.P1.RData'))
}

#### ####
##


######## WF Part 2. Read clustering into TransFrags
####  ####

if(is.lortia) {
  bam.all[,correct_tes := fifelse(grepl('correct', tag.l3) | grepl('correct', tag.r3), T, F)]
  bam.all[,correct_tss := fifelse(grepl('correct', tag.l5) | grepl('correct', tag.r5), T, F)]
}

## filter out reads with supplementary alignments, or
## these wer filtered out previously?
# bam.filt[qname]
bam.filt <- bam.all[!is.na(seqnames)]

rm('bam.all')

fwrite(bam.filt, paste0(outdir, '/bam.filt.tsv'), sep='\t')

bam.filt <- fread(paste0(outdir, '/bam.filt.tsv'), na.strings = '')


## make TransFrags
source('makeTRs.R')
TR.data <- TR.uni
TR.data[, TR_start := min(start), by =.(TR_ID)][, TR_end := max(end), by= .(TR_ID)]
TR.data <- dcast(TR.data, seqnames + strand + TR_ID + TR_start + TR_end + exon_combination ~ exon_rank, value.var = 'exon_ID')
TR.data <- merge(TR.data, TR.counts, by='TR_ID', all=T)

fwrite(TR.uni,  paste0(outdir, '/TR.uni.tsv'),  sep='\t')
fwrite(TR.data, paste0(outdir, '/TR.data.tsv'), sep='\t')

TR.data <- fread(paste0(outdir, '/TR.data.tsv'), na.strings = '')
TR.uni  <- fread(paste0(outdir, '/TR.uni.tsv'),  na.strings = '')

## CAGE
#CAGE.TR.data <- TR.data
#

## export read-TRs
TR.reads.gfffile <- paste0(outdir, '/TR.reads.gff2')
source('export.gffs.R')

TR.gff <- data.table(as.data.frame(rtracklayer::import.gff2(TR.reads.gfffile)))
TR.EX  <- fread(paste0(outdir, '/TR.EX.tsv'))


## Import REF TR-s
viral.ref <- "NC_001491.2.mRNA_corrected.gff3"
source('import.ref.TRs.R')
#fwrite(TR.Ref.data, paste0(outdir, '/TR.Ref.data.tsv'), sep = '\t')
#TR.Ref.data <- fread(paste0(outdir, '/TR.Ref.data.tsv'))

### remove
rm(list = c('dt', 'reads'))


## TR counts with adapter info per sample
aln.uni      <- unique(bam.TR[,.(seqnames, strand, aln_ID, TR_ID, TR_start, TR_end, correct_tss, correct_tes, sample)])

## Validate 3-primes based on reference transcripts TES
valid.TES.win <- 10
validate_TES  <- T

if (validate_TES) {

  TR.ref[,transcript_prime3 := fifelse(strand == '+', transcript_end, transcript_start)]
  TR.ref[,transcript_prime5 := fifelse(strand == '-', transcript_end, transcript_start)]

  valid.prime3 <- unique(TR.ref[,.(seqnames, transcript_start, transcript_end, transcript_prime3, transcript_prime5, strand)])
  valid.prime3 <- valid.prime3[,.(seqnames, strand, start = transcript_prime3 - valid.TES.win, end = transcript_prime3 + valid.TES.win)]
  valid.prime3 <- unique(valid.prime3)
  valid.prime3[,valid_tes := T]

  aln.uni[,TR_prime3 := fifelse(strand == '+', TR_end, TR_start )]
  aln.uni[,TR_prime5 := fifelse(strand == '-', TR_end, TR_start )]
  aln.uni[,start := TR_prime3]
  aln.uni[,end   := TR_prime3]

  prime3.TR.ov <- foverlaps2(aln.uni, valid.prime3, by.x=c('seqnames', 'strand', 'start', 'end'), by.y=c('seqnames', 'strand', 'start', 'end'), minoverlap = 1)
  prime3.TR.ov <- prime3.TR.ov[,.(seqnames, strand, TR_start, TR_end, TR_prime3, TR_prime5, correct_tss, correct_tes, aln_ID, TR_ID, start, end, sample)]
  prime3.TR.ov <- merge(prime3.TR.ov, valid.prime3, by=c('seqnames', 'strand', 'start', 'end'), all.x=T)

  prime3.TR.ov[, start  := NULL]
  prime3.TR.ov[, end    := NULL]
  aln.uni[, start       := NULL]
  aln.uni[, end         := NULL]


  prime3.TR.ov <- merge(prime3.TR.ov, aln.uni, by=colnames(aln.uni), all=T)
  prime3.TR.ov[,valid_tes := fifelse(is.na(valid_tes), F, T)]

  prime3.valid.corr.freq <- prime3.TR.ov[,.N, by=.(sample, correct_tes, valid_tes)]

  ggplot(prime3.valid.corr.freq) +
    geom_col(aes(x=sample, y=N, fill=correct_tes), color='black') +
    coord_flip() +
    facet_wrap(~valid_tes, nrow=1) +
    theme_bw() +
    ggtitle('3-prime end result, according to ref mRNA TES')

  ggsave(file.path(outdir, 'Ref_mRNA_TES_validation.jpg'), height = 12, width = 9)

  prime3.TR.ov[, correct_tes := fifelse(valid_tes == T | correct_tes == T, T, F)]
  prime3.TR.ov[, valid_tes   := NULL]

  keyby <- colnames(prime3.TR.ov)
  prime3.TR.ov <- unique(prime3.TR.ov[,]) ##.(), by=.()]

  if(nrow(aln.uni) == nrow(prime3.TR.ov)) {
    aln.uni      <- prime3.TR.ov
  } else { stop() }

}

TR.adapt.count <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_ID, TR_start, TR_end, correct_tss, correct_tes, sample)]
fwrite(TR.adapt.count, paste0(outdir, '/TR.adapt.count.tsv'), sep = '\t')

TR.adapt.count <- fread(paste0(outdir, '/TR.adapt.count.tsv'), na.strings = '')

TR.counts.sp <- dcast(TR.adapt.count, TR_ID + correct_tss + correct_tes ~ sample, value.var = 'count', fill=0)
##

#orf5 <- bam.TR[TR_end >= 46795 & TR_end <= 46995 & grepl('12h', sample) & strand == '-',]


#### ####
##

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P2.RData'))
  load(paste0(outdir, '.P2.RData'))
}


##
#### WF part 3. Ref Transcript count analysis ####

### Run GFF-compare on each ref TR separately,
## import results and calculate distances
source('run_GFF.COMPARE.R')

### Analyse GFF-compare results
## read in GFF-compare results
all.merged.result_gff.compare  <- fread(paste0(outdir, "/all.merged.result_gff.compare.tsv"), na.strings = '')

# Update the results as transcript IDs in the reference annotation was changed
#source('update_TR.ref.IDs.R')

## find the closest ref-TR for each query
## categorise non-equal matches
thresh.eq.prime5 <- 10
thresh.eq.prime3 <- 10
thresh.eq.junc   <- 2
source('analyse_GFF.COMPARE.R')
## import results
best.merged.result_gff.compare <- fread(paste0(outdir, "/best.merged.result_gff.compare.tsv"), na.strings = '')

### Summarise results
source('summarise_GFF.COMPARE.R')

TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"), na.strings = '')

### include LoRTIA adapter info to Transcripts
adapt.TR <- aln.uni[,.N,by=.(TR_ID, correct_tss, correct_tes)]
adapt.TR[,correct_tss := paste0('correct_tss::', as.character(correct_tss))]
adapt.TR[,correct_tes := paste0('correct_tes::', as.character(correct_tes))]
adapt.TR[,adapter := paste0(correct_tes, ';', correct_tss)]
adapt.TR.sp <- dcast(adapt.TR, TR_ID ~ adapter, value.var = 'N', fill=0)


##
all.merged.result_gff.compare      <- merge(all.merged.result_gff.compare,      adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID')
best.merged.result_gff.compare     <- merge(best.merged.result_gff.compare,     adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID')
TR.gff.compare.merged.TR.counts.gt <- merge(TR.gff.compare.merged.TR.counts.gt, adapt.TR.sp, by.x='transcript_id',  by.y='TR_ID')


#### ####
##

### remove
rm(list=c('all.merged.result_gff.compare', 'TR.gff.compare.merged.TR', 'TR.gff.compare.merged.TR.counts'))
gc()

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P3.RData'))
  load(paste0(outdir, '.P3.RData'))
}


##
#### WF part 4. Count genes and Transcripts ####

#source('_WF.part0.R')

res.dir <- outdir; try({ dir.create(res.dir) })

## Import canonic TSS and TES data
#gene.clusters.all <- fread('PRV.genes.TSS.TES.txt')
gene.clusters.all[,gene.region.start := as.integer(ifelse(strand == '+', TSS.canonic, TES.canonic))]
gene.clusters.all[,gene.region.end   := as.integer(ifelse(strand == '-', TSS.canonic, TES.canonic))]
### Dereplicate spliced genes !!!
gene.clusters.all <- unique(gene.clusters.all[,.(seqnames, TSS.canonic, TES.canonic, cluster_TES, strand, gene, gene_cluster, Kinetic_class, TSS.win.start, TSS.win.end, TES.win.start, TES.win.end)])



## filter some genes ?
genes_to_filter <- NULL # 'CTO-L'

## treat spliced transcripts differently?
check.spliced.TRs <- F

## Count genes based on TSS and/or TES
source('count.genes.v2.R')
## OK


## Count transcripts and combine with gene (and cluster) counts
source('count.transcripts.R')

######## PUT THIS INTO PLACE
source('CAGE_TRs.R')



gc()

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P4.RData'))
  load(paste0(outdir, '.P4.RData'))
}

##### FINISH #######




## Count coverages of GENE REGIONS (FROM TSS TO TES)
source('count.gene.cov.R')
cov.gene.ov.counts <- fread(paste0(outdir, '/cov.gene.ov.counts.tsv'), sep = '\t')


## Count coverages of CDS
source('count.CDS.cov.R')
cov.CDS.ov.counts <- fread(paste0(outdir, '/cov.CDS.ov.counts.tsv'), sep = '\t')


## Count overlaps of reference transcripts
source('count.overlaps.R')


## Finding ANY 10 Overlaps between exons and GENE REGIONS (FROM TSS TO TES)
source('count.generegions.any10.R')
TR.any10.generegion.counts <- fread(paste0(outdir, '/TR.any10.generegion.counts.tsv'), sep = '\t')


## Finding ANY 10-nt Overlaps between exons and ORFs
source('count.ORFs.any10.R')
TR.any10.ORF.counts <- fread(paste0(outdir, '/TR.any10.ORF.counts.tsv'), sep = '\t')


## Finding WITHIN Overlaps between exons and ORFs
source('count.ORFs.within.R')
TR.within.ORF.counts <- fread(paste0(outdir, '/TR.within.ORF.counts.tsv'), sep = '\t')

#### Analysis ready ! ####
##

gc()

### SAVE / LOAD IMAGE
if(save.images) {
  save.image(paste0(outdir, '.P4.RData'))
  load(paste0(outdir, '.P4.RData'))
}

##
#### Plotting ####

### GFF-compare results
source('plot_GFF.COMPARE.R')

### Reference transcript's overlaps
source('plot.overlaps.R')

### Venn diagram of shared transcripts in the cell lines
source('plot.TR.Venn.R')

#### ####
##


##
#### Coverage and read end plots ####

### CAGE
CAGE.TR.data <- fread('CAGE/TR.data.tsv')
CAGE.TR.data <- merge(CAGE.TR.data, meta.cage, by='sample')
setnames(CAGE.TR.data, new=c('start', 'end'), old=c('TR_start', 'TR_end'),skip_absent=TRUE)

CAGE.TR.data[,prime3 := ifelse(strand == '+', end,   start)]
CAGE.TR.data[,prime5 := ifelse(strand == '+', start, end)]

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime5')
prime5.counts <- CAGE.TR.data[count>0, .(count=sum(count)), by=cols_to_group][order(seqnames, strand, prime5)]
setnames(prime5.counts, 'prime5', 'pos')
prime5.counts[,endtype := 'prime5']

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime3')
prime3.counts <- CAGE.TR.data[count>0, .(count=sum(count)), by=cols_to_group][order(seqnames, strand, prime3)]
setnames(prime3.counts, 'prime3', 'pos')
prime3.counts[,endtype := 'prime3']

prime.counts  <- rbind(prime5.counts, prime3.counts)
prime.counts.CAGE <- prime.counts
#prime.counts.CAGE$group <- paste0('CAGE_', prime.counts.CAGE$cell_line, '_', prime.counts.CAGE$hpi)

fwrite(prime.counts.CAGE, paste0(outdir, '/prime.counts.CAGE.tsv'))

### Read end counts from read-transcripts
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"))

TR.gff.compare.merged.TR.counts.gt[,prime3 := fifelse(strand == '+', end,   start)]
TR.gff.compare.merged.TR.counts.gt[,prime5 := fifelse(strand == '+', start, end)]

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime5')
prime5.counts <- TR.gff.compare.merged.TR.counts.gt[count>0,.(count=sum(count)),by=cols_to_group][order(seqnames, strand, prime5)]
setnames(prime5.counts, 'prime5', 'pos')
prime5.counts[,endtype := 'prime5']

cols_to_group <- c(metacols, 'seqnames', 'strand', 'prime3')
prime3.counts <- TR.gff.compare.merged.TR.counts.gt[count>0,.(count=sum(count)),by=cols_to_group][order(seqnames, strand, prime3)]
setnames(prime3.counts, 'prime3', 'pos')
prime3.counts[,endtype := 'prime3']

#### OR ! From LoRTIA outpu, to include adapter info !!!

### Read end counts from bamfiles, separated by adapters
bam.counts <- bam.TR
bam.counts[,prime3 := fifelse(strand == '+', end,   start)]
bam.counts[,prime5 := fifelse(strand == '+', start, end)]
bam.counts[,correct_tss := fifelse(grepl('correct', tag.r5) | grepl('correct', tag.l5), T,   F)]
bam.counts[,correct_tes := fifelse(grepl('correct', tag.r3) | grepl('correct', tag.l3), T,   F)]

bam.counts[,TR_prime3 := fifelse(strand == '+', TR_end,   TR_start)]
bam.counts[,TR_prime5 := fifelse(strand == '+', TR_start, TR_end)]

cols_to_group <- c('sample', 'seqnames', 'TR_ID', 'aln_ID', 'strand', 'TR_start', 'TR_end', 'TR_prime3', 'TR_prime5', 'correct_tss', 'correct_tes')

aln.uni <- unique(bam.counts[,..cols_to_group])

prime5.bam.counts <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_prime5, correct_tss, sample)]; setnames(prime5.bam.counts, 'TR_prime5', 'pos')
prime3.bam.counts <- aln.uni[,.(count=.N), by=.(seqnames, strand, TR_prime3, correct_tes, sample)]; setnames(prime3.bam.counts, 'TR_prime3', 'pos')

prime5.bam.counts <- merge(prime5.bam.counts, metafilt[,metacols], by='sample')
prime5.bam.counts[, endtype := 'prime5']

prime3.bam.counts <- merge(prime3.bam.counts, metafilt[,metacols], by='sample')
prime3.bam.counts[, endtype := 'prime3']

##
prime3.counts <- prime3.bam.counts
prime5.counts <- prime5.bam.counts

fwrite(prime3.counts, paste0(outdir, '/prime3.counts.tsv'))
fwrite(prime5.counts, paste0(outdir, '/prime5.counts.tsv'))


### Combine
prime.counts  <- rbind(prime5.counts, prime3.counts)
fwrite(prime.counts, paste0(outdir, '/prime.counts.tsv'))


### Combine w CAGE
prime5.counts <-  rbind(prime.counts[endtype == 'prime5'], prime.counts.CAGE[endtype == 'prime5'])
prime3.counts <-  rbind(prime.counts[endtype == 'prime3'] ) #, prime.counts.CAGE[endtype == 'prime3'])
prime.counts  <-  rbind(prime.counts, prime.counts.CAGE)


#### Make 5- and 3-prime end plots
prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')
prime3.counts <- fread(paste0(outdir, '/prime3.counts.tsv'), na.strings = '')

## Filter out false adaptered reads?
#prime5.counts <- prime5.counts[correct_tss == T,]
#prime3.counts <- prime3.counts[correct_tes == T,]

## OR !! use ref_mRNAs to accept 3-primes?
TR.ref[,transcript_prime3 := fifelse(strand == '+', transcript_end, transcript_start)]
TR.ref[,transcript_prime5 := fifelse(strand == '-', transcript_end, transcript_start)]

valid.prime3 <- unique(TR.ref[,.(seqnames, transcript_start, transcript_end, transcript_prime3, transcript_prime5, strand)])
valid.prime3 <- valid.prime3[,.(seqnames, strand, start = transcript_prime3 - 10, end = transcript_prime3 + 10)]

prime3.counts[,start := pos]
prime3.counts[,end   := pos]
prime3.TR.ov <- foverlaps2(prime3.counts, valid.prime3, by.x=c('seqnames', 'strand', 'start', 'end'), by.y=c('seqnames', 'strand', 'start', 'end'), minoverlap = 1)
prime3.TR.ov <- unique(prime3.TR.ov[,.(seqnames, strand, pos, correct_tes, endtype, sample, hpi, Time, cell_line, group, count)])

valid.prime3 <- unique(prime3.TR.ov[,.(seqnames, strand, pos, valid_tes=T)])

prime3.counts <- merge(prime3.counts, valid.prime3, by=c('seqnames', 'strand', 'pos'), all.x=T)
prime3.counts[,valid_tes := fifelse(is.na(valid_tes), F, T)]

prime3.valid.corr.freq <- prime3.counts[,.(sum_count = sum(count)), by=.(correct_tes, valid_tes)]

prime3.counts[,correct_tes := fifelse(valid_tes == T | correct_tes == T, T, F)]
prime3.counts[,valid_tes := NULL]


## Overwrite?
fwrite(prime3.counts, paste0(outdir, '/prime3.counts.tsv'), sep = '\t')
fwrite(prime5.counts, paste0(outdir, '/prime5.counts.tsv'), sep = '\t')

## Ovewrite TR count table?
valid.prime3.TR <- unique(merge(valid.prime3, TR.uni, by.x=c('seqnames', 'strand', 'pos'), by.y=c('seqnames', 'strand', 'prime3'), all.x=T))
valid.prime3.TR <- valid.prime3.TR[,.(seqnames, strand, TR_ID, valid_tes)]
TR.adapt.count  <- merge(TR.adapt.count, valid.prime3.TR, by=c('seqnames', 'strand', 'TR_ID'), all.x=T)
TR.adapt.count[,valid_tes   := fifelse(is.na(valid_tes), F, T)]
TR.adapt.count[,correct_tes := fifelse(valid_tes == T | correct_tes == T, T, F)]
TR.adapt.count[,valid_tes   := NULL]

TR.adapt.count <- TR.adapt.count[,.(count=sum(count)), by=.(seqnames, strand, TR_ID, correct_tss, correct_tes, sample)]

fwrite(TR.adapt.count, paste0(outdir, '/TR.adapt.count.tsv'), sep = '\t')

## include adapter counts
TR.counts.sp <- dcast(TR.adapt.count, TR_ID + correct_tss + correct_tes ~ sample, value.var = 'count', fill = 0)


## use CAGE to accept 5-prime sites?
#     Answer: No.


##### NEW COUNT TABLE

## prime5 and prime3 counts from TR-table, considering adapters

prime.counts  <- merge(TR.adapt.count, metafilt, by='sample')
prime.counts[,start  := TR_start]
prime.counts[,end    := TR_end  ]
prime.counts[,prime5 := fifelse(strand == '+', start, end)]
prime.counts[,prime3 := fifelse(strand == '+', end,   start)  ]

prime3.counts <- prime.counts[,.(count=sum(count)), by=.(seqnames,	strand,	correct_tes, prime3, sample, hpi, Time,	cell_line,	group)]
prime3.counts[, endtype := 'prime3'][, pos := prime3]
prime3.counts[,start  := pos]
prime3.counts[,end    := pos]


prime5.counts <- prime.counts[,.(count=sum(count)), by=.(seqnames,	strand,	correct_tss, prime5, sample, hpi, Time,	cell_line,	group)]
prime5.counts[, endtype := 'prime5'][, pos := prime5]
prime5.counts[,start  := pos]
prime5.counts[,end    := pos]

## Mean coverage from stranded only bamfiles directly
cov.counts <- merged_cov[,.(count=mean(count)), by=.(seqnames,	strand,	pos, sample, hpi, Time,	cell_line,	group)]
cov.counts[,start  := pos]
cov.counts[,end    := pos]


## write
fwrite(prime3.counts, paste0(outdir, '/prime3.counts.tsv'), sep = '\t')
fwrite(prime5.counts, paste0(outdir, '/prime5.counts.tsv'), sep = '\t')
fwrite(cov.counts,    paste0(outdir, '/cov.counts.tsv'),    sep = '\t')




##### START PLOTTING FROM HERE
#cov.counts    <- fread(paste0(outdir, '/cov.counts.tsv'),    na.strings = '')

## read back in
prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')
prime3.counts <- fread(paste0(outdir, '/prime3.counts.tsv'), na.strings = '')



## filter now
correct.only <- T
if(correct.only) {
  prime5.counts <- prime5.counts[correct_tss == T,]
  prime3.counts <- prime3.counts[correct_tes == T,]
}


### normalize for viral read counts?
## based on correct adaptered-only
norm <- F
if(norm) {
  colsby <- colnames(prime5.counts)[!colnames(prime5.counts) %in% c('pos', 'start', 'end', 'prime5', 'count')]
  prime5.counts[, sum_count := sum(count), by=colsby]
  prime5.counts[, count     := round((count/ sum_count) * 100, 3)]
  prime5.counts[, .(sum_count = sum(count)), by=colsby]
  prime5.counts[, sum_count := NULL]

  colsby <- colnames(prime3.counts)[!colnames(prime3.counts) %in% c('pos', 'start', 'end', 'prime3', 'count')]
  prime3.counts[, sum_count := sum(count), by=colsby]
  prime3.counts[, count     := round((count/ sum_count) * 100, 3)]
  prime3.counts[, .(sum_count = sum(count)), by=colsby]
  prime3.counts[, sum_count := NULL]

}

## check if all the samples have the same library size
all(round(prime5.counts[,.(sum_count = (sum(count))), by=.(sample)][,sum_count], 0) == 1)


##calculate the mean? otherwise, the plot will SUM the counts based on the "group" column
calc.mean <- F
##
if(calc.mean) {
  prime3.counts <- prime3.counts[,.(count = round(mean(count), 2)), by=.(seqnames, strand, pos, correct_tes, hpi, Time, cell_line, group, endtype)]
  prime5.counts <- prime5.counts[,.(count = round(mean(count), 2)), by=.(seqnames, strand, pos, correct_tss, hpi, Time, cell_line, group, endtype)]

  #prime5.counts <- unique(prime5.counts[,.(seqnames, strand, pos, correct_tss, count, hpi, Time, cell_line, group, endtype)])
  #prime3.counts <- unique(prime3.counts[,.(seqnames, strand, pos, correct_tes, count, hpi, Time, cell_line, group, endtype)])

}

##
prime5.counts[,prime5 := pos][,prime5 := NULL]
prime3.counts[,prime3 := pos][,prime3 := NULL]

##### plotting settings

##
plot.all.together <- T

## grouping
combine.groups <- 'hpi' ## 'sample' ## OR:

## adapter settings
adapter.setting <- 'v4'

##
tolower_gene_name <- F



######## REVIEW

#### Normalized
file_name_suffix <- 'Mean_correct_norm'
##
ylim       <- NULL #c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)
#ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW_ori.tif'), plot = FigY, height = 20, width = 36, limitsize = F, dpi=300)
Fig1B <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

Fig2B <- Fig2

ggsave(paste0(fig.dir, '/', 'Figure 2B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)



#### ylim 0-500
file_name_suffix <- 'Mean_correct_ylim500'
##
ylim       <- c(0, 500) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'Figure 1A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig1A <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

ggsave(paste0(fig.dir, '/', 'Figure 2A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig2A <- Fig2




#### ylim 0-50
file_name_suffix <- 'Mean_correct_norm'
##
ylim       <- c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 1B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)
#ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW_ori.tif'), plot = FigY, height = 20, width = 36, limitsize = F, dpi=300)
Fig1B <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

Fig2B <- Fig2

ggsave(paste0(fig.dir, '/', 'SuppFig 2B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)


#### ylim 0-5000
file_name_suffix <- 'Mean_correct_ylim500'
##
ylim       <- c(0, 5000) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 1A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig1A <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 2A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig2A <- Fig2




######### ORIGINAL

#### 0-5000
file_name_suffix <- 'Mean_correct_norm'
##
ylim       <- c(0, 5000) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 1C_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)
#ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW_ori.tif'), plot = FigY, height = 20, width = 36, limitsize = F, dpi=300)
Fig1B <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

Fig2B <- Fig2

ggsave(paste0(fig.dir, '/', 'SuppFig 2C_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)


#### ylim500
file_name_suffix <- 'Mean_correct_ylim500'
##
ylim       <- c(0, 500) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 1D_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig1A <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 2D_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig2A <- Fig2











####### Supp Figs

#### ylim50
file_name_suffix <- 'Mean_correct_ylim50'
##
ylim       <- c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
source('plot_prime5.R')

SFig1A <- Fig1

## prime3
source('plot_prime3_area_settings.R')
source('plot_prime3.R')

SFig2A <- Fig2


#### ylim5000
file_name_suffix <- 'Mean_correct_ylim5000'
##
ylim       <- c(0, 5000) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
source('plot_prime5.R')

SFig1B <- Fig1

## prime3
source('plot_prime3_area_settings.R')
source('plot_prime3.R')

SFig2B <- Fig2

#
#### Coverage plots

##
file_name_suffix <- 'Mean_correct_ylim5000'
ylim       <- c(-5000, 5000) # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3A <- Fig3

##
file_name_suffix <- 'Mean_correct_ylim500'
ylim       <- c(-500, 500) # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3B <- Fig3

##
file_name_suffix <- 'Mean_correct_ylim50'
ylim       <- c(-50, 50) # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3C <- Fig3

##
file_name_suffix <- 'Mean_correct_norm'
ylim       <- NULL # c(-1000, 1000)

source('plot_coverage_area_settings.R')
source('plot_coverage.R')

SFig3D <- Fig3

#

## done
######


#### Small Figures for the article
fig.width  <- 30
fig.height <- 20*2

Fig1   <- cowplot::plot_grid(Fig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             Fig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 1.small.jpg', Fig1, height = fig.height, width = fig.width, limitsize = F)

Fig2  <- cowplot::plot_grid(Fig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            Fig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 2.small.jpg', Fig2, height = fig.height, width = fig.width, limitsize = F)


SFig1 <- cowplot::plot_grid(SFig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 1.small.jpg', SFig1, height = fig.height, width = fig.width, limitsize = F)


SFig2 <- cowplot::plot_grid(SFig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 2.small.jpg', SFig2, height = fig.height, width = fig.width, limitsize = F)



#### Normal Figures for the article
fig.width  <- 55
fig.height <- 20*2

Fig1   <- cowplot::plot_grid(Fig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             Fig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 1.jpg',  Fig1, height = fig.height, width = fig.width, limitsize = F)

Fig2  <- cowplot::plot_grid(Fig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            Fig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Figure 2.jpg',  Fig2, height = fig.height, width = fig.width, limitsize = F)


SFig1 <- cowplot::plot_grid(SFig1A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig1B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 1.jpg', SFig1, height = fig.height, width = fig.width, limitsize = F)


SFig2 <- cowplot::plot_grid(SFig2A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            SFig2B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                            ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 2.jpg', SFig2, height = fig.height, width = fig.width, limitsize = F)


##### Sup Fig 3
fig.width  <- 30
fig.height <- 20*3

SFig3  <- cowplot::plot_grid(SFig3A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3C + draw_label("c", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 3.small.jpg', SFig3, height = fig.height, width = fig.width, limitsize = F)

fig.width  <- 55
fig.height <- 20*3

SFig3  <- cowplot::plot_grid(SFig3A + draw_label("a", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3B + draw_label("b", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             SFig3C + draw_label("c", x = 0.025, y = 1, hjust = 0.1, vjust = 1, fontface = 'bold', size = 24),
                             ncol=1)
ggsave('../EHV-1 dynamic article/Figures/Supp Fig 3.jpg', SFig3, height = fig.height, width = fig.width, limitsize = F)























## CTO vicinity
source('plot_prime5_bar_settings.R')
source('plot_prime5.R')

source('plot_prime3_bar_settings.R')
source('plot_prime3.R')


###

## coverage
cov.counts <- data.table(data.frame(win.cov.sum[,.(pos = round(mean(c(window_end, window_start)),0), count=sum_coverage),
                                                by=.(seqnames, sample, strand, group, hpi, Time, window_start )]))

source('plot_coverage_area_bin50_settings.R')
source('plot_prime3_area_bin50.R')


#### CAGE anyalysis
source('CAGE.R')


#### DBSCAN clustering ####

TR.prime3 <- prime3.counts
TR.prime3 <- TR.prime3[,.(seqnames, strand, sample, count, TR.prime3 = pos, hpi, cell_line, group)]

##### Cluster the counts of the prime3, in th hpi12 samples with count as weight
## OMIT CAGE FOR NOW
TR.prime3.sum <- unique(TR.prime3[cell_line == 'PK-15' & hpi == '12h' & !grepl('CAGE', group),
                                  .(sum_count = sum(count)),
                                  by=.(seqnames, strand, TR.prime3, hpi)])

fwrite(TR.prime3.sum, paste0(outdir, '/DBScan.prime3.cluster.data.tsv'), sep='\t')

TR.prime3.sum <- fread(paste0(outdir, '/DBScan.prime3.cluster.data.tsv'), na.strings = '')


TR.prime5 <- prime5.counts
TR.prime5 <- TR.prime5[,.(seqnames, strand, sample, count, TR.prime5 = pos, hpi, cell_line, group)]

##### Cluster the counts of the prime5, in th hpi12 samples with count as weight
## OMIT CAGE FOR NOW
TR.prime5.sum <- unique(TR.prime5[cell_line == 'PK-15' & hpi == '12h' & !grepl('CAGE', group),
                                  .(sum_count = sum(count)),
                                  by=.(seqnames, strand, TR.prime5, hpi)])

fwrite(TR.prime5.sum, paste0(outdir, '/DBScan.prime5.cluster.data.tsv'), sep='\t')

TR.prime5.sum <- fread(paste0(outdir, '/DBScan.prime5.cluster.data.tsv'), na.strings = '')



#### Initial Clustering
## Parameters
eps    <- 20
minPts <- 5

DT <- unique(TR.prime3.sum[,.(seqnames, strand, position=TR.prime3, count=sum_count)])

source('DBScan_primary.R')


#### Secondary Clustering
## Parameters
eps_sec    <- 10
minPts_sec <- 5

DT.clust     <- fread(paste0(outdir, "/DBSCAN_primary_clusters_", "eps.", eps, "_minPts.", minPts, '.tsv'), na.strings = '')
DT.clust.uni <- fread(paste0(outdir, "/DBSCAN_UNIQ_primary_clusters_", "eps.", eps, "_minPts.", minPts, '.tsv'), na.strings = '')
DT.clust.uni <- DT.clust.uni[order(cluster_center)]


source('DBScan_secondary.R')

DT.sec.clust <- fread(paste0(outdir, "/DBSCAN_secondary_clusters_", "eps.", eps_sec, "_minPts.", minPts_sec, '.tsv'), na.strings = '')



#### ####
##


##
#### Write outputs ####

save.image('PRV.LoRTIA.all.RData')

#### ####
##






stop()



### Make the plots
source('knit.docs.R')

## 1.) Coverage + TES, TSS windows
source('cluster.plot1_coverages.and.windows.R')

## 2.) Coverage + TR-reads ## -->> How to plot reads in each cell line?
#source('cluster.plot2_coverages.and.TR.reads.R')

## 3.) Coverage + TR-annot
source('cluster.plot3_coverages.and.TR.ref.R')

## 4.) TR-annot + TR-reads ## -->> How to plot reads in each cell line?
#source('cluster.plot3_coverages.and.TR.ref.R')

## 5.) TES, TSS windows + TR-annot ##
source('cluster.plot5_windows.and.TR.ref.R')


## OK

source('TSS_abund.LoRTIA.Rmd')
## KNIT !

#### ####
##













### host reads
viral_reads <- melt(gene.sample_count.sp, variable.name = 'sample', value.name = 'count')
viral_reads <- viral_reads[,.(viral_read_count = sum(count)), by=sample]

total_reads <- fread('../fastq_read_counts.tsv')

all_read_counts <- merge(total_reads, viral_reads, by='sample')
all_read_counts[,host_read_count := read_count - viral_read_count]

fwrite(all_read_counts, 'all_read_counts.tsv', sep = '\t')

