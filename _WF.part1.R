


#### Import bamfiles and rename alignments ####

require(BiocParallel)
safeBPParam <- function(nworkers) {
  if (.Platform$OS.type=="windows") {
    BiocParallel::SerialParam(bpparam())
  } else {
    BiocParallel::MulticoreParam(workers=nworkers, tasks = 20, progressbar = TRUE)
  }
}
if (.Platform$OS.type=="windows") { NULL } else {safeBPParam(nproc)}




## Import bamfiles
if (.Platform$OS.type=="windows") {

  bam.all.list <- purrr::map(bamfiles[],
                             dt.from.bam,
                             pattern=pattern,
                             flag=flag,
                             is.lortia = is.lortia,
                             crop.na.cigar = T,
                             rm.gaps.in.aln = rm.gaps.in.aln,
                             add.primes=T
  )
} else {
  bam.all.list <- bplapply(bamfiles[],
                           dt.from.bam,
                           pattern=pattern,
                           flag=flag,
                           is.lortia = is.lortia,
                           crop.na.cigar = T,
                           rm.gaps.in.aln = rm.gaps.in.aln,
                           add.primes=T
  )

}


bam.all <- rbindlist(bam.all.list)

## overwrite aln ID which was unique only for each sample
bam.all[,aln_ID := paste0(qname, '_', aln_nr)]

### IF MORE HOST GENOMES
if (rename_host_contigs) {
  ### Modify seqnames because of multiple genomes
  ## overwrite seqnames with genome_ID
  bam.all[!grepl('PK', sample) & !is.na(seqnames), seqnames:=paste0(seqnames, '_Rnor')]
  bam.all[ grepl('PK', sample) & !is.na(seqnames), seqnames:=paste0(seqnames, '_Sscrofa')]

  bam.all[,.N,by=.(sample, seqnames, strand)]
}
### fix viral contig
if (fix.viral.contigs) {
  bam.all[grepl(genome, seqnames),seqnames := genome]
}


### calculate mapped and non-mapped read counts
reads <- rbind( unique(bam.all[!is.na(seqnames),.(sample, qname)])[,is.na:=F],
                unique(bam.all[ is.na(seqnames),.(sample, qname)])[,is.na:=T])
readcounts <- reads[,.N,by=.(sample, is.na)]
readcounts <- dcast(readcounts, sample ~ is.na)

colnames(readcounts)[2] <- 'count'
readcounts[order(sample),]

print(readcounts)

#### ####
##



#### COVERAGE ANALYSIS ####
message('Analysing coverages ...')
source('coverage_virus.R')

### coverage summary plots
if(make.plots) {
  source('cov_stat_plots.R')
}
#### ####
##



#### write outputs ####

fwrite(readcounts, paste0(outdir, '/readcounts.tsv'), sep = '\t')

if (write.all) {
  fwrite(bam.all, paste0(outdir, '/bam.all.tsv'), sep = '\t')
  fwrite(mapped.cov, paste0(outdir, '/mapped.cov.tsv'), sep = '\t')
  fwrite(merged_cov, paste0(outdir, '/merged_cov.tsv'), sep = '\t')
}

#### ####
##

