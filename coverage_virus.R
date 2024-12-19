##### Viral coverages from LoRTIA



bam.cov.list <- bplapply(bamfiles,
                         cov.from.bam,
                         pattern=pattern,
                         param=param
)


bam.cov <- rbindlist(bam.cov.list)

cov.outfile <- paste0(outdir, '/bam.cov.rds')
saveRDS(bam.cov, cov.outfile)

bam.cov         <- readRDS(cov.outfile)


##### Analyse coverages

mapped.cov      <- merge(bam.cov, data.table(metafilt[,metacols]), by='sample', all=F)
#mapped.cov[,sample:=combine.groups]

### Normalize coverage
normalise_coverage <- function(mapped.cov, 
                               method = 'sum', norm_to = 1e6, 
                               contig_sizes = NULL, fasta.ref = NULL,
                               return = 'norm_cov') {
  
  # Calculate total counts per sample
  total_counts    <- mapped.cov[, .(total = sum(count)), by = sample]
  
  # Normalize coverage
  norm_mapped.cov <- mapped.cov[total_counts, on = 'sample']
  
  if (method == 'sum') {
    norm_mapped.cov[, norm_count := count / total * norm_to]
    
  } else if (method == 'genome') {
    
    # Check if contig_sizes is provided, otherwise use fasta.ref
    if (is.null(contig_sizes)) {
      if (is.null(fasta.ref)) {
        stop("Either contig_sizes or fasta.ref must be provided.")
      }
      
      # Load the reference FASTA file
      fasta  <- seqinr::read.fasta(fasta.ref)
      # Get contig sizes
      contig_sizes <- unlist(lapply(fasta, length))
    }
    contig_sizes <- data.table(seqnames=names(contig_sizes), contig_size=contig_sizes)
    
    norm_mapped.cov <- merge(norm_mapped.cov, contig_sizes, by='seqnames')
    
    norm_mapped.cov[, norm_count := count / (total / contig_size)]
    
  } else {
    stop("Invalid method. Currently, only 'sum' is supported.")
    
  }
  
  if (return == 'norm_cov') {
    norm_mapped.cov[, contig_size := NULL]
    
    # Drop the total column and return the result
    norm_mapped.cov[, total := NULL]
    norm_mapped.cov[, count := norm_count]
    norm_mapped.cov[, norm_count := NULL]
    
    return(norm_mapped.cov)
  } else if (return == 'summary') {
    
    colstouni <- setdiff(colnames(norm_mapped.cov), 
                         c('strand', 'count', 'norm_count', 'pos'))
    norm_cov_summary <- unique(norm_mapped.cov[,
                                               ..colstouni])
    norm_cov_summary[, average_coverage := total / contig_size]
  }
}

## Get the normalization base (total genome coverage / genome length) for each sample
norm_cov_summary <- normalise_coverage(mapped.cov, method = 'genome', fasta.ref = fasta.ref, return = 'summary')

## Normalise each position
norm.cov         <- normalise_coverage(mapped.cov, method = 'genome', fasta.ref = fasta.ref)


### Associate each position with a window

window.cov <- function(mapped_cov, contig_sizes = NULL, fasta.ref = NULL, window_size, window_step) {
  require(data.table)
  require(Biostrings)
  
  # Check if contig_sizes is provided, otherwise use fasta.ref
  if (is.null(contig_sizes)) {
    if (is.null(fasta.ref)) {
      stop("Either contig_sizes or fasta.ref must be provided.")
    }
    
    # Load the reference FASTA file
    fasta  <- seqinr::read.fasta(fasta.ref)
    # Get contig sizes
    contig_sizes <- unlist(lapply(fasta, length))
  }
  
  mapped_cov <- mapped_cov[,.(sample, seqnames, strand, pos, count)] %>% 
    spread(sample, count,  fill=0)  ;  setDT(mapped_cov)
  
  # Create a data.table with all positions 
  all_windows <- data.table(seqnames= names(contig_sizes), pos=1:contig_sizes)
  all_windows[, window_start := floor((pos - 1) / window_step) * window_step + 1]
  
  # Calculate window end
  all_windows[, window_end := pmin(window_start + window_size - 1, contig_sizes[seqnames])]
  
  all_windows <- data.table( plyr::rbind.fill(
    data.table(all_windows, strand = '+'), data.table(all_windows, strand = '-') ) )
  setDT(all_windows)
  
  # window start
  mapped_cov[, window_start := floor((pos - 1) / window_step) * window_step + 1]
  mapped_cov[, window_end := pmin(window_start + window_size - 1, contig_sizes[seqnames])]
  
  # Merge with mapped_cov
  merged_cov <- merge(all_windows, mapped_cov, 
                      by = c("seqnames", "strand", "pos", "window_start", "window_end"), 
                      all = TRUE, allow.cartesian=F)
  
  return(merged_cov)
}

### Normalized counts
merged_cov_norm <- window.cov(norm.cov, 
                              #contig_sizes = contig_sizes,
                              fasta.ref=fasta.ref,
                              window_size = window_size,
                              window_step = window_step
)

merged_cov_norm <- melt(merged_cov_norm,
                        #id.vars = setdiff(names(TR.gff.compare.merged.TR.counts), metafilt$sample),
                        measure.vars = unique(mapped.cov[,sample]),
                        variable.name = "sample",
                        value.name = "count")

merged_cov_norm <- merge(merged_cov_norm, metafilt, by='sample')
merged_cov_norm[is.na(count), count := 0]


### Raw counts
merged_cov <- window.cov(mapped.cov, 
                         #contig_sizes = contig_sizes,
                         fasta.ref=fasta.ref,
                         window_size = window_size,
                         window_step = window_step
)

merged_cov <- melt(merged_cov,
                   #id.vars = setdiff(names(TR.gff.compare.merged.TR.counts), metafilt$sample),
                   measure.vars = unique(mapped.cov[,sample]),
                   variable.name = "sample",
                   value.name = "count")

merged_cov <- merge(merged_cov, metafilt, by='sample')
merged_cov[is.na(count), count := 0]


### Summarise the coverages in each window (mean, sd, sum, etc) in each sample

win_coverage_summary <- function(merged_cov,  stat_by=metacols) {
  
  stat_by <- c('sample', 'seqnames', 'strand', 'window_start', "window_end")
  
  win.cov.sum <- 
    merged_cov[, 
               .(sum_coverage    = sum (count,   na.rm = TRUE),
                 mean_coverage   = mean(count,   na.rm = TRUE),
                 sd_coverage     = sd  (count,   na.rm = TRUE),
                 median_coverage = median(count, na.rm = TRUE),
                 min_coverage    = min (count,   na.rm = TRUE),
                 max_coverage    = max (count,   na.rm = TRUE)),
               by = .(sample, seqnames, strand, window_start, window_end)]
  
  win.cov.sum[,varcoeff_coverage := sd_coverage / mean_coverage]
  
  return(win.cov.sum)
  
}

win.cov.sum   <- win_coverage_summary(merged_cov)

win.cov.sum   <- merge(win.cov.sum, metafilt, by='sample')



### Statistics for the window coverage sums
win.cov.hpi.sum <- win.cov.sum[,
                               .(sum_coverage.hpi    = sum (sum_coverage),
                                 mean_coverage.hpi   = mean(sum_coverage,   na.rm = TRUE),
                                 sd_coverage.hpi     = sd  (sum_coverage,   na.rm = TRUE),
                                 median_coverage.hpi = median(sum_coverage, na.rm = TRUE),
                                 min_coverage.hpi    = min (sum_coverage,   na.rm = TRUE),
                                 max_coverage.hpi    = max (sum_coverage,   na.rm = TRUE)),
                               by = .(seqnames, strand, window_start, window_end, group, hpi, Time, cell_line)]

win.cov.hpi.sum[,varcoeff_coverage.hpi := sd_coverage.hpi / mean_coverage.hpi]


win.cov.varcoeff.stats <- win.cov.hpi.sum[,.(mean_varcoeff= mean(varcoeff_coverage.hpi, na.rm = TRUE),
                                             sd_varcoeff  = sd(varcoeff_coverage.hpi,   na.rm = TRUE) ),
                                          by = .(seqnames, strand, group, hpi, Time, cell_line)]

### Summarise the window coverages (sum) per hpi per genome
cov.hpi.sum.mean <- win.cov.hpi.sum[,
                                    .(mean_sum      = mean(sum_coverage.hpi,      na.rm = TRUE),
                                      mean_mean     = mean(mean_coverage.hpi,     na.rm = TRUE),
                                      mean_sd       = mean(sd_coverage.hpi,       na.rm = TRUE),
                                      mean_min      = mean(min_coverage.hpi,      na.rm = TRUE),
                                      mean_max      = mean(max_coverage.hpi,      na.rm = TRUE)),
                                    by=.(seqnames, strand, group, hpi, Time, cell_line)]

cov.hpi.sum.mean[,varcoeff_coverage := mean_sd / mean_mean]


### Summarise the non-windowed coverages (per genome) per hpi

norm_cov_summary.hpi <- norm_cov_summary[,
                                         .(sum_coverage.hpi    = sum (average_coverage),
                                           mean_coverage.hpi   = mean(average_coverage,   na.rm = TRUE),
                                           sd_coverage.hpi     = sd  (average_coverage,   na.rm = TRUE),
                                           median_coverage.hpi = median(average_coverage, na.rm = TRUE),
                                           min_coverage.hpi    = min (average_coverage,   na.rm = TRUE),
                                           max_coverage.hpi    = max (average_coverage,   na.rm = TRUE)),
                                         by = .(seqnames, group, hpi, Time, cell_line)]

norm_cov_summary.hpi[,varcoeff_coverage := sd_coverage.hpi / mean_coverage.hpi]

######### Get the statistics of the means of the normalized coverage for each window


#### mean the mean normalised count for each hpi
win.cov.hpi.mean <- win.cov.sum[,
                                .(sum_coverage.hpi    = sum (mean_coverage),
                                  mean_coverage.hpi   = mean(mean_coverage,   na.rm = TRUE),
                                  sd_coverage.hpi     = sd  (mean_coverage,   na.rm = TRUE),
                                  median_coverage.hpi = median(mean_coverage, na.rm = TRUE),
                                  min_coverage.hpi    = min (mean_coverage,   na.rm = TRUE),
                                  max_coverage.hpi    = max (mean_coverage,   na.rm = TRUE)),
                                by = .(seqnames, strand, window_start, window_end, group, hpi, Time, cell_line)]

win.cov.hpi.mean[,varcoeff_coverage.hpi := sd_coverage.hpi / mean_coverage.hpi]

#### summarise the normalised count for each hpi per genome
cov.hpi.mean <- win.cov.hpi.mean[,
                                 .(mean_varcoeff = mean(varcoeff_coverage.hpi,   na.rm = TRUE),
                                   mean_mean     = mean(mean_coverage.hpi,   na.rm = TRUE),
                                   mean_sd       = mean(sd_coverage.hpi,   na.rm = TRUE),
                                   min_sd        = mean(min_coverage.hpi,   na.rm = TRUE),
                                   max_sd        = mean(max_coverage.hpi,   na.rm = TRUE)),
                                 by=.(seqnames, strand, group, hpi, Time, cell_line)]

cov.hpi.mean[,varcoeff_coverage := mean_sd / mean_mean]

### Write outputs

fwrite(norm_cov_summary,     paste0(outdir, '/norm.cov.summary.tsv'), sep='\t')
fwrite(norm_cov_summary.hpi, paste0(outdir, '/norm.cov.summary.hpi.tsv'), sep='\t')


fwrite(win.cov.sum,      paste0(outdir, '/win.cov.sum.tsv'), sep='\t')
fwrite(win.cov.hpi.sum,  paste0(outdir, '/win.cov.hpi.sum.tsv'), sep='\t')
fwrite(cov.hpi.sum.mean, paste0(outdir, '/cov.hpi.sum.tsv'), sep='\t')


fwrite(win.cov.hpi.mean, paste0(outdir, '/win.cov.hpi.mean.tsv'), sep='\t')
fwrite(cov.hpi.mean,     paste0(outdir, '/cov.hpi.mean.tsv'), sep='\t')

