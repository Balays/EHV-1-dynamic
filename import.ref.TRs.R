

require(data.table)
#require(tidyverse)
require(dplyr)
require(data.table)
require(rtracklayer)


#### Import Toram et al annotation transferred from Kaplan
#viral.ref <- "NC_001491.2.mRNA_corrected.gff3"
viral.ref <- data.table(as.data.frame(rtracklayer::import.gff(viral.ref, version='3')))
#if(nrow(viral.ref[ID != transcript_id]) != 0 ) { stop('The IDs do not match for these transcripts: ', viral.ref[ID != transcript_id])}
viral.ref[type=='exon', exon_number := 1:.N, by = .(transcript_id)]
viral.mrna <- viral.ref[,.(seqnames, source, type, phase, strand, start, end, transcript_id, Name, ID, gene_id, ORF_id)]

export.gff2(viral.ref, "NC_001491.2.mRNA_corrected.gff2")
export.gff (viral.ref, "NC_001491.2.mRNA_corrected.gtf")

## genes of ref transcripts
TR.genes <- unique(viral.mrna[,.(seqnames, strand, gene_id, ORF_id, transcript_id)]); tr.dup <- dup(TR.genes$transcript_id)
if(length(tr.dup != 0)) {message(
  'These Reference Transcripts have some problem with them: ', paste(dup(TR.genes$transcript_id), collapse='; '))
}

## add CDS annotation as gene (because transcripts are annotated as genes in this one)
viral.CDS <- feature.df
viral.ref$ORF_id[!viral.ref$ORF_id %in% viral.CDS$gene]
viral.ref$gene_id[!viral.ref$gene_id %in% viral.CDS$gene]

viral.ref <- data.table(plyr::rbind.fill(viral.mrna, viral.CDS))
export.gff2(viral.ref, "NC_001491.2.mRNA_corrected_w_CDS.genes.gff2")
export.gff3(viral.ref, "NC_001491.2.mRNA_corrected_w_CDS.genes.gff3")

#### CDS instead of genes
viral.CDS.ref <- viral.ref
viral.CDS.ref[type == 'gene', type := 'CDS']
export.gff(viral.CDS.ref, "NC_001491.2.mRNA_corrected_w_CDS.gtf")

#### Annotated transcripts reference
TR.ref <- viral.mrna
TR.ref[,transcript_start:=min(start), by=.(transcript_id)][,transcript_end:=max(end), by=.(transcript_id)]
TR.ref[type == 'mRNA', type:='transcript']
## Exclude duplicate transcripts
dups   <- TR.ref[type != 'exon'][duplicated(transcript_id)][,Name]


#### Merged transcripts and exons reference table
TR.ref <- viral.mrna
TR.ref[,transcript_start:=min(start), by=.(transcript_id)][,transcript_end:=max(end), by=.(transcript_id)]
TR.ref[type == 'mRNA', type:='transcript']
TR.ref[type == 'exon', exon_number := 1:.N, by = .(transcript_id)]

TR.merged.data    <- TR.ref
TR.merged.data.TR <- TR.merged.data[type=='transcript', .(seqnames, strand, transcript_id, start, end)]
setnames(TR.merged.data.TR, old=c('start', 'end'), new=c('start.TR', 'end.TR') )

TR.merged.data.ex <- TR.merged.data[type=='exon',       .(seqnames, strand, transcript_id, start, end, exon_number)]
setnames(TR.merged.data.ex, old=c('start', 'end'), new=c('start.exon', 'end.exon') )

TR.merged.data <- merge(TR.merged.data.TR, TR.merged.data.ex, by=c('seqnames', 'strand', 'transcript_id') )
TR.merged.data[,start := start.TR][,end := end.TR]

TR.merged.data[,prime5.TR := ifelse(strand == '+', start.TR,   end.TR)]
TR.merged.data[,prime5.ex := ifelse(strand == '+', start.exon, end.exon)]

TR.merged.data[,prime3.TR := ifelse(strand == '+', end.TR, start.TR)]
TR.merged.data[,prime3.ex := ifelse(strand == '+', end.exon, start.exon)]

TR.merged.data[,last_exon := ifelse(prime3.TR == prime3.ex, T, F)]

TR.merged.data[,strand := factor(strand, levels=c('+', '-', '*'))]

## Source ????
TR.merged.data[,source:=ifelse(grepl('MSTRG', transcript_id), 'StringTie', 'Torma_et_al')]


TR.Ref.data   <- TR.merged.data


###
colsby <- colnames(TR.merged.data.ex)
TR.ref.ex <- TR.merged.data.ex[,.(exon_pos = paste0(start.exon, '-', end.exon)), 
                               by=colsby]
TR.ref.ex[,exon_composition := paste(exon_pos, collapse = ';'), by=transcript_id]
TR.ref.ex <- unique(TR.ref.ex[,.(seqnames, strand, transcript_id, exon_composition)])

TR.ref.ex[,exon_composition_freq := .N, by=exon_composition]
stopifnot(length(unique(TR.ref.ex$exon_composition)) == nrow(TR.ref.ex))


### DONT RUN
dontrun <- T
if (dontrun) {
  return()
} else {
  
  TR.ref <- fread('Torma_NC_001491.2.tsv', na.strings = '') #("NC_001491.2.mRNA_corrected.tsv")
  colnames(TR.ref)[] <- c('seqnames', 'source', 'type', 'start', 'end', 'phase', 'strand', 'score', 'transcript_id', 'tr_category', 'gene_id', 'ORF_id')
  TR.ref[,source := 'LoRTIA']
  TR.ref[,transcript_id := gsub('transcript_id=', '', transcript_id)]
  TR.ref[,Name := transcript_id]
  TR.ref[,ID := transcript_id]
  export.gff3(TR.ref, "NC_001491.2.mRNA_corrected.gff3")
  
  ## previoous
  TR.ref <- fread("NC_001491.2.mRNA_corrected.tsv", na.strings = '')
  export.gff3(TR.ref, "NC_001491.2.mRNA_corrected.gff3")
}