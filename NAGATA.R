
nagata.fwd <- data.table(as.data.frame(rtracklayer::import.gff3('./NAGATA/Final_cluster.NAGATA.fwd.gff3')))

nagata.rev <- data.table(as.data.frame(rtracklayer::import.gff3('./NAGATA/Final_cluster.NAGATA.rev.gff3')))


nagata.gff <- rbind(nagata.fwd, nagata.rev)


nagata.gff <- nagata.gff[type %in% c('exon', 'mRNA'),]

nagata.gff[,Parent := as.character(Parent)]

nagata.gff[type == 'mRNA', Parent := ID]

nagata.gff[,ID := Parent]

nagata.gff[,Parent := NULL][,Name := NULL]


nagata.gff[ID %in% dup(nagata.gff$ID), ID := paste0(ID, '::', strand)]


nagata.mrna <- nagata.gff[type == 'mRNA',]



stopifnot(
  nrow(nagata.mrna) == length(unique(nagata.mrna[,ID]))
)




nagata.ex   <- nagata.gff[type == 'exon',]




rtracklayer::export.gff3(as.data.frame(nagata.gff), 'NAGATA_dRNA.gff3')
