
melted_z_scores <- melt(z_scores_data, variable.name = 'hpi', id.vars = 'gene')
melted_z_scores <- melted_z_scores[grepl('z_', hpi), ]
melted_z_scores[,hpi := gsub('z_', '', hpi)]
melted_z_scores[,hpi := factor(hpi, levels = time_points)]

# Plot the heatmap using ggplot2
p <- ggplot(melted_z_scores, aes(x = hpi, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(#low = "#d73027", mid = "white", high = "#4575b4",
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = 0, #median(melted_dists$distance, na.rm = TRUE),
    limit    = c(-1, 1), #range(melted_z_scores$value, na.rm = TRUE),
    name="z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title  = element_blank(),
        panel.spacing = unit(0.5, 'mm'))

print(p)


# Plot the heatmap using ggplot2
p <- ggplot(viral_tss_data, aes(x = hpi, y = gene, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(#low = "#d73027", mid = "white", high = "#4575b4",
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = range(viral_tss_data$mean, na.rm = TRUE) / 2,
    limit    = c(0, 1), #range(melted_z_scores$value, na.rm = TRUE),
    name="z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title  = element_blank(),
        panel.spacing = unit(0.5, 'mm'))

print(p)



library(DESeq2)


orf.sum  <- TR.ref.sum

#### Keep only those reads, which could be assigned to a gene (based on '=' ref transcripts or 5-prime overlaps)
orf.sum <- orf.perc[!is.na(gene),]

fwrite(orf.sum, 'TR.orf.sum.tsv', sep = '\t')


#### Summarize the counts on Gene
orf.sum <- orf.sum[,.(read_count  = sum(read_count)),
         by=.(seqnames, gene, sample, group, hpi, Time, cell_line)]

orf.sum.sp     <- dcast.data.table(orf.sum, gene~sample, value.var = 'read_count')

orf.sum.sp.mat <- as.matrix(data.frame(orf.sum.sp[,-1], row.names = orf.sum.sp$gene))
colnames(orf.sum.sp.mat) <- gsub('\\.', '-', colnames(orf.sum.sp.mat))

colData = data.frame(metafilt, row.names = metafilt$sample)

orf.sum.sp.mat <- orf.sum.sp.mat[,rownames(colData)]

deseq <- DESeqDataSetFromMatrix(orf.sum.sp.mat, colData = colData, ~ hpi)
deseq <- estimateSizeFactors(deseq, type =  "poscounts")
deseq <- DESeq(deseq)
results(deseq)

vsd   <- varianceStabilizingTransformation(deseq)
vsd_data <- assay(vsd)
vsd_data <- data.table(gene = rownames(vsd_data), vsd_data)
vsd_melt <- melt.data.table(vsd_data, variable.name = 'sample')
vsd_melt <- merge(vsd_melt, metafilt, by='sample')

rlog      <- rlog(deseq)
rlog_data <- assay(rlog)
rlog_data <- data.table(gene = rownames(rlog_data), rlog_data)
rlog_melt <- melt.data.table(rlog_data, variable.name = 'sample')
rlog_melt <- merge(rlog_melt, metafilt, by='sample')


# Plot the heatmap using ggplot2
p <- ggplot(vsd_melt, aes(x = sample, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(#low = "#d73027", mid = "white", high = "#4575b4",
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = mean(range(vsd_melt$value, na.rm = TRUE)),
    limit    = range(vsd_melt$value, na.rm = TRUE),
    name="VST-normalized expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title  = element_blank(),
        panel.spacing = unit(0.5, 'mm')) +
  facet_nested(cols=vars(hpi), scales = 'free')

ggsave('TR.ORF.Heatmap.jpg', p, height = 12, width = 16)
print(p)



# Plot the heatmap using ggplot2
p <- ggplot(rlog_melt, aes(x = sample, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(#low = "#d73027", mid = "white", high = "#4575b4",
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = mean(range(rlog_melt$value, na.rm = TRUE)),
    limit    = range(rlog_melt$value, na.rm = TRUE),
    name="rlog-normalized expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title  = element_blank(),
        panel.spacing = unit(0.5, 'mm')) +
  facet_nested(cols=vars(hpi), scales = 'free')


print(p)



#####
# Plot the heatmap using ggplot2
p <- ggplot(orf.perc.mean.melt, aes(x = variable, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(#low = "#d73027", mid = "white", high = "#4575b4",
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = mean(range(orf.perc.mean.melt$value, na.rm = TRUE)),
    limit    = range(orf.perc.mean.melt$value, na.rm = TRUE),
    name="Viral read count normalized expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title  = element_blank(),
        panel.spacing = unit(0.5, 'mm')) +
  facet_nested(cols=vars(variable), scales = 'free')

#ggsave('TR.ORF.Heatmap.jpg', p, height = 12, width = 16)
print(p)
