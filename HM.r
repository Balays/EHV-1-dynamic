library(DESeq2)


orf.sum  <- TR.ref.sum

#### Keep only those reads, which could be assigned to a gene (based on '=' ref transcripts or 5-prime overlaps)
orf.sum <- orf.perc[!is.na(gene),]

fwrite(orf.sum, 'TR.orf.sum.tsv', sep = '\t')


orf.sum <- fread(file.path(outdir, 'TR.orf.sum.tsv'), na.strings = '')

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


kin_classes_and_clusters <- merge(kin_classes, cluster_labels, by.x='gene', by.y='ORF', all=T)
kin_classes_and_clusters[,gene := gsub('_.*', '', gene)]

vsd_melt <- merge(kin_classes_and_clusters, vsd_melt, by='gene', all.y=T)
kin_classes_and_clusters
vsd_melt <- unique(vsd_melt)

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

# Ensure the order of genes is preserved based on vsd_melt data
vsd_melt$gene <- factor(vsd_melt$gene, levels = unique(vsd_melt$gene))

# Create a ggplot for the kinetic_class and best_label annotations
annotation_data <- unique(vsd_melt[, .(gene, Kinetic_class, best_label)])

# Plot annotations for kinetic class
p_kinetic_class <- ggplot(annotation_data, aes(x = 1, y = gene, fill = Kinetic_class)) +
  geom_tile() +
  scale_fill_manual(values = c("L" = "skyblue", "E" = "lightgreen", "IE" = "yellow", "unknown" = "grey"),
                    na.value = "grey", guide = guide_legend(ncol = 1)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.box   = "vertical",
        #axis.text.y  = element_text(angle = 90, hjust = 1, vjust = 0.5), #element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
        #,plot.margin  = unit(c(5,0,5,15), 'mm')
        ) +
  labs(fill = "Kinetic Class")

# Plot annotations for best_label
p_best_label <- ggplot(annotation_data, aes(x = 1, y = gene, fill = as.factor(best_label))) +
  geom_tile() +
  scale_fill_brewer(palette = "Set1",
                    na.value = "grey", guide = guide_legend(ncol = 1)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.box   = "vertical",
        axis.text.y  = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
        #,plot.margin  = unit(c(5,0,5,0), 'mm')
        ) +
  labs(fill = "Best Label")

# Plot the heatmap for gene expression
p_heatmap <- ggplot(vsd_melt, aes(x = sample, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4",
                       midpoint = mean(range(vsd_melt$value, na.rm = TRUE)),
                       limit = range(vsd_melt$value, na.rm = TRUE),
                       name = "VST-normalized expression") +
  theme_void() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(), # Hide y-axis text
        axis.title  = element_blank(),
        panel.spacing = unit(0.5, 'mm'),
        #,plot.margin   = unit(c(5,0,5,5), 'mm')
        ) +
  facet_nested(cols = vars(hpi), scales = 'free')


combined_plot <- plot_grid(
  p_kinetic_class, p_best_label, p_heatmap,
  ncol = 3, align = 'vh', axis = 'tlb', rel_widths = c(0.1, 0.1, 1)
)

# Save the combined plot
ggsave('TR.ORF.Heatmap_Annotated.jpg', combined_plot, height = 12, width = 18)

print(combined_plot)

#### complexhaeatmap
library(ComplexHeatmap); library(circlize); library(RColorBrewer); library(viridis)
# Convert vsd_melt to a matrix format suitable for ComplexHeatmap
vsd_matrix <- as.data.frame(dcast(vsd_melt, gene ~ sample, value.var = "value"))
rownames(vsd_matrix) <- vsd_matrix$gene

# Ordering
vsd_matrix <- vsd_matrix[,metafilt$sample]

# Convert vsd_matrix to a numeric matrix
vsd_matrix <- as.matrix(vsd_matrix)

# Create row annotations for Kinetic Class and Best Label
row_annotation <- rowAnnotation(
  Kinetic_class = vsd_melt$Kinetic_class[match(rownames(vsd_matrix), vsd_melt$gene)],
  Best_label = vsd_melt$best_label[match(rownames(vsd_matrix), vsd_melt$gene)],
  col = list(Kinetic_class = c("L" = "skyblue", "E" = "lightgreen", "IE" = "yellow", "unknown" = "grey"),
             Best_label = colorRamp2(c(1, 2, 3, 4, 5), brewer.pal(5, "Set1"))),
  na_col = "grey"
)

# Create a color scale for the top annotation (hpi) using viridis inferno
hpi_levels <- unique(vsd_melt$hpi)
hpi_colors <- inferno(length(hpi_levels))

top_annotation <- HeatmapAnnotation(
  hpi = anno_simple(vsd_melt$hpi[1:ncol(vsd_matrix)],
                    col = setNames(hpi_colors, hpi_levels)),
  annotation_name_side = "left"
)

# Create the heatmap with row ordering but without column ordering
heatmap <- Heatmap(
  vsd_matrix,
  name = "VST-normalized expression",
  row_order = NULL, # Allow automatic clustering for row ordering
  column_order = 1:ncol(vsd_matrix), # Keep the original column order (no reordering)
  top_annotation = HeatmapAnnotation(hpi = vsd_melt$hpi[1:ncol(vsd_matrix)]), #top_annotation,
  left_annotation = row_annotation,
  show_row_names = TRUE,
  show_column_names = FALSE,
  clustering_distance_rows = "euclidean", # Clustering by similarity (e.g., Euclidean distance)
  clustering_method_rows = "complete",
  cluster_columns = FALSE, # Do not cluster columns
  col = colorRamp2(c(min(vsd_matrix, na.rm = TRUE),
                     mean(range(vsd_matrix, na.rm = TRUE)),
                     max(vsd_matrix, na.rm = TRUE)),
                   c("#d73027", "white", "#4575b4"))
)

# Draw the heatmap
draw(heatmap)





######
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
