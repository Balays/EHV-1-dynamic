library(mclust)
library(data.table)
library(ggplot2)

Mclust_genes <- function(countData, use_weights = FALSE, plot_results = FALSE, G = NULL) {
  # countData: data.table with rows as genes and columns as samples
  # use_weights: if TRUE, include additional weighting for clustering
  # plot_results: if TRUE, plot the clustering results
  # G: Number of mixture components (clusters). If NULL, let Mclust select the best number.

  # Check for required structure
  if (!is.data.table(countData)) {
    stop("countData must be a data.table with genes as rows and samples as columns.")
  }

  # Transpose the data to cluster genes (rows become columns)
  gene_ids <- countData$gene
  countData <- as.matrix(countData[, -1, with = FALSE])  # Exclude gene IDs for clustering

  if (!use_weights) {
    message("Weights will not be used; all genes will have equal importance.")
    weights <- NULL
  } else {
    # Example weights: sum of counts across all samples
    weights <- rowSums(countData)
    message("Using total expression as weights for clustering.")
  }

  # Perform clustering with Mclust
  mclust_result <- Mclust(data = countData, G = G, weights = weights)

  # Add clustering results to gene data
  cluster_ids <- mclust_result$classification
  result_dt <- data.table(gene = gene_ids, cluster = cluster_ids)

  # Optionally plot the clustering results (PCA projection)
  if (plot_results) {
    library(ggplot2)
    pca_result <- prcomp(countData, scale. = TRUE)
    pca_dt <- data.table(pca_result$x[, 1:2], cluster = factor(cluster_ids), gene = gene_ids)
    ggplot(pca_dt, aes(x = PC1, y = PC2, color = cluster, label = gene)) +
      geom_point(size = 2) +
      geom_text(hjust = 0.5, vjust = -0.5) +
      theme_minimal() +
      labs(title = "Gene Clustering (Mclust)", x = "PC1", y = "PC2", color = "Cluster")
  }

  return(list(clustering_result = mclust_result, result_table = result_dt))
}

# Load your data

countData <- fread('LoRTIA_virus/TSS_abund.norm_LoRTIA/viral_read.count_TSS_abund.cast.tsv')[,-c(2,3)]

#countData <- fread("viral_read.count_TSS_abund.cast.tsv")

# Run the clustering function
mclust_results <- Mclust_genes(
  countData = countData,
  use_weights = TRUE,
  plot_results = TRUE,
  G = NULL  # Let Mclust determine the optimal number of clusters
)

# Access results
clustering_result <- mclust_results$clustering_result
result_table <- mclust_results$result_table
