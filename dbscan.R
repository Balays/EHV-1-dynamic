library(dbscan)
library(data.table)
library(ggplot2)

DBSCAN_genes <- function(countData, eps = 0.5, minPts = 5, plot_results = FALSE) {
  # countData: data.table with rows as genes and columns as samples
  # eps: The epsilon parameter for DBSCAN (neighborhood size)
  # minPts: Minimum number of points required to form a cluster
  # plot_results: If TRUE, plot the clustering results (PCA visualization)

  # Check for required structure
  if (!is.data.table(countData)) {
    stop("countData must be a data.table with genes as rows and samples as columns.")
  }

  # Extract gene IDs and normalize the count data
  gene_ids <- countData$gene
  countData <- as.matrix(countData[, -1, with = FALSE])  # Exclude the gene column
  countData <- scale(countData)  # Scale the data for clustering

  # Run DBSCAN
  dbscan_result <- dbscan(countData, eps = eps, minPts = minPts)

  # Add clustering results to gene data
  result_dt <- data.table(gene = gene_ids, cluster = dbscan_result$cluster)

  # Optionally plot the clustering results (PCA projection)
  if (plot_results) {
    pca_result <- prcomp(countData, scale. = TRUE)
    pca_dt <- data.table(pca_result$x[, 1:2], cluster = factor(dbscan_result$cluster), gene = gene_ids)
    ggplot(pca_dt, aes(x = PC1, y = PC2, color = cluster, label = gene)) +
      geom_point(size = 2) +
      geom_text(hjust = 0.5, vjust = -0.5) +
      theme_minimal() +
      labs(title = "Gene Clustering (DBSCAN)", x = "PC1", y = "PC2", color = "Cluster")
  }

  return(result_dt)
}

# Example Usage:

countData <- fread('LoRTIA_virus/TSS_abund.norm_LoRTIA/viral_read.count_TSS_abund.cast.tsv')[,-c(2,3)]

# Run the DBSCAN clustering function
dbscan_results <- DBSCAN_genes(
  countData = countData,
  eps = 1,  # Adjust epsilon based on your dataset
  minPts = 2,  # Minimum points per cluster
  plot_results = TRUE
)

# View results
dbscan_results[,.N,cluster]
