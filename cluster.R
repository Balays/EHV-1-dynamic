# Example: correlation-based distance
dist_matrix <- as.dist(1 - cor(t(vsd_mat[ , -1]), method="pearson"))
hc <- hclust(dist_matrix, method="complete")
plot(hc)
# Cut the tree where it makes sense biologically
clusters <- cutree(hc, k = 5) # for example, 5 clusters


vsd_mat <- vsd_mat[,metafilt$sample]

# Extract time points from column names
col_names <- colnames(vsd_mat)
# Assuming format like "EHV-1_2h_1", we can parse out the hour:
hours <- sapply(col_names, function(x) str_extract(x, "\\d+h"))
time_points <- as.numeric(str_replace(hours, "h", ""))

# Average replicates at each time point to get a single expression value per time point per gene
unique_times <- sort(unique(time_points))
avg_mat <- sapply(unique_times, function(tp) {
  rowMeans(vsd_mat[, time_points == tp, drop=FALSE])
})

rownames(avg_mat) <- rownames(vsd_mat)
colnames(avg_mat) <- paste0("Time_", unique_times, "h")


library(dtwclust)

# Perform DTW-based clustering
# k = number of clusters (example: 5)
set.seed(123)
res <- tsclust(avg_mat,
               type = "partitional",
               k = 5,
               distance = "dtw_basic",
               centroid = "pam",
               seed = 123)

# Inspect results
print(res)
# res@cluster gives cluster assignments for each gene
cluster_assignments <- res@cluster

# Plot the clusters
plot(res)




library(GRENITS)

# GRENITS requires a matrix: genes in rows, time points in columns
TimeSeries <- avg_mat  # already genes x time
# Ensure TimeSeries is numeric and log/vst transformed; should be from DESeq2 VST.

# Run GRENITS time-series inference
# Adjust nIterations, burnin, and thinning as needed.
# Larger nIterations = more accurate but slower.
results <- ReplicatesNet_gauss(
  data = TimeSeries,
  resultsFolder = "GRENITS_output",
  nIterations = 50000,
  burnin = 10000,
  thinning = 10
)

# After this finishes, analyze results:
plotResults("GRENITS_output", what="network")

# Extract the inferred network (average.over = TRUE gives averaged adjacency)
network <- extractNetwork("GRENITS_output", average.over=TRUE)
head(network)
