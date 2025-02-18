

#countData <- vsd_mat#[,metafilt$sample]


# correlation-based distance
dist_matrix <- as.dist(1 - cor(t(countData[ ,]), method="pearson"))
hc <- hclust(dist_matrix, method="complete")
plot(hc)
# Cut the tree where it makes sense biologically
clusters <- cutree(hc, k = cluster_num) # for example, 5 clusters


cluster_dt <- data.table(cluster=clusters, gene=rownames(countData))






