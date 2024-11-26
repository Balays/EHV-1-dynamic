CustomClustering <- function(DT,
                             max_distance_within_cluster = 10,
                             max_cluster_size = Inf,
                             sep_strands = TRUE,
                             include_zeros = FALSE,
                             plot_results = TRUE) {
  # Ensure only positions with counts > 0 are included
  if (!include_zeros) {
    DT <- DT[count > 0, ]
    message('Positions with a count > 0 will be used only.')
  }

  if (sep_strands) {
    DT_list <- split(DT, DT$strand)
  } else {
    DT_list <- list(DT)
  }

  result_list <- lapply(DT_list, function(DT_subset) {
    if (nrow(DT_subset) == 0) return(NULL)

    # Sort positions
    DT_subset <- DT_subset[order(position)]

    # Initialize variables
    cluster_id <- 1
    point_in_cluster <- 1  # Counter for points in the current cluster

    DT_subset[, cluster_ID := {
      cluster_ids <- integer(.N)
      cluster_ids[1] <- cluster_id
      for (i in 2:.N) {
        distance <- position[i] - position[i - 1]
        # Check if adding the point exceeds max_distance_within_cluster or max_cluster_size
        if (distance > max_distance_within_cluster || point_in_cluster >= max_cluster_size) {
          cluster_id <<- cluster_id + 1
          point_in_cluster <<- 1  # Reset point counter
        } else {
          point_in_cluster <<- point_in_cluster + 1
        }
        cluster_ids[i] <- cluster_id
      }
      cluster_ids
    }]
    DT_subset[, cluster_ID := paste0('cluster_', cluster_ID)]

    # Calculate cluster statistics
    DT_subset[, `:=`(
      cluster_start  = min(position),
      cluster_end    = max(position),
      cluster_width  = max(position) - min(position) + 1,
      cluster_center = round(mean(position)),
      cluster_peak   = position[which.max(count)]
    ), by = cluster_ID]

    return(DT_subset)
  })

  # Combine results
  DT_result <- rbindlist(result_list)

  # Plot results if requested
  if (plot_results) {
    ggplot(DT_result, aes(x = position, y = count, color = cluster_ID)) +
      geom_point(size = 2) +
      facet_wrap(~ strand) +
      theme_minimal() +
      labs(title = 'Custom Clustering Results with Max Cluster Size',
           x = 'Position',
           y = 'Count',
           color = 'Cluster ID')
  }

  return(DT_result)
}

if(!dontrun) {
  DT_clust_custom <- CustomClustering(
    DT, # = tx_uni_trs_nov[position >= 1908 & position <= 2012 & strand == '-'],
    max_distance_within_cluster = 5,  # Adjust as needed
    sep_strands = FALSE,
    include_zeros = FALSE,
    plot_results = TRUE
  )
}
