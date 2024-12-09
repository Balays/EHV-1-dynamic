
#install.packages('data.table')
library(data.table)

countData <- fread('LoRTIA_virus/TSS_abund.norm_LoRTIA/viral_read.count_TSS_abund.cast.tsv')

#install.packages('keras')

library(reticulate)
use_condaenv("keras", required = TRUE) # Replace "keras" with your environment name
py_config()

library(keras)
library(tensorflow)
tf_config()

library(cluster)

# Preprocess your count data (normalize, log-transform if necessary)
# Define and train the autoencoder
encoder <- keras_model_sequential() %>%
  layer_dense(units = 128, activation = "relu", input_shape = ncol(countData)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 10, activation = "relu") # Latent space with 10 dimensions

# Use GPU for faster computation
install_keras(tensorflow = "gpu")

# After training, extract the encoded data
encoded_data <- predict(encoder, normalized_count_data)

# Cluster using K-means or other methods
kmeans_result <- kmeans(encoded_data, centers = 5)
