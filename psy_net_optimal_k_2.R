# Load required libraries
library(dplyr)
library(bootnet)
library(cluster) # For silhouette calculation
library(viridis) # For color scales
library(future.apply) # For parallel processing
#install.packages("NetworkComparisonTest")
library(NetworkComparisonTest)
#install.packages("DiscreteGapStatistic")
library(DiscreteGapStatistic)
#install.packages("NbClust")
library(NbClust)
library(qgraph)
library(igraph)


# Initialize result containers
cor_results_sliding <- data.frame()
chi_results_sliding <- data.frame()
cross_tables_sliding <- list()


# ---------------------


# Load required libraries
library(dplyr)
library(ggplot2)




# Function to process a single window and evaluate all k values
process_window <- function(start_val) {
  
  # Define the current window
  end_val <- start_val + window_size - 1
  rango_puntaje <- start_val:end_val
  
  # Filter data for the current window
  sub_bloque <- base_igi_bin %>%
    filter(puntaje_total %in% rango_puntaje)
  
  # Skip if there are insufficient observations
  if (nrow(sub_bloque) < min_cluster_size) return(NULL)
  
  # Loop over candidate k values (parallelized)
  results_k <- future_lapply(k_values, function(k) {
    
    print("start_val: ", start_val, "k_value", k )
    
    # Perform clustering with k clusters
    matriz_binaria <- as.matrix(sub_bloque[, descriptivo_grupal])
    distancias <- stats::dist(matriz_binaria, method = "binary")
    clust_hier <- stats::hclust(distancias, method = "ward.D2")
    sub_bloque$sub_cluster <- cutree(clust_hier, k = k)
    
    # Skip if any cluster has fewer than the minimum required observations
    cluster_sizes <- table(sub_bloque$sub_cluster)
    if (any(cluster_sizes < min_cluster_size)) return(NULL)
    
    # Evaluate cluster quality (e.g., silhouette score)
    silhouette_score <- mean(cluster::silhouette(
      sub_bloque$sub_cluster, distancias)[, "sil_width"])
    
    # Estimate Ising networks for each cluster and compare them
    network_models <- list()
    for (i in unique(sub_bloque$sub_cluster)) {
      cluster_data <- sub_bloque[sub_bloque$sub_cluster == i, descriptivo_grupal]
      network_models[[i]] <- bootnet::estimateNetwork(cluster_data, default = "IsingFit", tuning = 0.25)
    }
    
    # Compare networks using Network Comparison Test (NCT)
    pairwise_comparisons <- combn(unique(sub_bloque$sub_cluster), 2, simplify = FALSE)
    nct_pvals <- sapply(pairwise_comparisons, function(pair) {
      model_1 <- network_models[[pair[1]]]
      model_2 <- network_models[[pair[2]]]
      NetworkComparisonTest::NCT(model_1, model_2, it = 1000)$nwinv.pval
    })
    
    avg_nct_pval <- mean(nct_pvals) # Average p-value across comparisons
    
    # Calculate composite score (lower is better)
    composite_score <- (1 - silhouette_score) + avg_nct_pval
    
    return(data.frame(
      start_val = start_val,
      end_val = end_val,
      k = k,
      silhouette_score = silhouette_score,
      avg_nct_pval = avg_nct_pval,
      composite_score = composite_score,
      cluster_sizes = paste(cluster_sizes, collapse = ", ")
    ))
  })
  
  # Combine results for all k values into a single data frame
  do.call(rbind, results_k)
}

# Define parameters
window_size <- 10 # Size of the sliding window
min_cluster_size <- 50 # Minimum cluster size for stable network estimation
k_values <- 2:6 # Range of cluster numbers to explore
external_variables <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# Initialize result container
results_grid <- data.frame()

# Set up parallel processing

time_init <- Sys.time()

plan(multisession, workers = availableCores() - 1) # Use all but one core

# Parallelized loop over sliding windows
results_list <- future_lapply(0:(max(base_igi_bin$puntaje_total) - window_size + 1), process_window)

# Combine results into a single data frame
results_grid <- do.call(rbind, results_list)

time_fin <- Sys.time()
time_tot <- difftime(time_fin, time_init, units = "auto")  # Tiempo en segundos

# Step: Visualize results with a heatmap
ggplot(results_grid, aes(x = start_val, y = k, fill = composite_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Composite Score",
    title = "Optimal Parameters Heatmap"
  )

