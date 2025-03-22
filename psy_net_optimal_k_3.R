# --- Configuración de computación paralela ---
library(doParallel)
library(dplyr)
library(ggplot2)
library(bootnet)
library(cluster) # For silhouette calculation
library(viridis) # For color scales
library(NetworkComparisonTest)
library(DiscreteGapStatistic)
library(NbClust)
library(qgraph)
library(igraph)

base_igi_bin <- read.csv("psy_net_recidivism_files/base_igi_bin.csv")
descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR24","PAR25","PAR26","PAR27",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

# Define parameters
window_size <- 10 # Size of the sliding window
min_cluster_size <- 50 # Minimum cluster size for stable network estimation
k_values <- 2:3 # Range of cluster numbers to explore
external_variables <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# --- Función para procesar una ventana ---
process_window_parallel <- function(base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values) {
  
  # Parallelized loop over sliding windows
  resultados <- foreach(start_val = 0:(max(base_igi_bin$puntaje_total) - window_size + 1), 
                        .combine = "rbind",
                        .packages = c("dplyr", "bootnet", "cluster"),
                        .export = c("window_size", "min_cluster_size", "k_values", "base_igi_bin", "descriptivo_grupal")) %dopar% {
                          
                          tryCatch({
                            
                            # Define the current window
                            end_val <- start_val + window_size - 1
                            rango_puntaje <- start_val:end_val
                            
                            # Filter data for the current window
                            sub_bloque <- base_igi_bin %>%
                              filter(puntaje_total %in% rango_puntaje)
                            
                            # Skip if there are insufficient observations
                            if (nrow(sub_bloque) < min_cluster_size) {
                              return(data.frame(start_val = start_val, error_message = "Insufficient observations"))
                            }
                            
                            # Loop over candidate k values
                            results_k <- lapply(k_values, function(k) {
                              
                              message(paste("Processing start_val:", start_val, "k_value:", k))
                              
                              # Perform clustering with k clusters
                              matriz_binaria <- as.matrix(sub_bloque[, descriptivo_grupal])
                              distancias <- stats::dist(matriz_binaria, method = "binary")
                              clust_hier <- stats::hclust(distancias, method = "ward.D2")
                              sub_bloque$sub_cluster <- cutree(clust_hier, k = k)
                              
                              # Skip if any cluster has fewer than the minimum required observations
                              cluster_sizes <- table(sub_bloque$sub_cluster)
                              if (any(cluster_sizes < min_cluster_size)) {
                                return(data.frame(start_val = start_val, k = k, error_message = "Cluster too small"))
                              }
                              
                              # Evaluate cluster quality (e.g., silhouette score)
                              silhouette_score <- mean(cluster::silhouette(
                                sub_bloque$sub_cluster, distancias)[, "sil_width"])
                              
                              # Estimate Ising networks for each cluster and compare them
                              network_models <- list()
                              for (i in unique(sub_bloque$sub_cluster)) {
                                cluster_data <- sub_bloque[sub_bloque$sub_cluster == i, descriptivo_grupal]
                                network_models[[i]] <- bootnet::estimateNetwork(cluster_data, default = "IsingFit", tuning = 0.25)
                              }
                              
                              # Compare networks using Pearson correlation (simplified metric)
                              pairwise_comparisons <- combn(unique(sub_bloque$sub_cluster), 2, simplify = FALSE)
                              pears_cor <- sapply(pairwise_comparisons, function(pair) {
                                model_1 <- network_models[[pair[1]]]
                                model_2 <- network_models[[pair[2]]]
                                
                                cor(as.vector(model_1$graph), as.vector(model_2$graph), use = "complete.obs")
                              })
                              
                              # Calculate composite score (lower is better)
                              composite_score <- (1 - silhouette_score) + mean(pears_cor)
                              
                              return(data.frame(
                                start_val = start_val,
                                end_val = end_val,
                                k = k,
                                silhouette_score = silhouette_score,
                                composite_score = composite_score,
                                cluster_sizes = paste(cluster_sizes, collapse = ", ")
                              ))
                            })
                            
                            do.call(rbind, results_k)
                            
                          }, error = function(e) {
                            message(paste("Error in window:", start_val))
                            return(data.frame(start_val = start_val, error_message = e$message))
                          })
                        }
  
  resultados
}


# Set up parallel backend for M1/M2 (4 performance cores + 4 efficiency)
cl <- makeCluster(6, type = "FORK")  # Uses 6 cores (4P + 4E)
registerDoParallel(cl)

# --- Ejecución del análisis paralelo ---
time_init <- Sys.time()

results_grid_parallel <- process_window_parallel(base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)

time_fin <- Sys.time()
time_total_parallel <- difftime(time_fin, time_init, units = "auto")  # Time in seconds

# --- Limpieza y almacenamiento de resultados ---
stopCluster(cl)

write.csv(results_grid_parallel, "results_grid_parallel.csv")

# Visualize results with a heatmap
ggplot(results_grid_parallel, aes(x = start_val, y = k, fill = composite_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Composite Score",
    title = "Optimal Parameters Heatmap"
  )
