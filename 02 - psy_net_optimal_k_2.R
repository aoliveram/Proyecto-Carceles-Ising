# Barrido paramétrico en ventana puntaje y k, para seleccionar un valor de k óptimo.

# --- Configuración de computación paralela ---

library(dplyr)
library(tidyr)
library(ggplot2)
library(bootnet)
library(cluster) # For silhouette calculation
library(viridis) # For color scales
library(qgraph)
library(igraph)
library(doParallel)

descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR24","PAR25","PAR26","PAR27",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

# --- Función para procesar una ventana ---
process_window_parallel <- function(start_val_in, base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values) {
  
  resultados <- foreach(
    #start_val = 0:(max(base_igi_bin$puntaje_total) - window_size + 1), 
    start_val = start_val_in, 
    .combine = "rbind",
    .packages = c("dplyr", "bootnet", "cluster")
    #.export = c("window_size", "min_cluster_size", "k_values", "base_igi_bin", "descriptivo_grupal") 
    # not necessary if the variables are defined globally and not inside the function
    ) %dopar% {
      
      tryCatch({
        
        # Define the current window
        end_val <- start_val + window_size - 1
        rango_puntaje <- start_val:end_val
        
        # Filter data for the current window
        sub_bloque <- base_igi_bin %>%
          filter(puntaje_total %in% rango_puntaje)
        
        # Skip if there are insufficient observations
        if (nrow(sub_bloque) < min_cluster_size) {
          return(data.frame(
            start_val = start_val, end_val = end_val, k = NA, silhouette_score = NA,
            composite_score = NA, mean_pears_cor = NA, cluster_sizes = NA,
            error_message = "Insufficient observations"))
        }
        
        # Loop over candidate k values
        results_k <- lapply(k_values, function(k) {
          
          #message(paste("Processing start_val:", start_val, "k_value:", k))
          
          # Perform clustering with k clusters
          matriz_binaria <- as.matrix(sub_bloque[, descriptivo_grupal])
          distancias <- stats::dist(matriz_binaria, method = "binary")
          clust_hier <- stats::hclust(distancias, method = "ward.D2")
          sub_bloque$sub_cluster <- cutree(clust_hier, k = k)
          
          # Skip if any cluster has fewer than the minimum required observations
          cluster_sizes <- table(sub_bloque$sub_cluster)
          if (any(cluster_sizes < min_cluster_size)) {
            return(data.frame(
              start_val = start_val, end_val = end_val, k = k, silhouette_score = NA,
              composite_score = NA, mean_pears_cor = NA, cluster_sizes = NA,
              error_message = "Cluster too small"))
          }
          
          # Evaluate cluster quality (silhouette width: close 1: obsv. well clusteres
          # around 0: obsv. lies between two clusters,  negative: placed in the wrong cluster)
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
            mean_pears_cor = mean(pears_cor),
            cluster_sizes = paste(cluster_sizes, collapse = ", ")
          ))
        })
        
        do.call(rbind, results_k)
        
      }, error = function(e) {
        warning(paste("Error in window:", start_val, "-", e$message))
        return(data.frame(
          start_val = start_val, end_val = NA, k = NA,
          silhouette_score = NA, composite_score = NA, mean_pears_cor = NA,
          cluster_sizes = NA, error_message = e$message))
      })
    }
  
  resultados
}

# -------------------- Formato {3,2,1},{0} -> {0},{1} --------------------------

base_igi_bin <- read.csv("psy_net_files/base_igi_bin_1.csv")

# Define parameters
window_size <- 10 # Size of the sliding window
min_cluster_size <- 50 # Minimum cluster size for stable network estimation
k_values <- 2:6 # Range of cluster numbers to explore
external_variables <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

start_val_vec = 0:(max(base_igi_bin$puntaje_total) - window_size + 1)

# See how many obsv we have per sub_bloque
# for (start_val in start_val_vec) {
#   end_val <- start_val + window_size - 1
#   rango_puntaje <- start_val:end_val
#   
#   # Filter data for the current window
#   sub_bloque <- base_igi_bin %>%
#     filter(puntaje_total %in% rango_puntaje)
#   
#   print(paste("start_val:", start_val, ", nrow(sub_bloque):", nrow(sub_bloque)))
# }

# -- BY BLOCKS

cl <- makeCluster(12, type = "FORK")  
registerDoParallel(cl)

time_init <- Sys.time()
results_grid_parallel_1 <- process_window_parallel(start_val_vec[1:10], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel <- difftime(time_fin, time_init, units = "auto") # 3.85 k=2:6

time_init <- Sys.time()
results_grid_parallel_2 <- process_window_parallel(start_val_vec[11:15], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_2 <- difftime(time_fin, time_init, units = "auto") # 5.14 min

time_init <- Sys.time()
results_grid_parallel_3 <- process_window_parallel(start_val_vec[16:20], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_3 <- difftime(time_fin, time_init, units = "auto") # 3.98 min

time_init <- Sys.time()
results_grid_parallel_4 <- process_window_parallel(start_val_vec[21:25], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_4 <- difftime(time_fin, time_init, units = "auto") # 2.3 min

time_init <- Sys.time()
results_grid_parallel_5 <- process_window_parallel(start_val_vec[26:28], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_5 <- difftime(time_fin, time_init, units = "auto") # 1.8 min

time_init <- Sys.time()
results_grid_parallel_6 <- process_window_parallel(start_val_vec[29:30], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_6

stopCluster(cl)

write.csv(results_grid_parallel_1, "psy_net_files/optimal_k_grid_parallel_1.csv", row.names = FALSE)
write.csv(results_grid_parallel_2, "psy_net_files/optimal_k_grid_parallel_2.csv", row.names = FALSE)
write.csv(results_grid_parallel_3, "psy_net_files/optimal_k_grid_parallel_3.csv", row.names = FALSE)
write.csv(results_grid_parallel_4, "psy_net_files/optimal_k_grid_parallel_4.csv", row.names = FALSE)
write.csv(results_grid_parallel_5, "psy_net_files/optimal_k_grid_parallel_5.csv", row.names = FALSE)
write.csv(results_grid_parallel_6, "psy_net_files/optimal_k_grid_parallel_6.csv", row.names = FALSE)

results_grid_parallel_1 <- read.csv("psy_net_files/optimal_k_grid_parallel_1.csv")
results_grid_parallel_2 <- read.csv("psy_net_files/optimal_k_grid_parallel_2.csv")
results_grid_parallel_3 <- read.csv("psy_net_files/optimal_k_grid_parallel_3.csv")
results_grid_parallel_4 <- read.csv("psy_net_files/optimal_k_grid_parallel_4.csv")
results_grid_parallel_5 <- read.csv("psy_net_files/optimal_k_grid_parallel_5.csv")
results_grid_parallel_6 <- read.csv("psy_net_files/optimal_k_grid_parallel_6.csv")

# -- ALL AT ONCE (DONT WORK FINE)

cl <- makeCluster(12, type = "FORK")  # Uses 6 cores (4P + 4E)
registerDoParallel(cl)
time_init <- Sys.time()

results_grid_parallel_all <- process_window_parallel(start_val_vec, base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)

time_fin <- Sys.time()
time_total_parallel_all <- difftime(time_fin, time_init, units = "auto") # 

stopCluster(cl)

write.csv(results_grid_parallel_all, "psy_net_files/optimal_k_grid_parallel_all.csv", row.names = FALSE)

# -- EXPLORE RESULTS

# Combine all data frames into one
results_grid_parallel_all <- rbind(
  results_grid_parallel_1,
  results_grid_parallel_2,
  results_grid_parallel_3,
  results_grid_parallel_4,
  results_grid_parallel_5,
  results_grid_parallel_6
)

png("psy_net_plots/optimal_k_mean_cor_1.png", width = 800*1.1, height = 600*1.1)
ggplot(results_grid_parallel_all, aes(x = start_val, y = k, fill = mean_pears_cor)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Score (lower is better)",
    title = "Mean Pears Cor Heatmap Parallel {3,2,1},{0} -> {0},{1}"
  )
dev.off()

png("psy_net_plots/optimal_k_silhouette_score_1.png", width = 800*1.1, height = 600*1.1)
ggplot(results_grid_parallel_all, aes(x = start_val, y = k, fill = silhouette_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Score (higher is better)",
    title = "Silhouette Score Heatmap Parallel {3,2,1},{0} -> {0},{1}"
  )
dev.off()

png("psy_net_plots/optimal_k_composite_score_1.png", width = 800*1.1, height = 600*1.1)
ggplot(results_grid_parallel_all, aes(x = start_val, y = k, fill = composite_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Score (lower is better)",
    title = "Composite Score Heatmap Parallel {3,2,1},{0} -> {0},{1}"
  )
dev.off()

# Comparamos con el original

cor_results_sliding <- read.csv("psy_net_files/cor_results_sliding_1.csv")
cor_results_sliding$mean_cor <- rowMeans(cor_results_sliding[, c("cor_1_2", "cor_1_3", "cor_2_3")])

# Creamos objeto para plot
cor_long <- cor_results_sliding %>%
  select(min_puntaje, max_puntaje, cor_1_2, cor_1_3, cor_2_3) %>%
  pivot_longer(
    cols = c("cor_1_2", "cor_1_3", "cor_2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )
ggplot(cor_long, aes(x = min_puntaje, y = correlacion)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ clusters_comparados, ncol = 1) +
  theme_minimal() + 
  ylim(0,1) +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Correlación de la red",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  )

png("psy_net_plots/optimal_k_cor_original.png", width = 800*1.1, height = 600*1.1/5)
ggplot(cor_results_sliding, aes(x = min_puntaje, y = 3, fill = mean_cor)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Number of Clusters (k)",
    fill = "Cor (lower is better)",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  )
dev.off()


# -------------------- Formato {3,2},{0} -> {0,1},{1} --------------------------

base_igi_bin <- read.csv("psy_net_files/base_igi_bin_2.csv")

# Define parameters
window_size <- 10 # Size of the sliding window
min_cluster_size <- 50 # Minimum cluster size for stable network estimation
k_values <- 2:6 # Range of cluster numbers to explore
external_variables <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

start_val_vec = 0:(max(base_igi_bin$puntaje_total) - window_size + 1)

# -- BY BLOCKS

cl <- makeCluster(12, type = "FORK")  
registerDoParallel(cl)

time_init <- Sys.time()
results_grid_parallel_1 <- process_window_parallel(start_val_vec[2:10], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_1 <- difftime(time_fin, time_init, units = "auto") # 2.76

time_init <- Sys.time()
results_grid_parallel_2 <- process_window_parallel(start_val_vec[11:15], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_2 <- difftime(time_fin, time_init, units = "auto") # 3.01 min

time_init <- Sys.time()
results_grid_parallel_3 <- process_window_parallel(start_val_vec[16:20], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_3 <- difftime(time_fin, time_init, units = "auto") # 3.91 min

time_init <- Sys.time()
results_grid_parallel_4 <- process_window_parallel(start_val_vec[21:25], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_4 <- difftime(time_fin, time_init, units = "auto") # 3.50 min

time_init <- Sys.time()
results_grid_parallel_5 <- process_window_parallel(start_val_vec[26:30], base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_5 <- difftime(time_fin, time_init, units = "auto") # 2.84 min

stopCluster(cl)

write.csv(results_grid_parallel_1, "psy_net_files/optimal_k_grid_parallel_2_1.csv", row.names = FALSE)
write.csv(results_grid_parallel_2, "psy_net_files/optimal_k_grid_parallel_2_2.csv", row.names = FALSE)
write.csv(results_grid_parallel_3, "psy_net_files/optimal_k_grid_parallel_2_3.csv", row.names = FALSE)
write.csv(results_grid_parallel_4, "psy_net_files/optimal_k_grid_parallel_2_4.csv", row.names = FALSE)
write.csv(results_grid_parallel_5, "psy_net_files/optimal_k_grid_parallel_2_5.csv", row.names = FALSE)

results_grid_parallel_2_1 <- read.csv("psy_net_files/optimal_k_grid_parallel_2_1.csv")
results_grid_parallel_2_2 <- read.csv("psy_net_files/optimal_k_grid_parallel_2_2.csv")
results_grid_parallel_2_3 <- read.csv("psy_net_files/optimal_k_grid_parallel_2_3.csv")
results_grid_parallel_2_4 <- read.csv("psy_net_files/optimal_k_grid_parallel_2_4.csv")
results_grid_parallel_2_5 <- read.csv("psy_net_files/optimal_k_grid_parallel_2_5.csv")


# -- ALL AT ONCE (DONT WORK FINE)

base_igi_bin <- read.csv("psy_net_files/base_igi_bin_2.csv")
window_size <- 10 # Size of the sliding window
min_cluster_size <- 50 # Minimum cluster size for stable network estimation
k_values <- 2:6 # Range of cluster numbers to explore
external_variables <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")
start_val_vec = 0:(max(base_igi_bin$puntaje_total) - window_size + 1)

cl <- makeCluster(12, type = "FORK") 
registerDoParallel(cl)
time_init <- Sys.time()
results_grid_parallel_all <- process_window_parallel(start_val_vec, base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)
time_fin <- Sys.time()
time_total_parallel_all <- difftime(time_fin, time_init, units = "auto") # 

stopCluster(cl)

write.csv(results_grid_parallel_all, "psy_net_files/optimal_k_grid_parallel_all_2.csv", row.names = FALSE)

# -- EXPLORE RESULTS 2 

# Combine all data frames into one
results_grid_parallel_all_2 <- rbind(
  results_grid_parallel_2_1,
  results_grid_parallel_2_2,
  results_grid_parallel_2_3,
  results_grid_parallel_2_4,
  results_grid_parallel_2_5
)

png("psy_net_plots/optimal_k_mean_cor_2.png", width = 800*1.1, height = 600*1.1)
ggplot(results_grid_parallel_all_2, aes(x = start_val, y = k, fill = mean_pears_cor)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Score (lower is better)",
    title = "Mean Pears Cor Heatmap Parallel {3,2},{0,1} -> {0},{1}"
  )
dev.off()

png("psy_net_plots/optimal_k_silhouette_score_2.png", width = 800*1.1, height = 600*1.1)
ggplot(results_grid_parallel_all_2, aes(x = start_val, y = k, fill = silhouette_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Score (higher is better)",
    title = "Silhouette Score Heatmap Parallel {3,2},{0,1} -> {0},{1}"
  )
dev.off()

png("psy_net_plots/optimal_k_composite_score_2.png", width = 800*1.1, height = 600*1.1)
ggplot(results_grid_parallel_all_2, aes(x = start_val, y = k, fill = composite_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Score (lower is better)",
    title = "Composite Score Heatmap Parallel {3,2},{0,1} -> {0},{1}"
  )
dev.off()

# Comparamos con el original

cor_results_sliding <- read.csv("psy_net_files/cor_results_sliding_1.csv")
cor_results_sliding$mean_cor <- rowMeans(cor_results_sliding[, c("cor_1_2", "cor_1_3", "cor_2_3")])

# Creamos objeto para plot
cor_long <- cor_results_sliding %>%
  select(min_puntaje, max_puntaje, cor_1_2, cor_1_3, cor_2_3) %>%
  pivot_longer(
    cols = c("cor_1_2", "cor_1_3", "cor_2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )
ggplot(cor_long, aes(x = min_puntaje, y = correlacion)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ clusters_comparados, ncol = 1) +
  theme_minimal() + 
  ylim(0,1) +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Correlación de la red",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  )

png("psy_net_plots/optimal_k_cor_original.png", width = 800*1.1, height = 600*1.1/5)
ggplot(cor_results_sliding, aes(x = min_puntaje, y = 3, fill = mean_cor)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Number of Clusters (k)",
    fill = "Cor (lower is better)",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  )
dev.off()




##
##
##

# No paralelo ----------------------------------------------------------------

# --- Carga de librerías necesarias ---
library(dplyr)
library(ggplot2)
library(bootnet)
library(cluster) # Para el cálculo de silhouette
library(viridis) # Para escalas de color
library(qgraph)
library(igraph)

# --- Carga de datos ---
base_igi_bin <- read.csv("psy_net_files/base_igi_bin.csv")
descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR24","PAR25","PAR26","PAR27",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

# --- Función para procesar una ventana ---
process_window <- function(base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values) {
  
  resultados <- data.frame()
  
  for (start_val in 0:(max(base_igi_bin$puntaje_total) - window_size)) {
  #for (start_val in 10:14) {
    tryCatch({
      # Define la ventana actual
      end_val <- start_val + window_size - 1
      rango_puntaje <- start_val:end_val
      
      # Filtra datos para la ventana actual
      sub_bloque <- base_igi_bin %>%
        filter(puntaje_total %in% rango_puntaje)
      
      # Salta si no hay suficientes observaciones
      if (nrow(sub_bloque) < min_cluster_size) {
        resultados <- rbind(resultados, data.frame(start_val = start_val, error_message = "Insufficient observations"))
        next
      }
      
      # Bucle sobre los valores candidatos de k
      for (k in k_values) {
        message(paste("Processing start_val:", start_val, "k_value:", k))
        
        # Realiza clustering con k clusters
        matriz_binaria <- as.matrix(sub_bloque[, descriptivo_grupal])
        distancias <- stats::dist(matriz_binaria, method = "binary")
        clust_hier <- stats::hclust(distancias, method = "ward.D2")
        sub_bloque$sub_cluster <- cutree(clust_hier, k = k)
        
        # Salta si algún cluster tiene menos que el mínimo requerido de observaciones
        cluster_sizes <- table(sub_bloque$sub_cluster)
        if (any(cluster_sizes < min_cluster_size)) {
          resultados <- rbind(resultados, data.frame(start_val = start_val, k = k, error_message = "Cluster too small"))
          next
        }
        
        # Evalúa la calidad del cluster (e.g., puntuación silhouette)
        silhouette_score <- mean(cluster::silhouette(
          sub_bloque$sub_cluster, distancias)[, "sil_width"])
        
        # Estima redes de Ising para cada cluster y las compara
        network_models <- list()
        for (i in unique(sub_bloque$sub_cluster)) {
          cluster_data <- sub_bloque[sub_bloque$sub_cluster == i, descriptivo_grupal]
          network_models[[i]] <- bootnet::estimateNetwork(cluster_data, default = "IsingFit", tuning = 0.25)
        }
        
        # Compara redes usando correlación de Pearson (métrica simplificada)
        pairwise_comparisons <- combn(unique(sub_bloque$sub_cluster), 2, simplify = FALSE)
        pears_cor <- sapply(pairwise_comparisons, function(pair) {
          model_1 <- network_models[[pair[1]]]
          model_2 <- network_models[[pair[2]]]
          
          cor(as.vector(model_1$graph), as.vector(model_2$graph), use = "complete.obs")
        })
        
        # Calcula puntuación compuesta (menor es mejor)
        composite_score <- (1 - silhouette_score) + mean(pears_cor)
        
        resultados <- rbind(resultados, data.frame(
          start_val = start_val,
          end_val = end_val,
          k = k,
          silhouette_score = silhouette_score,
          composite_score = composite_score,
          cluster_sizes = paste(cluster_sizes, collapse = ", ")
        ))
      }
    }, error = function(e) {
      message(paste("Error in window:", start_val))
      resultados <- rbind(resultados, data.frame(start_val = start_val, error_message = e$message))
    })
  }
  
  resultados
}

# Define parámetros
window_size <- 10 # Tamaño de la ventana deslizante
min_cluster_size <- 50 # Tamaño mínimo del cluster para una estimación de red estable
k_values <- 2:3 # Rango de números de cluster a explorar

# --- Ejecución del análisis secuencial ---
time_init <- Sys.time()

results_grid <- process_window(base_igi_bin, descriptivo_grupal, window_size, min_cluster_size, k_values)

time_fin <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")  # 3.91 mins

write.csv(results_grid, "psy_net_files/optimal_k_grid_1.csv", row.names = FALSE)
results_grid <- read.csv("psy_net_files/optimal_k_grid_1.csv")

# Visualize results with a heatmap
ggplot(results_grid, aes(x = start_val, y = k, fill = composite_score)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x = "Window Start",
    y = "Number of Clusters (k)",
    fill = "Composite Score",
    title = "Optimal Parameters Heatmap {3,2,1},{0} -> {0},{1}"
  )

