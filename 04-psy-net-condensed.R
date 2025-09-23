# =============================================================================
# Proyecto Cárceles - Ising
# Script condensado: selección de k por ventana y análisis de sub-redes
# Entrada principal: psy_net_files/base_igi_bin_1.csv (esquema 1: {3,2,1}->{0}; {0}->{1})
#
# Qué hace:
# 1) Recorta extremos de puntaje IGI a P10–P90 para mayor estabilidad.
# 2) Barrido ventana × k (k=2..6): calcula silhouette por dos distancias (binary-matching y
#    Jaccard), media de similitud entre redes (correlación de Pearson entre matrices de aristas),
#    y composite = (1 - silhouette) + mean_cor.
#    - Entrega k óptimo por ventana y un k global sugerido, por distancia.
#    - Genera heatmaps (PDF) separados por distancia, con cruces rojas donde mean_cor > umbral.
# 3) Ciclo analítico para TODOS los k del rango (por ventana):
#    - Estima redes Ising por subcluster y calcula 4 indicadores de similitud:
#      pearson_signed, pearson_abs, sign_agreement, jaccard_pos/neg.
#    - Métricas de nodo: fuerza (abs y signed), deltas por par y correlaciones de ranking (Spearman).
#    - Chi-cuadrado + V de Cramér por variables externas (FDR-BH global).
# 4) Heatmap de activación por puntaje (P10–P90) con etiquetas de ítems coloreadas por comunidad
#    (Louvain en red global de referencia).
#
# Salidas (todas como .csv y .rds):
# - psy_net_files/optimal_k_grid_by_window_1.csv|rds        (ventana × k × distancia: silhouette, mean_cor, composite)
# - psy_net_files/optimal_k_by_window_1.csv|rds             (k óptimo por ventana + k global sugerido, por distancia)
# - psy_net_files/network_similarity_by_window.csv|rds      (4 métricas por ventana×k×par)
# - psy_net_files/node_deltas_by_window.csv|rds             (deltas de fuerza y correlaciones de ranking)
# - psy_net_files/chi_and_cramersV_by_window.csv|rds        (chi2, df, p, fisher_p, V, p_adj_BH)
#
# Gráficos (PDF, alta calidad) en psy_net_plots/:
# - optimal_k_mean_cor_1_binary.pdf  y optimal_k_mean_cor_1_jaccard.pdf
# - optimal_k_silhouette_score_1_binary.pdf y _jaccard.pdf
# - optimal_k_composite_score_1_binary.pdf y _jaccard.pdf
# - optimal_k_by_window_1_binary.pdf y _jaccard.pdf (línea de k óptimo vs start_val)
# - heatmap_prop_trimmed_by_comm.pdf (activación por puntaje con etiquetas por comunidad)
#
# Paralelización:
# - Por ventanas, con makeCluster(8, type="FORK"); evita paralelo anidado.
# - Se cronometra cada tramo con time_init/time_fin/difftime.
# =============================================================================

# --------- Parámetros --------------------------------------------------------
set.seed(123)

# Intenta evitar mensajes de macOS relacionados a MallocStackLogging
if (nzchar(Sys.getenv("MallocStackLogging"))) Sys.unsetenv("MallocStackLogging")
if (nzchar(Sys.getenv("MallocStackLoggingNoCompact"))) Sys.unsetenv("MallocStackLoggingNoCompact")

input_csv           <- "psy_net_files/base_igi_bin_1.csv"
out_dir_files       <- "psy_net_files"
out_dir_plots       <- "psy_net_plots"

# Ventaneo y clustering
window_size         <- 5 # 10
k_values            <- 2:6           # Rango de k a explorar y luego analizar
min_cluster_size    <- 50 # 75           # Mínimo por cluster y por ventana

# Trimming P10–P90
use_trim            <- TRUE
trim_q              <- c(0.10, 0.90)

# Distancia para el ciclo analítico completo (clustering de subredes); el barrido k calcula ambas
analysis_distance_method <- "binary"  # "binary" (matching normalizada) o "jaccard"

# Paralelización
n_cores             <- 8             # M4 Pro: 8 P-cores
use_parallel        <- TRUE          # TRUE para paralelizar por ventanas
cluster_type        <- "PSOCK"       # "PSOCK" recomendado en macOS/RStudio; "FORK" si ejecutas via Rscript

# Umbral para marcar mean_cor “alto” en heatmap (cruz roja)
cor_threshold       <- 0.7

# Variables externas para perfilamiento
external_variables  <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# Vector de ítems (43)
descriptivo_grupal <- c(
  "HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
  "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
  "FAM18","FAM19","FAM20","FAM21",
  "UTL22","UTL23",
  "PAR24","PAR25","PAR26","PAR27",
  "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
  "PRO36","PRO37","PRO38","PRO39",
  "PAT40","PAT41","PAT42","PAT43"
)

# --------- Librerías ---------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(viridis)
library(ggtext)      # para etiquetas coloreadas en eje
library(bootnet)     # estimateNetwork (IsingFit)
library(qgraph)
library(igraph)
library(cluster)     # silhouette
library(foreach)
library(doParallel)

# (sin dependencias externas para distancias; Jaccard y Binary-matching se implementan abajo)

# --------- Preparación de rutas ----------------------------------------------
# dir.create(out_dir_files, showWarnings = FALSE, recursive = TRUE)
# dir.create(out_dir_plots, showWarnings = FALSE, recursive = TRUE)

# --------- Utilidades --------------------------------------------------------
safe_write <- function(df, basepath) {
  # Guarda CSV y RDS
  readr::write_csv(df, paste0(basepath, ".csv"))
  saveRDS(df, paste0(basepath, ".rds"))
}

# Fuerzas de nodo a partir de matriz de pesos
node_strength <- function(W) {
  stopifnot(nrow(W) == ncol(W))
  s_abs <- rowSums(abs(W), na.rm = TRUE)
  s_sig <- rowSums(W, na.rm = TRUE)
  tibble(node = rownames(W) %||% as.character(seq_len(nrow(W))),
         strength_abs = s_abs,
         strength_signed = s_sig)
}

# 4 métricas de similitud entre dos grafos (mismas dimensiones y orden)
edge_metrics <- function(W1, W2, exclude_both_zero = FALSE) {
  stopifnot(all(dim(W1) == dim(W2)))
  ut <- upper.tri(W1, diag = FALSE)
  w1 <- as.numeric(W1[ut]); w2 <- as.numeric(W2[ut])
  if (exclude_both_zero) {
    keep <- !(w1 == 0 & w2 == 0)
    w1 <- w1[keep]; w2 <- w2[keep]
  }
  s1 <- sign(w1); s2 <- sign(w2)
  pearson_signed <- suppressWarnings(cor(w1, w2))
  pearson_abs    <- suppressWarnings(cor(abs(w1), abs(w2)))
  sign_agreement <- mean(s1 == s2)
  
  jaccard <- function(a, b) {
    inter <- sum(a & b)
    union <- sum(a | b)
    if (union == 0) return(NA_real_)
    inter/union
  }
  pos1 <- w1 > 0; pos2 <- w2 > 0
  neg1 <- w1 < 0; neg2 <- w2 < 0
  jaccard_pos <- jaccard(pos1, pos2)
  jaccard_neg <- jaccard(neg1, neg2)
  
  tibble(pearson_signed = pearson_signed,
         pearson_abs    = pearson_abs,
         sign_agreement = sign_agreement,
         jaccard_pos    = jaccard_pos,
         jaccard_neg    = jaccard_neg)
}

# Distancias para clustering (vectorizadas y rápidas)
# - Binary (matching normalizada): distancia = Hamming(p)/p
compute_distance_binary <- function(mat_bin) {
  M <- as.matrix(mat_bin)
  M[is.na(M)] <- 0
  M <- (M > 0) + 0L
  p <- ncol(M)
  rs <- as.numeric(rowSums(M))
  inter <- M %*% t(M)                         # |A ∩ B|
  hamming <- outer(rs, rs, "+") - 2*inter    # b + c
  D <- as.matrix(hamming / max(p, 1L))        # normaliza por #ítems
  diag(D) <- 0
  D
}

# - Jaccard (asymétrica): 1 - |A ∩ B| / |A ∪ B|
compute_distance_jaccard <- function(mat_bin) {
  M <- as.matrix(mat_bin)
  M[is.na(M)] <- 0
  M <- (M > 0) + 0L
  rs <- as.numeric(rowSums(M))
  inter <- M %*% t(M)                         # |A ∩ B|
  union <- outer(rs, rs, "+") - inter       # |A ∪ B|
  sim <- inter / pmax(union, 1)               # si union==0 -> 0
  sim[union == 0] <- 1                        # filas sin 1s: similitud 1 (dist 0)
  D <- 1 - as.matrix(sim)
  diag(D) <- 0
  D
}

# Selección de una distancia para el ciclo completo
compute_distance_single <- function(mat_bin, method = c("binary", "jaccard")) {
  method <- match.arg(method)
  if (method == "binary") compute_distance_binary(mat_bin) else compute_distance_jaccard(mat_bin)
}

# Silhouette medio
mean_silhouette <- function(cluster_labels, dist_matrix) {
  as.numeric(mean(cluster::silhouette(cluster_labels, as.dist(dist_matrix))[, "sil_width"]))
}

# Cramér's V
cramers_v <- function(tabla) {
  chi <- suppressWarnings(chisq.test(tabla, correct = FALSE))
  n <- sum(tabla)
  r <- nrow(tabla); c <- ncol(tabla)
  if (min(r, c) < 2) return(NA_real_)
  as.numeric(sqrt(chi$statistic / (n * min(r - 1, c - 1))))
}

# Fisher seguro: exacto solo 2x2; para tablas más grandes usa Monte Carlo para evitar FEXACT error 6
safe_fisher_p <- function(tabla, B = 10000L) {
  out <- list(p = NA_real_, method = NA_character_)
  # Intento exacto para 2x2
  if (nrow(tabla) == 2L && ncol(tabla) == 2L) {
    res <- tryCatch(stats::fisher.test(tabla), error = function(e) NULL)
    if (!is.null(res)) {
      out$p <- as.numeric(res$p.value)
      out$method <- "exact_2x2"
      return(out)
    }
  }
  # Monte Carlo para r×c grandes para evitar problemas de memoria/espacio
  res_mc <- tryCatch(stats::fisher.test(tabla, simulate.p.value = TRUE, B = B), error = function(e) NULL)
  if (!is.null(res_mc)) {
    out$p <- as.numeric(res_mc$p.value)
    out$method <- sprintf("monte_carlo_B=%d", B)
  }
  out
}

# --------- Carga de datos ----------------------------------------------------
base <- read.csv(input_csv, stringsAsFactors = FALSE)
if (!all(c(descriptivo_grupal, "puntaje_total") %in% names(base))) {
  stop("Faltan columnas esperadas en base_igi_bin_1.csv (ítems y/o puntaje_total).")
}

# --------- Trimming P10–P90 --------------------------------------------------
if (use_trim) {
  q <- quantile(base$puntaje_total, trim_q, na.rm = TRUE, type = 1)
  lo <- as.integer(q[1]); hi <- as.integer(q[2])
  base_trim <- base %>% filter(puntaje_total >= lo, puntaje_total <= hi)
} else {
  lo <- min(base$puntaje_total, na.rm = TRUE)
  hi <- max(base$puntaje_total, na.rm = TRUE)
  base_trim <- base
}
message(sprintf("Rango de puntaje utilizado: [%d, %d] (use_trim=%s).", lo, hi, use_trim))

# --------- Red de referencia y comunidades -----------------------------------
# Red global en base_trim (solo ítems), para colorear etiquetas por comunidad
ref_data <- base_trim[, descriptivo_grupal, drop = FALSE]
ref_net  <- bootnet::estimateNetwork(ref_data, default = "IsingFit", tuning = 0.25)
G        <- graph_from_adjacency_matrix(ref_net$graph, mode = "undirected", weighted = TRUE)
E(G)$weight <- abs(E(G)$weight)
comm     <- cluster_louvain(G, weights = E(G)$weight)
memb     <- membership(comm)  # named vector: item -> id_comunidad

# --------- Heatmap de activación (P10–P90) con etiquetas por comunidad -------
# Proporción de activación por puntaje IGI (dentro del rango lo..hi)
valores <- sort(unique(base_trim$puntaje_total))
# Matriz [items x puntajes]
prop_mat <- matrix(NA_real_, nrow = length(descriptivo_grupal), ncol = length(valores),
                   dimnames = list(descriptivo_grupal, as.character(valores)))
for (j in seq_along(valores)) {
  casos <- base_trim$puntaje_total == valores[j]
  denom <- sum(casos)
  if (denom > 0) {
    for (i in seq_along(descriptivo_grupal)) {
      z <- sum(base_trim[casos, descriptivo_grupal[i]], na.rm = TRUE) / denom
      prop_mat[i, j] <- z
    }
  }
}
prop_df <- as.data.frame(as.table(prop_mat))
colnames(prop_df) <- c("Variable", "PuntajeTotal", "Proporcion")
prop_df$PuntajeTotal <- as.integer(as.character(prop_df$PuntajeTotal))

# Informatividad por ítem
var_info <- prop_df %>%
  group_by(Variable) %>%
  summarise(sd_prop = sd(Proporcion, na.rm = TRUE),
            range_prop = max(Proporcion, na.rm = TRUE) - min(Proporcion, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(comunidad = as.integer(factor(memb[Variable], exclude = NULL)))

# Orden: (comunidad, -range_prop)
ord_vars <- var_info %>% arrange(comunidad, desc(range_prop)) %>% pull(Variable)
prop_df  <- prop_df %>%
  left_join(var_info %>% select(Variable, comunidad), by = "Variable") %>%
  mutate(Variable = factor(Variable, levels = ord_vars))

# Paleta por comunidad y etiquetas coloreadas (ggtext)
comms <- sort(unique(na.omit(var_info$comunidad)))
pal   <- scales::hue_pal()(length(comms))
col_map <- setNames(pal, comms)
label_map <- var_info %>%
  mutate(label = ifelse(is.na(comunidad),
                        Variable,
                        paste0("<span style='color:", col_map[as.character(comunidad)], ";'>", Variable, "</span>"))) %>%
  select(Variable, label)
prop_df <- prop_df %>% left_join(label_map, by = "Variable")

pdf(file.path(out_dir_plots, "heatmap_prop_trimmed_by_comm.pdf"), width = 10, height = 8)
print(
  ggplot(prop_df, aes(x = PuntajeTotal, y = label, fill = Proporcion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = sprintf("Proporción de activación de ítems (IGI %d–%d), etiquetas por comunidad", lo, hi),
      x = "Puntaje IGI", y = "Ítem", fill = "Proporción") +
    theme_minimal() +
    theme(axis.text.y = ggtext::element_markdown(),
    axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

# --------- Ventanas a evaluar -------------------------------------------------
start_val_vec <- seq(from = lo, to = hi - window_size + 1, by = 1)

# --------- Función: procesa una ventana (barrido k) --------------------------
process_window_kgrid <- function(start_val) {
  end_val <- start_val + window_size - 1
  rango   <- start_val:end_val
  sub     <- base_trim %>% filter(puntaje_total %in% rango)
  
  if (nrow(sub) < min_cluster_size) {
    grid_empty <- tidyr::crossing(k = k_values, distance = c("binary", "jaccard"))
    return(grid_empty %>%
             mutate(start_val = start_val, end_val = end_val,
                    silhouette_score = NA_real_, mean_pears_cor = NA_real_,
                    composite_score = NA_real_, cluster_sizes = NA_character_,
                    error_message = "Insufficient observations"))
  }
  
  mat_bin <- as.matrix(sub[, descriptivo_grupal, drop = FALSE])
  D_bin   <- compute_distance_binary(mat_bin)
  D_jac   <- compute_distance_jaccard(mat_bin)
  
  do_one <- function(k, D_here, dist_name) {
    cl <- hclust(as.dist(D_here), method = "ward.D2")
    sub$sub_cluster <- cutree(cl, k = k)
    sizes <- table(sub$sub_cluster)
    if (any(sizes < min_cluster_size)) {
      return(tibble(start_val = start_val, end_val = end_val, k = k, distance = dist_name,
                    silhouette_score = NA_real_, mean_pears_cor = NA_real_,
                    composite_score = NA_real_, cluster_sizes = paste(sizes, collapse = ", "),
                    error_message = "Cluster too small"))
    }
    sil <- mean_silhouette(sub$sub_cluster, D_here)
    # Redes por subcluster
    models <- map(sort(unique(sub$sub_cluster)), function(cc) {
      dat_cc <- sub %>% filter(sub_cluster == cc) %>% select(all_of(descriptivo_grupal))
      bootnet::estimateNetwork(dat_cc, default = "IsingFit", tuning = 0.25)
    })
    # Correlaciones par-a-par entre grafos (vectorizados)
    pairs <- combn(seq_along(models), 2, simplify = FALSE)
    corrs <- map_dbl(pairs, function(p) {
      cor(as.vector(models[[p[1]]]$graph), as.vector(models[[p[2]]]$graph), use = "complete.obs")
    })
    mean_cor <- mean(corrs)
    tibble(start_val = start_val, end_val = end_val, k = k, distance = dist_name,
           silhouette_score = sil, mean_pears_cor = mean_cor,
           composite_score = (1 - sil) + mean_cor,
           cluster_sizes = paste(sizes, collapse = ", "),
           error_message = NA_character_)
  }
  out_bin <- map_dfr(k_values, ~do_one(.x, D_bin, "binary"))
  out_jac <- map_dfr(k_values, ~do_one(.x, D_jac, "jaccard"))
  out <- bind_rows(out_bin, out_jac)
  out
}

# --------- BARRIDO VENTANA × k (paralelizado por ventanas) -------------------
if (use_parallel) {
  cl <- makeCluster(n_cores, type = cluster_type)
  registerDoParallel(cl)
}
time_init <- Sys.time()
if (use_parallel) {
  kgrid <- foreach(v = start_val_vec, .combine = "rbind",
                   .packages = c("dplyr", "tibble", "tidyr", "purrr", "bootnet", "cluster"),
                   .export = c(
                     "process_window_kgrid", "compute_distance_binary", "compute_distance_jaccard",
                     "mean_silhouette", "edge_metrics", "node_strength", "cramers_v",
                     "descriptivo_grupal", "k_values", "min_cluster_size", "window_size",
                     "base_trim", "cor_threshold"
                   )) %dopar% {
                     process_window_kgrid(v)
                   }
} else {
  kgrid <- map_dfr(start_val_vec, process_window_kgrid)
}
time_fin <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")
if (use_parallel) stopCluster(cl)
message(sprintf("Tiempo barrido ventana×k: %s", as.character(time_total)))

# Guardar grid por ventana×k
kgrid <- as_tibble(kgrid)
safe_write(kgrid, file.path(out_dir_files, "optimal_k_grid_by_window_1"))

# k óptimo por ventana (mínimo composite); k global sugerido
opt_by_win <- kgrid %>%
  group_by(distance, start_val, end_val) %>%
  filter(!is.na(composite_score)) %>%
  slice_min(order_by = composite_score, n = 1, with_ties = TRUE) %>%
  mutate(k_opt = k) %>%
  select(distance, start_val, end_val, k_opt, silhouette_score, mean_pears_cor, composite_score, cluster_sizes)

# k global por dos criterios
k_global_by_mode <- opt_by_win %>%
  group_by(distance) %>%
  count(k_opt, name = "freq") %>%
  group_by(distance) %>%
  slice_max(order_by = freq, n = 1, with_ties = TRUE) %>%
  summarise(distance = first(distance), k_global = paste(sort(k_opt), collapse = "/"), .groups = "drop") %>%
  mutate(method = "mode_of_window_optima")

k_global_by_mean <- kgrid %>%
  group_by(distance, k) %>%
  summarise(mean_composite = mean(composite_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(distance) %>%
  slice_min(order_by = mean_composite, n = 1, with_ties = TRUE) %>%
  summarise(distance = first(distance), k_global = paste(sort(k), collapse = "/"), .groups = "drop") %>%
  mutate(method = "min_mean_composite")

opt_summary <- bind_rows(k_global_by_mode, k_global_by_mean) %>%
  select(method, distance, k_global)

safe_write(opt_by_win, file.path(out_dir_files, "optimal_k_by_window_1"))
safe_write(opt_summary, file.path(out_dir_files, "optimal_k_global_summary_1"))

# --------- Heatmaps (PDF) con overlay de cruces rojas (> cor_threshold) ------
for (dist_name in c("binary", "jaccard")) {
  kd <- kgrid %>% filter(distance == dist_name)
  overlay_pts <- kd %>% filter(!is.na(mean_pears_cor), mean_pears_cor > cor_threshold)
  pdf(file.path(out_dir_plots, sprintf("optimal_k_mean_cor_1_%s.pdf", dist_name)), width = 10, height = 7)
  print(
    ggplot(kd, aes(x = start_val, y = k, fill = mean_pears_cor)) +
      geom_tile() +
      geom_point(data = overlay_pts, aes(x = start_val, y = k),
                 inherit.aes = FALSE, shape = 4, color = "red", stroke = 0.9, size = 1.8) +
                 scale_fill_viridis_c(option = "plasma", direction = -1) +
      labs(x = "Window Start (IGI)", y = "k",
           fill = "Mean Corr (lower=better)",
           title = sprintf("Mean Pearson Correlation (lower=better) - %s", dist_name)) +
      theme_minimal()
  )
  dev.off()
  
  pdf(file.path(out_dir_plots, sprintf("optimal_k_silhouette_score_1_%s.pdf", dist_name)), width = 10, height = 7)
  print(
    ggplot(kd, aes(x = start_val, y = k, fill = silhouette_score)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", direction = -1) +
      labs(x = "Window Start (IGI)", y = "k",
           fill = "Silhouette (higher=better)",
           title = sprintf("Silhouette by Window × k - %s", dist_name)) +
      theme_minimal()
  )
  dev.off()
  
  pdf(file.path(out_dir_plots, sprintf("optimal_k_composite_score_1_%s.pdf", dist_name)), width = 10, height = 7)
  print(
    ggplot(kd, aes(x = start_val, y = k, fill = composite_score)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", direction = -1) +
      labs(x = "Window Start (IGI)", y = "k",
           fill = "Composite (lower=better)",
           title = sprintf("Composite Score (1 - silhouette) + mean_cor - %s", dist_name)) +
      theme_minimal()
  )
  dev.off()
  
  # k óptimo por ventana (línea)
  od <- opt_by_win %>% filter(distance == dist_name)
  pdf(file.path(out_dir_plots, sprintf("optimal_k_by_window_1_%s.pdf", dist_name)), width = 10, height = 4)
  print(
    ggplot(od, aes(x = start_val, y = k_opt)) +
      geom_line(color = "steelblue") +
      geom_point(color = "steelblue") +
      scale_y_continuous(breaks = sort(unique(k_values))) +
      labs(x = "Window Start (IGI)", y = "k óptimo (por ventana)",
           title = sprintf("k óptimo por ventana (mínimo composite) - %s", dist_name)) +
      theme_minimal()
  )
  dev.off()
}

# --------- CICLO ANALÍTICO: todos los k, todas las ventanas ------------------
process_window_full <- function(start_val) {
  end_val <- start_val + window_size - 1
  rango   <- start_val:end_val
  sub     <- base_trim %>% filter(puntaje_total %in% rango)
  
  if (nrow(sub) < min_cluster_size) {
    return(list(
      metrics = tibble(), nodes = tibble(), chi = tibble()
    ))
  }
  
  mat_bin <- as.matrix(sub[, descriptivo_grupal, drop = FALSE])
  D       <- compute_distance_single(mat_bin, method = analysis_distance_method)
  
  metrics_rows <- list()
  node_rows    <- list()
  chi_rows     <- list()
  
  for (k in k_values) {
    cl <- hclust(as.dist(D), method = "ward.D2")
    sub$sub_cluster <- cutree(cl, k = k)
    sizes <- table(sub$sub_cluster)
    if (any(sizes < min_cluster_size)) next
    
    # Redes por subcluster
    clusters <- sort(unique(sub$sub_cluster))
    models <- map(clusters, function(cc) {
      dat_cc <- sub %>% filter(sub_cluster == cc) %>% select(all_of(descriptivo_grupal))
      bootnet::estimateNetwork(dat_cc, default = "IsingFit", tuning = 0.25)
    })
    names(models) <- paste0("c", clusters)
    
    # Métricas de aristas (4 indicadores) por par
    pairs <- combn(names(models), 2, simplify = FALSE)
    for (p in pairs) {
      W1 <- models[[p[1]]]$graph
      W2 <- models[[p[2]]]$graph
      m  <- edge_metrics(W1, W2)
      metrics_rows[[length(metrics_rows) + 1]] <- m %>%
        mutate(start_val = start_val, end_val = end_val, k = k, pair = paste(p, collapse = "-"))
      # Métricas de nodo y deltas
      ns1 <- node_strength(W1); colnames(ns1)[-1] <- paste0(colnames(ns1)[-1], "_1")
      ns2 <- node_strength(W2); colnames(ns2)[-1] <- paste0(colnames(ns2)[-1], "_2")
      nd  <- ns1 %>% inner_join(ns2, by = "node") %>%
        mutate(delta_strength_abs    = strength_abs_1    - strength_abs_2,
               delta_strength_signed = strength_signed_1 - strength_signed_2)
      rank_corr_abs    <- suppressWarnings(cor(rank(ns1$strength_abs), rank(ns2$strength_abs), method = "spearman"))
      rank_corr_signed <- suppressWarnings(cor(rank(ns1$strength_signed), rank(ns2$strength_signed), method = "spearman"))
      node_rows[[length(node_rows) + 1]] <- nd %>%
        mutate(start_val = start_val, end_val = end_val, k = k, pair = paste(p, collapse = "-"),
               rank_corr_abs = rank_corr_abs, rank_corr_signed = rank_corr_signed)
    }
    
    # Chi-cuadrado + Cramér V por variables externas (contra sub_cluster)
    for (var_externa in external_variables) {
      if (!var_externa %in% names(sub)) next
      tab <- table(sub[[var_externa]], sub$sub_cluster)
      if (all(dim(tab) >= 2L) && sum(tab) > 0) {
        chi  <- suppressWarnings(chisq.test(tab))
        X2   <- as.numeric(chi$statistic); df <- as.integer(chi$parameter)
        pval <- as.numeric(chi$p.value)
        # Fisher como test alternativo cuando hay esperados pequeños
        fisher_p <- NA_real_
        fisher_method <- NA_character_
        if (any(chi$expected < 5)) {
          sf <- safe_fisher_p(tab)
          fisher_p <- sf$p
          fisher_method <- sf$method
        }
        V <- cramers_v(tab)
        chi_rows[[length(chi_rows) + 1]] <- tibble(
          start_val = start_val, end_val = end_val, k = k,
          variable = var_externa, X2 = X2, df = df, p_value = pval,
          fisher_p = fisher_p, fisher_method = fisher_method, V = V, n = sum(tab)
        )
      }
    }
  }
  
  list(
    metrics = bind_rows(metrics_rows),
    nodes   = bind_rows(node_rows),
    chi     = bind_rows(chi_rows)
  )
}

# Versión segura: captura errores por ventana y devuelve tibbles vacíos
process_window_full_safe <- function(start_val) {
  tryCatch(
    process_window_full(start_val),
    error = function(e) {
      warning(sprintf("process_window_full falló en ventana %d: %s", start_val, conditionMessage(e)))
      list(metrics = tibble(), nodes = tibble(), chi = tibble())
    }
  )
}

if (use_parallel) {
  cl <- makeCluster(n_cores, type = cluster_type)
  registerDoParallel(cl)
}
time_init <- Sys.time()
if (use_parallel) {
  full_results <- foreach(v = start_val_vec, .combine = "c",
                          .packages = c("dplyr", "tibble", "tidyr", "purrr", "bootnet", "cluster", "igraph"),
                          .export = c(
                            "process_window_full", "process_window_full_safe", "compute_distance_single", "analysis_distance_method",
                            "descriptivo_grupal", "k_values", "min_cluster_size", "window_size",
                            "external_variables", "cramers_v", "edge_metrics", "node_strength",
                            "base_trim"
                          )) %dopar% {
                            list(process_window_full_safe(v))
                          }
} else {
  full_results <- map(start_val_vec, process_window_full_safe)
}
time_fin <- Sys.time()
time_total <- difftime(time_fin, time_init, units = "auto")
if (use_parallel) stopCluster(cl)
message(sprintf("Tiempo ciclo analítico completo: %s", as.character(time_total)))

# Unir resultados
metrics_results <- bind_rows(lapply(full_results, `[[`, "metrics"))
node_deltas     <- bind_rows(lapply(full_results, `[[`, "nodes"))
chi_results     <- bind_rows(lapply(full_results, `[[`, "chi"))

# Ajuste FDR-BH global por variable externa
if (nrow(chi_results) > 0) {
  chi_results <- chi_results %>%
    group_by(variable) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
           fisher_p_adj_BH = if (all(is.na(fisher_p))) NA_real_ else p.adjust(fisher_p, method = "BH")) %>%
    ungroup()
}

# Guardar salidas analíticas
safe_write(metrics_results, file.path(out_dir_files, "network_similarity_by_window"))
safe_write(node_deltas,     file.path(out_dir_files, "node_deltas_by_window"))
safe_write(chi_results,     file.path(out_dir_files, "chi_and_cramersV_by_window"))

# --------- Resumen final en consola ------------------------------------------
message("Resumen k global sugerido:")
print(opt_summary)
