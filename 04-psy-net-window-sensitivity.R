# =============================================================================
# 05-psy-net-window-sensitivity.R
# Sensibilidad al tamaño de ventana: window = 3..10
# - Excluye ítems de baja varianza: HD6 y PAR26
# - Para cada window:
#   * Barrido ventana × k (k=2..6) con distancias binary y jaccard
#   * Guarda CSV/RDS y heatmaps por distancia en carpetas específicas
#   * Selecciona una ventana representativa por k (mínimo composite con distancia de análisis)
#   * Genera:
#     - Un heatmap de activación (P10–P90) con etiquetas de ítems coloreadas por asignación a la subred (k-colores)
#     - Un PDF por k con k componentes (una subred por panel), todos con el mismo layout
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(bootnet)
  library(qgraph)
  library(igraph)
  library(cluster)
  library(foreach)
  library(doParallel)
  library(grid)
  library(gridExtra)
})

set.seed(123)

# Entradas / Salidas base
input_csv      <- "psy_net_files/base_igi_bin_1.csv"
root_files_dir <- "psy_net_files"
root_plots_dir <- "psy_net_plots"

# Parámetros
window_sizes          <- 3:10
k_values              <- 2:6
min_cluster_size      <- 50          # ajusta si hace falta
use_trim              <- TRUE
trim_q                <- c(0.10, 0.90)
analysis_distance_method <- "jaccard" # para seleccionar la ventana representativa por k
n_cores               <- 8
use_parallel          <- TRUE
cluster_type          <- "PSOCK"
cor_threshold         <- 0.70

# Ítems, excluyendo HD6 y PAR26
descriptivo_grupal_all <- c(
  "HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
  "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
  "FAM18","FAM19","FAM20","FAM21",
  "UTL22","UTL23",
  "PAR24","PAR25","PAR26","PAR27",
  "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
  "PRO36","PRO37","PRO38","PRO39",
  "PAT40","PAT41","PAT42","PAT43"
)
drop_items <- c("HD6","PAR24")
descriptivo_grupal <- setdiff(descriptivo_grupal_all, drop_items)

# Utilidades
safe_write <- function(df, basepath) {
  readr::write_csv(df, paste0(basepath, ".csv"))
  saveRDS(df, paste0(basepath, ".rds"))
}

# Distancias para clustering (vectorizadas rápidas)
compute_distance_binary <- function(mat_bin) {
  M <- as.matrix(mat_bin)
  M[is.na(M)] <- 0
  M <- (M > 0) + 0L
  p <- ncol(M)
  rs <- as.numeric(rowSums(M))
  inter <- M %*% t(M)
  hamming <- outer(rs, rs, "+") - 2*inter
  D <- as.matrix(hamming / max(p, 1L))
  diag(D) <- 0
  D
}
compute_distance_jaccard <- function(mat_bin) {
  M <- as.matrix(mat_bin)
  M[is.na(M)] <- 0
  M <- (M > 0) + 0L
  rs <- as.numeric(rowSums(M))
  inter <- M %*% t(M)
  union <- outer(rs, rs, "+") - inter
  sim <- inter / pmax(union, 1)
  sim[union == 0] <- 1
  D <- 1 - as.matrix(sim)
  diag(D) <- 0
  D
}
compute_distance_single <- function(mat_bin, method = c("binary","jaccard")) {
  method <- match.arg(method)
  if (method == "binary") compute_distance_binary(mat_bin) else compute_distance_jaccard(mat_bin)
}

mean_silhouette <- function(cluster_labels, dist_matrix) {
  as.numeric(mean(cluster::silhouette(cluster_labels, as.dist(dist_matrix))[, "sil_width"]))
}

node_strength <- function(W) {
  stopifnot(nrow(W) == ncol(W))
  rn <- rownames(W)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(W)))
  tibble(
    node = rn,
    strength_abs = rowSums(abs(W), na.rm = TRUE),
    strength_signed = rowSums(W, na.rm = TRUE)
  )
}

edge_metrics <- function(W1, W2) {
  ut <- upper.tri(W1, diag = FALSE)
  w1 <- as.numeric(W1[ut]); w2 <- as.numeric(W2[ut])
  tibble(
    pearson_signed = suppressWarnings(cor(w1, w2)),
    pearson_abs    = suppressWarnings(cor(abs(w1), abs(w2)))
  )
}

# Carga base y trimming
stopifnot(file.exists(input_csv))
base <- read.csv(input_csv, stringsAsFactors = FALSE)
stopifnot(all(c(descriptivo_grupal, "puntaje_total") %in% names(base)))

if (use_trim) {
  q <- quantile(base$puntaje_total, trim_q, na.rm = TRUE, type = 1)
  lo <- as.integer(q[1]); hi <- as.integer(q[2])
  base_trim <- base %>% filter(puntaje_total >= lo, puntaje_total <= hi)
} else {
  lo <- min(base$puntaje_total, na.rm = TRUE)
  hi <- max(base$puntaje_total, na.rm = TRUE)
  base_trim <- base
}
message(sprintf("Rango IGI usado: [%d, %d] (use_trim=%s).", lo, hi, use_trim))

# Heatmap de activación (P10–P90) auxiliar: proporción por puntaje
compute_activation_df <- function(df_items, puntajes, items) {
  valores <- sort(unique(puntajes))
  prop_mat <- matrix(NA_real_, nrow = length(items), ncol = length(valores),
                     dimnames = list(items, as.character(valores)))
  for (j in seq_along(valores)) {
    casos <- puntajes == valores[j]
    denom <- sum(casos)
    if (denom > 0) {
      for (i in seq_along(items)) {
        prop_mat[i, j] <- sum(df_items[casos, items[i]], na.rm = TRUE)/denom
      }
    }
  }
  prop_df <- as.data.frame(as.table(prop_mat))
  names(prop_df) <- c("Variable","PuntajeTotal","Proporcion")
  prop_df$PuntajeTotal <- as.integer(as.character(prop_df$PuntajeTotal))
  prop_df
}

# Para el heatmap coloreado por subred: asigna cada ítem al subcluster donde tiene mayor fuerza_abs
# Asigna cada ítem a la subred donde tiene mayor fuerza absoluta.
# Devuelve un data.frame con columnas:
# - "Variable": nombre del ítem
# - "asignacion": entero en 1..k indicando subred
assign_items_to_subnet <- function(models_list) {
  # models_list: lista de objetos bootnet (una por subcluster)
  strengths <- lapply(models_list, function(m) {
    ns <- node_strength(m$graph)
    ns <- ns %>% select(node, strength_abs)
    ns
  })
  df <- reduce(strengths, function(a,b) full_join(a,b, by = "node"))
  names(df) <- c("node", paste0("s", seq_len(length(models_list))))
  # Argmax por fila
  best <- apply(as.matrix(df[,-1]), 1, function(v) {
    if (all(is.na(v))) NA_integer_ else which.max(v)
  })
  tibble(Variable = df$node, asignacion = best)
}

# K-grid para una ventana start_val y window_size w
process_window_kgrid <- function(start_val, w) {
  end_val <- start_val + w - 1
  sub <- base_trim %>% filter(puntaje_total >= start_val, puntaje_total <= end_val)
  if (nrow(sub) < min_cluster_size) {
    return(tidyr::crossing(
      k = k_values,
      distance = c("binary","jaccard")
    ) %>% mutate(
      start_val = start_val, end_val = end_val,
      silhouette_score = NA_real_, mean_pears_cor = NA_real_,
      composite_score = NA_real_, cluster_sizes = NA_character_,
      error_message = "Insufficient observations"
    ))
  }
  mat_bin <- as.matrix(sub[, descriptivo_grupal, drop = FALSE])
  D_bin <- compute_distance_binary(mat_bin)
  D_jac <- compute_distance_jaccard(mat_bin)

  do_one <- function(k, D_here, dist_name) {
    cl <- hclust(as.dist(D_here), method = "ward.D2")
    sub$sub_cluster <- cutree(cl, k = k)
    sizes <- table(sub$sub_cluster)
    if (any(sizes < min_cluster_size)) {
      return(tibble(
        start_val = start_val, end_val = end_val, k = k, distance = dist_name,
        silhouette_score = NA_real_, mean_pears_cor = NA_real_,
        composite_score = NA_real_, cluster_sizes = paste(sizes, collapse=", "),
        error_message = "Cluster too small"
      ))
    }
    sil <- mean_silhouette(sub$sub_cluster, D_here)
    # Redes por subcluster
    clusters <- sort(unique(sub$sub_cluster))
    models <- lapply(clusters, function(cc) {
      dat_cc <- sub %>% filter(sub_cluster == cc) %>% select(all_of(descriptivo_grupal))
      bootnet::estimateNetwork(dat_cc, default = "IsingFit", tuning = 0.25)
    })
    pairs <- combn(seq_along(models), 2, simplify = FALSE)
    corrs <- vapply(pairs, function(p) {
      cor(as.vector(models[[p[1]]]$graph), as.vector(models[[p[2]]]$graph), use = "complete.obs")
    }, numeric(1))
    mean_cor <- mean(corrs)
    tibble(
      start_val = start_val, end_val = end_val, k = k, distance = dist_name,
      silhouette_score = sil, mean_pears_cor = mean_cor,
      composite_score = (1 - sil) + mean_cor,
      cluster_sizes = paste(sizes, collapse = ", "),
      error_message = NA_character_
    )
  }

  out_bin <- map_dfr(k_values, ~do_one(.x, D_bin, "binary"))
  out_jac <- map_dfr(k_values, ~do_one(.x, D_jac, "jaccard"))
  bind_rows(out_bin, out_jac)
}

# Ploteo: heatmaps por distancia
plot_heatmaps_distance <- function(kgrid_w, plots_dir) {
  for (dist_name in c("binary","jaccard")) {
    kd <- kgrid_w %>% filter(distance == dist_name)
    if (nrow(kd) == 0) next
    overlay_pts <- kd %>% filter(!is.na(mean_pears_cor), mean_pears_cor > cor_threshold)

    pdf(file.path(plots_dir, sprintf("optimal_k_mean_cor_%s.pdf", dist_name)), width = 10, height = 7)
    print(
      ggplot(kd, aes(x = start_val, y = k, fill = mean_pears_cor)) +
        geom_tile() +
        geom_point(data = overlay_pts, aes(x = start_val, y = k),
                   inherit.aes = FALSE, shape = 4, color = "red", stroke = 0.9, size = 1.8) +
        scale_fill_viridis_c(option = "plasma", direction = -1) +
        labs(x = "Window Start (IGI)", y = "k",
             fill = "Mean Corr (lower=better)",
             title = sprintf("Mean Pearson Correlation - %s", dist_name)) +
        theme_minimal()
    )
    dev.off()

    pdf(file.path(plots_dir, sprintf("optimal_k_silhouette_%s.pdf", dist_name)), width = 10, height = 7)
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

    pdf(file.path(plots_dir, sprintf("optimal_k_composite_%s.pdf", dist_name)), width = 10, height = 7)
    print(
      ggplot(kd, aes(x = start_val, y = k, fill = composite_score)) +
        geom_tile() +
        scale_fill_viridis_c(option = "plasma", direction = -1) +
        labs(x = "Window Start (IGI)", y = "k",
             fill = "Composite (lower=better)",
             title = sprintf("Composite (1 - silhouette) + mean_cor - %s", dist_name)) +
        theme_minimal()
    )
    dev.off()
  }
}

# Red con mismo layout para k subredes en una ventana concreta (base graphics)
plot_k_networks_onepage <- function(models, plot_path, title = NULL) {
  # Layout: usa matriz de pesos promedio (abs) para fijar posiciones
  k <- length(models)
  if (k == 0) return(invisible(NULL))
  Ws <- lapply(models, function(m) abs(m$graph))
  W_avg <- Reduce("+", Ws) / k
  # qgraph en modo 'DoNotPlot' para calcular layout común
  q0 <- qgraph::qgraph(W_avg, DoNotPlot = TRUE)
  coords <- q0$layout

  # Abrir dispositivo PDF y dibujar con base graphics en paneles
  pdf(plot_path, width = 3.5 * k, height = 4.2)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, k), mar = c(1, 1, 1, 1), oma = c(0, 0, 2, 0))

  for (i in seq_along(models)) {
    ok <- TRUE
    tryCatch({
      qgraph::qgraph(models[[i]]$graph,
                     layout = coords,
                     labels = colnames(models[[i]]$graph),
                     label.cex = 0.6,
                     edge.color = "grey50",
                     color = "skyblue",
                     vsize = 4,
                     esize = 3,
                     details = FALSE)
    }, error = function(e) {
      ok <<- FALSE
    })
    if (!ok) {
      plot.new(); box(); title(main = "(panel vacío)", cex.main = 0.8)
    }
  }
  if (!is.null(title)) graphics::mtext(title, outer = TRUE, line = 0.1, font = 2)
}

# Heatmap de activación con etiquetas coloreadas por asignación a subred (para k seleccionado)
plot_activation_colored_by_assignment <- function(base_trim, items, assign_df, plots_dir, w, k_sel, start_sel, end_sel) {
  prop_df <- compute_activation_df(base_trim[, items, drop = FALSE], base_trim$puntaje_total, items)
  # Colores por asignación
  k <- max(assign_df$asignacion, na.rm = TRUE)
  pal <- scales::hue_pal()(k)
  col_map <- setNames(pal, as.character(seq_len(k)))
  label_map <- assign_df %>%
    mutate(label = ifelse(is.na(asignacion),
                          Variable,
                          paste0("<span style='color:", col_map[as.character(asignacion)], ";'>", Variable, "</span>"))) %>%
    select(Variable, label)
  # Ordena por asignación y por rango de proporción
  var_info <- prop_df %>%
    group_by(Variable) %>%
    summarise(range_prop = max(Proporcion, na.rm = TRUE) - min(Proporcion, na.rm = TRUE), .groups = "drop") %>%
    left_join(assign_df, by = c("Variable")) %>%
    arrange(asignacion, desc(range_prop))
  prop_df2 <- prop_df %>%
    left_join(label_map, by = "Variable") %>%
    mutate(Variable = factor(Variable, levels = var_info$Variable))

  pdf(file.path(plots_dir, sprintf("heatmap_prop_by_assignment_w%d_k%d_start%d_end%d.pdf", w, k_sel, start_sel, end_sel)),
      width = 11, height = 8.5)
  print(
    ggplot(prop_df2, aes(x = PuntajeTotal, y = label, fill = Proporcion)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(title = sprintf("Proporción de activación (IGI %d–%d), etiquetas por subred (k=%d)", min(base_trim$puntaje_total), max(base_trim$puntaje_total), k_sel),
           x = "Puntaje IGI", y = "Ítem", fill = "Proporción") +
      theme_minimal() +
      theme(axis.text.y = ggtext::element_markdown(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

# -----------------------------------------------------------------------------#
# Loop principal por tamaños de ventana
# -----------------------------------------------------------------------------#

if (use_parallel) {
  cl <- makeCluster(n_cores, type = cluster_type)
  registerDoParallel(cl)
}

for (w in window_sizes) {
  message(sprintf("=== window = %d ===", w))
  files_dir <- file.path(root_files_dir, sprintf("files_window_%d", w))
  plots_dir <- file.path(root_plots_dir, sprintf("plots_window_%d", w))
  dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  start_vals <- seq(from = lo, to = hi - w + 1, by = 1)

  # k-grid (por distancia) para todas las ventanas de este tamaño
  time_init <- Sys.time()
  if (use_parallel) {
    kgrid_w <- foreach(v = start_vals, .combine = "rbind",
                       .packages = c("dplyr","tibble","tidyr","purrr","bootnet","cluster"),
                       .export = c("base_trim","descriptivo_grupal","k_values","min_cluster_size",
                                   "compute_distance_binary","compute_distance_jaccard","process_window_kgrid",
                                   "mean_silhouette")) %dopar% {
      process_window_kgrid(v, w)
    }
  } else {
    kgrid_w <- map_dfr(start_vals, ~process_window_kgrid(.x, w))
  }
  time_fin <- Sys.time()
  message(sprintf("k-grid (w=%d) listo en %s", w, as.character(difftime(time_fin, time_init, units = "auto"))))

  kgrid_w <- as_tibble(kgrid_w)
  safe_write(kgrid_w, file.path(files_dir, "optimal_k_grid_by_window"))

  # k óptimo por ventana y k global sugeridos (por distancia)
  opt_by_win <- kgrid_w %>%
    group_by(distance, start_val, end_val) %>%
    filter(!is.na(composite_score)) %>%
    slice_min(order_by = composite_score, n = 1, with_ties = TRUE) %>%
    mutate(k_opt = k) %>%
    select(distance, start_val, end_val, k_opt, silhouette_score, mean_pears_cor, composite_score, cluster_sizes) %>%
    ungroup()
  safe_write(opt_by_win, file.path(files_dir, "optimal_k_by_window"))

  opt_summary <- bind_rows(
    opt_by_win %>%
      count(distance, k_opt, name = "freq") %>%
      group_by(distance) %>%
      slice_max(order_by = freq, n = 1, with_ties = TRUE) %>%
      summarise(distance = first(distance), k_global = paste(sort(k_opt), collapse="/"), method = "mode_of_window_optima", .groups = "drop"),
    kgrid_w %>%
      group_by(distance, k) %>%
      summarise(mean_composite = mean(composite_score, na.rm = TRUE), .groups = "drop") %>%
      group_by(distance) %>%
      slice_min(order_by = mean_composite, n = 1, with_ties = TRUE) %>%
      summarise(distance = first(distance), k_global = paste(sort(k), collapse="/"), method = "min_mean_composite", .groups = "drop")
  ) %>% select(method, distance, k_global)
  safe_write(opt_summary, file.path(files_dir, "optimal_k_global_summary"))

  # Heatmaps por distancia
  plot_heatmaps_distance(kgrid_w, plots_dir)

  # Para los plots avanzados (heatmap coloreado por asignación y redes por k),
  # elegimos una ventana "representativa" por k con la distancia de análisis:
  # - Tomamos la ventana con composite mínimo para (distance=analysis_distance_method, k fijo).
  kd_sel <- kgrid_w %>% filter(distance == analysis_distance_method)

  for (k_sel in k_values) {
    cand <- kd_sel %>% filter(k == k_sel, !is.na(composite_score))
    if (nrow(cand) == 0) next
    row_best <- cand %>% slice_min(order_by = composite_score, n = 1, with_ties = TRUE) %>% slice(1)
    start_sel <- row_best$start_val[1]; end_sel <- row_best$end_val[1]

    # Estimar subredes en la ventana seleccionada
    sub <- base_trim %>% filter(puntaje_total >= start_sel, puntaje_total <= end_sel)
    if (nrow(sub) < min_cluster_size) next

    D <- compute_distance_single(as.matrix(sub[, descriptivo_grupal, drop = FALSE]), method = analysis_distance_method)
    cl <- hclust(as.dist(D), method = "ward.D2")
    sub$sub_cluster <- cutree(cl, k = k_sel)
    sizes <- table(sub$sub_cluster)
    if (any(sizes < min_cluster_size)) next

    clusters <- sort(unique(sub$sub_cluster))
    models <- lapply(clusters, function(cc) {
      dat_cc <- sub %>% filter(sub_cluster == cc) %>% select(all_of(descriptivo_grupal))
      bootnet::estimateNetwork(dat_cc, default = "IsingFit", tuning = 0.25)
    })
    names(models) <- paste0("c", clusters)

    # 4.a Heatmap de activación con etiquetas coloreadas por asignación a subred (k colores)
    assign_df <- assign_items_to_subnet(models)
    # Asegura el nombre de columna consistente para join posteriores
    if ("node" %in% names(assign_df) && !("Variable" %in% names(assign_df))) {
      assign_df <- dplyr::rename(assign_df, Variable = node)
    }
    plot_activation_colored_by_assignment(
      base_trim = base_trim, items = descriptivo_grupal, assign_df = assign_df,
      plots_dir = plots_dir, w = w, k_sel = k_sel, start_sel = start_sel, end_sel = end_sel
    )

    # 5. Un único PDF con k componentes (k subredes), mismo layout
    plot_k_networks_onepage(
      models = models,
      plot_path = file.path(plots_dir, sprintf("networks_k%d_w%d_start%d_end%d.pdf", k_sel, w, start_sel, end_sel)),
      title = sprintf("k=%d (w=%d) | ventana [%d,%d] | distancia=%s", k_sel, w, start_sel, end_sel, analysis_distance_method)
    )
  }

  message(sprintf("Archivos de w=%d en:\n- %s\n- %s", w, files_dir, plots_dir))
}

if (use_parallel) try(parallel::stopCluster(cl), silent = TRUE)

message("Listo. Sensibilidad por tamaños de ventana completada.")

