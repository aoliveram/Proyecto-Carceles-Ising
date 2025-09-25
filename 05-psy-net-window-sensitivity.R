# =============================================================================
# 05-psy-net-window-sensitivity.R
# Sensibilidad al tamaño de ventana: window = 3..10
# - Excluye ítems de baja varianza: HD6 y PAR24
# - Para cada window:
#   * Barrido ventana × k (k=2..6) con distancias binary y jaccard
#   * Guarda CSV/RDS y heatmaps por distancia en carpetas específicas
#   * Selecciona una ventana representativa por k (mínimo mean_pears_cor con distancia de análisis)
#   * Genera:
#     - Un heatmap de activación (P10–P90) con etiquetas de ítems coloreadas por asignación a la subred (k-colores)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(bootnet)
  library(cluster)
  library(foreach)
  library(doParallel)
  library(scales)
  library(grid)
})

set.seed(123)

# Entradas / Salidas base
input_csv      <- "psy_net_files/base_igi_bin_1.csv"
root_out_dir <- "psy-net-window-sensitivity"  

# Parámetros
window_sizes          <- 3:10
k_values              <- 2:6
min_cluster_size      <- 100
use_trim              <- TRUE
trim_q                <- c(0.00, 1.00)
analysis_distance_method <- "jaccard"
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
      mean_pears_cor = NA_real_,
      mean_frobenius = NA_real_,
      cluster_sizes = NA_character_,
      error_message = "Insufficient observations",
      cluster_variable_sizes = NA_character_
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
        mean_pears_cor = NA_real_,
        mean_frobenius = NA_real_,
        cluster_sizes = paste(sizes, collapse=", "),
        error_message = "Cluster too small",
        cluster_variable_sizes = NA_character_
      ))
    }
    # Redes por subcluster (Ising por cluster de personas)
    clusters <- sort(unique(sub$sub_cluster))
    models <- lapply(clusters, function(cc) {
      dat_cc <- sub %>% filter(sub_cluster == cc) %>% select(all_of(descriptivo_grupal))
      bootnet::estimateNetwork(dat_cc, default = "IsingFit", tuning = 0.25)
    })
    pairs <- combn(seq_along(models), 2, simplify = FALSE)
    # Pearson similarity across networks (lower = more different)
    corrs <- vapply(pairs, function(p) {
      cor(as.vector(models[[p[1]]]$graph), as.vector(models[[p[2]]]$graph), use = "complete.obs")
    }, numeric(1))
    mean_cor <- mean(corrs)

    # Frobenius distance across networks (higher = more different)
    frobs <- vapply(pairs, function(p) {
      W1 <- models[[p[1]]]$graph; W2 <- models[[p[2]]]$graph
      sqrt(sum((W1 - W2)^2, na.rm = TRUE))
    }, numeric(1))
    mean_frob <- mean(frobs)

    # Tamaño de grupos de variables (asignación por fuerza absoluta)
    assign_df <- assign_items_to_subnet(models)
    var_counts <- sapply(seq_along(models), function(i) sum(assign_df$asignacion == i, na.rm = TRUE))
    var_counts_str <- paste(as.integer(var_counts), collapse = ", ")

    tibble(
      start_val = start_val, end_val = end_val, k = k, distance = dist_name,
      mean_pears_cor = mean_cor,
      mean_frobenius = mean_frob,
      cluster_sizes = paste(sizes, collapse = ", "),
      cluster_variable_sizes = var_counts_str,
      error_message = NA_character_
    )
  }

  out_bin <- map_dfr(k_values, ~do_one(.x, D_bin, "binary"))
  out_jac <- map_dfr(k_values, ~do_one(.x, D_jac, "jaccard"))
  bind_rows(out_bin, out_jac)
}

# Ploteo: heatmaps por distancia
plot_heatmaps_distance <- function(kgrid_w, plots_dir, w) {
  for (dist_name in c("binary","jaccard")) {
    kd <- kgrid_w %>% filter(distance == dist_name)
    if (nrow(kd) == 0) next
    overlay_pts <- kd %>% filter(!is.na(mean_pears_cor), mean_pears_cor > cor_threshold)
    window_size <- unique(kd$end_val - kd$start_val + 1)

    pdf(file.path(plots_dir, sprintf("optimal_k_mean_cor_w%d_%s.pdf", w, dist_name)), width = 10, height = 7)
    print(
      ggplot(kd, aes(x = start_val, y = k, fill = mean_pears_cor)) +
        geom_tile() +
        geom_point(data = overlay_pts, aes(x = start_val, y = k),
                   inherit.aes = FALSE, shape = 4, color = "red", stroke = 0.9, size = 1.8) +
        scale_fill_viridis_c(option = "plasma", direction = -1) +
        labs(x = "Window Start (IGI)", y = "k",
             fill = "Mean Corr (lower=better)",
             title = sprintf("Mean Pearson Correlation - %s | Window size = %d", dist_name, window_size)) +
        theme_minimal()
    )
    dev.off()

    # Frobenius heatmap
    pdf(file.path(plots_dir, sprintf("optimal_k_mean_frobenius_w%d_%s.pdf", w, dist_name)), width = 10, height = 7)
    print(
      ggplot(kd, aes(x = start_val, y = k, fill = mean_frobenius)) +
        geom_tile() +
        scale_fill_viridis_c(option = "plasma", direction = 1) +
        labs(x = "Window Start (IGI)", y = "k",
             fill = "Mean Frobenius (higher=worse)",
             title = sprintf("Mean Frobenius Distance - %s | Window size = %d", dist_name, window_size)) +
        theme_minimal()
    )
    dev.off()
  }
}

# -----------------------------------------------------------------------------#
# Loop principal por tamaños de ventana
# -----------------------------------------------------------------------------#

if (use_parallel) {
  cl <- makeCluster(n_cores, type = cluster_type)
  registerDoParallel(cl)
}

time_init_all <- Sys.time()
for (w in window_sizes) {
  message(sprintf("=== window = %d ===", w))
  out_dir <- root_out_dir
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  files_dir <- out_dir
  plots_dir <- out_dir

  start_vals <- seq(from = lo, to = hi - w + 1, by = 1)

  # k-grid (por distancia) para todas las ventanas de este tamaño
  time_init <- Sys.time()
  kgrid_w <- foreach(v = start_vals, .combine = "rbind",
                      .packages = c("dplyr","tibble","tidyr","purrr","bootnet","cluster"),
                      .export = c("base_trim","descriptivo_grupal","k_values","min_cluster_size",
                                  "compute_distance_binary","compute_distance_jaccard","process_window_kgrid")) %dopar% {
    process_window_kgrid(v, w)
  }
  time_fin <- Sys.time()
  message(sprintf("k-grid (w=%d) listo en %s", w, as.character(difftime(time_fin, time_init))))

  kgrid_w <- as_tibble(kgrid_w)
  safe_write(kgrid_w, file.path(files_dir, sprintf("optimal_k_grid_by_window_w%d", w)))

  # k óptimo por ventana y k>2 (criterio: mean_pears_cor más bajo); además marcar si el mejor absoluto habría sido k=2
  opt_by_win <- kgrid_w %>%
    group_by(distance, start_val, end_val) %>%
    group_modify(~{
      cand_all <- dplyr::filter(.x, !is.na(mean_pears_cor))
      if (nrow(cand_all) == 0) {
        return(tibble::tibble(
          k_opt = NA_integer_,
          best_is_k2 = NA,
          mean_pears_cor = NA_real_,
          mean_frobenius = NA_real_,
          cluster_sizes = NA_character_,
          cluster_variable_sizes = NA_character_
        ))
      }
      # Mejor absoluto (incluye k=2) por mean_pears_cor (más bajo = mejor)
      best_overall_row <- cand_all %>%
        dplyr::slice_min(order_by = mean_pears_cor, n = 1, with_ties = TRUE) %>%
        dplyr::slice(1)
      best_is_k2 <- best_overall_row$k[1] == 2

      # Mejor entre k > 2
      cand_gt2 <- cand_all %>% dplyr::filter(k > 2)
      if (nrow(cand_gt2) == 0) {
        k_opt <- NA_integer_
        mcor <- NA_real_; mfrob <- NA_real_; cs <- NA_character_; cvs <- NA_character_
      } else {
        best_gt2_row <- cand_gt2 %>%
          dplyr::slice_min(order_by = mean_pears_cor, n = 1, with_ties = TRUE) %>%
          dplyr::slice(1)
        k_opt <- best_gt2_row$k[1]
        mcor <- best_gt2_row$mean_pears_cor[1]
        mfrob <- best_gt2_row$mean_frobenius[1]
        cs <- best_gt2_row$cluster_sizes[1]
        cvs <- best_gt2_row$cluster_variable_sizes[1]
      }

      tibble::tibble(
        k_opt = k_opt,
        best_is_k2 = best_is_k2,
        mean_pears_cor = mcor,
        mean_frobenius = mfrob,
        cluster_sizes = cs,
        cluster_variable_sizes = cvs
      )
    }) %>%
    ungroup()
  safe_write(opt_by_win, file.path(files_dir, sprintf("optimal_k_by_window_w%d", w)))

  # Heatmaps por distancia
  plot_heatmaps_distance(kgrid_w, plots_dir, w)

  # Para los plots avanzados (heatmap coloreado por asignación y redes por k),
  # elegimos una ventana "representativa" por k con la distancia de análisis:
  # - Tomamos la ventana con mean_pears_cor mínimo para (distance=analysis_distance_method, k fijo).
  kd_sel <- kgrid_w %>% filter(distance == analysis_distance_method)

  # Bump chart por window size
  # 1) Recolecta asignaciones Variable->cluster_raw para cada k (usando la ventana representativa por k)
  assignments_list <- list()
  for (k_sel in k_values) {
    cand <- kd_sel %>% filter(k == k_sel, !is.na(mean_pears_cor))
    if (nrow(cand) == 0) next
    row_best <- cand %>% slice_min(order_by = mean_pears_cor, n = 1, with_ties = TRUE) %>% slice(1)
    start_sel <- row_best$start_val[1]; end_sel <- row_best$end_val[1]

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

    assign_df <- assign_items_to_subnet(models)
    if ("node" %in% names(assign_df) && !("Variable" %in% names(assign_df))) {
      assign_df <- dplyr::rename(assign_df, Variable = node)
    }
    assignments_list[[as.character(k_sel)]] <- assign_df %>%
      dplyr::transmute(Variable, k = k_sel, cluster_raw = asignacion)
  }

  # 2) Alinea k clusters (herederos por máxima coincidencia)
  align_colors <- function(df_list) {
    ks <- sort(as.integer(names(df_list)))
    if (length(ks) == 0) return(NULL)
    # mapa: lista por k con vector nombres=cluster_raw, valores=color_id
    color_maps <- list()

    k0 <- ks[1]
    # colores iniciales: 1..k0
    init_clusters <- sort(unique(df_list[[as.character(k0)]]$cluster_raw))
    color_maps[[as.character(k0)]] <- setNames(seq_along(init_clusters), init_clusters)
    next_color_id <- length(init_clusters) + 1

    for (idx in 2:length(ks)) {
      k_prev <- ks[idx-1]; k_cur <- ks[idx]
      prev_df <- df_list[[as.character(k_prev)]]
      cur_df  <- df_list[[as.character(k_cur)]]
      prev_map <- color_maps[[as.character(k_prev)]]

      prev_clusters <- sort(unique(prev_df$cluster_raw))
      cur_clusters  <- sort(unique(cur_df$cluster_raw))

      # tabla de coincidencias |prev ∩ cur|
      tab <- table(factor(prev_df$cluster_raw, levels = prev_clusters),
                   factor(cur_df$cluster_raw,  levels = cur_clusters))[ , , drop = FALSE]

      # greedy matching por mayor solapamiento, uno-a-uno
      assigned_prev <- setNames(rep(FALSE, length(prev_clusters)), prev_clusters)
      assigned_cur  <- setNames(rep(FALSE, length(cur_clusters)),  cur_clusters)
      cur_color_map <- setNames(rep(NA_integer_, length(cur_clusters)), cur_clusters)

      overlaps <- as.data.frame(as.table(tab), stringsAsFactors = FALSE)
      names(overlaps) <- c("prev","cur","n")
      overlaps <- overlaps[order(-overlaps$n), ]
      for (r in seq_len(nrow(overlaps))) {
        pv <- as.character(overlaps$prev[r])
        cu <- as.character(overlaps$cur[r])
        if (overlaps$n[r] <= 0) break
        if (!assigned_prev[pv] && !assigned_cur[cu]) {
          cur_color_map[cu] <- prev_map[[pv]]
          assigned_prev[pv] <- TRUE
          assigned_cur[cu]  <- TRUE
        }
      }
      # clusters actuales no asignados reciben colores nuevos
      for (cu in names(cur_color_map)) {
        if (is.na(cur_color_map[cu])) {
          cur_color_map[cu] <- next_color_id
          next_color_id <- next_color_id + 1
        }
      }
      color_maps[[as.character(k_cur)]] <- cur_color_map
    }

    # Construye salida larga Variable-k-color_id
    out <- dplyr::bind_rows(lapply(ks, function(kk){
      dkk <- df_list[[as.character(kk)]]
      cmap <- color_maps[[as.character(kk)]]
      dkk$color_id <- unname(cmap[as.character(dkk$cluster_raw)])
      dkk
    }))
    out
  }

  assignments_long <- align_colors(assignments_list)
  if (!is.null(assignments_long) && nrow(assignments_long) > 0) {
    # Orden estable de variables (todas en el eje Y)
    var_levels <- sort(unique(assignments_long$Variable))
    grid_df <- assignments_long %>%
      dplyr::mutate(Variable = factor(Variable, levels = var_levels))

    # Paleta estable por color_id
    n_colors <- length(unique(grid_df$color_id))
    pal <- scales::hue_pal()(max(n_colors, 1))

    # Alto del PDF en función del # de variables (con límites razonables)
    n_vars <- length(var_levels)
    pdf_height <- max(8, min(0.28 * n_vars, 40))

    p_grid <- ggplot(grid_df, aes(x = k, y = Variable, fill = factor(color_id))) +
      geom_tile(width = 0.85, height = 0.85, show.legend = FALSE) +
      scale_fill_manual(values = pal) +
      scale_x_continuous(breaks = sort(unique(grid_df$k))) +
      scale_y_discrete(limits = var_levels) +
      labs(x = "k (número de subredes)", y = "Variable",
          title = sprintf("Asignación de variables por k (node assignment by absolute strength) — w=%d, dist=%s", 
                          w, analysis_distance_method)) +
      theme_minimal()

    pdf(file.path(plots_dir, sprintf("clusters_bump_w%d_%s.pdf", w, analysis_distance_method)),
        width = 14, height = pdf_height)
    print(p_grid)
    dev.off()
  }

  message(sprintf("Archivos de w=%d en: %s", w, out_dir))
}

if (use_parallel) try(parallel::stopCluster(cl), silent = TRUE)

time_fin_all <- Sys.time()
time_tot_all <- time_fin_all - time_init_all
print(time_tot_all)