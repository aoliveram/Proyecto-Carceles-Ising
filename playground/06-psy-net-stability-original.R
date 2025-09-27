# ======================================================================
# 06-psy-net-stability.R
# Estabilidad de "grupos de variables" con partición de personas fija
# - Recorre TODAS las ventanas dentro de [P10,P90] del puntaje IGI
# - w fijo (e.g., 5), k en 2:4
# - B réplicas bootstrap por subgrupo (paralelo FORK por réplica)
# - Métricas: p_hat(v), co-asignación, ARI medio + IC (por ventana)
# - Plots: UN PDF por (w,k), con k filas (subredes) x columnas=ventanas (top-25 por celda)
# - Salida: psy-net-stability/
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(bootnet)
  library(cluster)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(mclust)      # adjustedRandIndex
  library(stringr)
})

set.seed(123)

# ----------------------- Configuración --------------------------------
input_csv  <- "psy_net_files/base_igi_bin_1.csv"
out_dir    <- "psy-net-stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Ítems (excluir HD6 y PAR24)
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

# Parámetros del estudio
w                 <- 5
k_values          <- 2:4
B                 <- 500               # nº de réplicas bootstrap (ajusta si hace falta)
min_cluster_size  <- 100
analysis_distance <- "binary"         # "jaccard" o "binary"
trim_q            <- c(0.25, 0.75)     # percentiles para recorte de IGI

# Paralelo (sin cambios): FORK para réplicas y ensamblaje; ventanas en serie
n_cores           <- 8
cluster_type      <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")

# ------------------------- Utilidades ---------------------------------
compute_distance_binary <- function(mat_bin) {
  M <- as.matrix(mat_bin); M[is.na(M)] <- 0; M <- (M > 0) + 0L
  p <- ncol(M); rs <- rowSums(M)
  inter <- M %*% t(M)
  hamming <- outer(rs, rs, "+") - 2*inter
  D <- hamming / max(p, 1L); diag(D) <- 0; as.matrix(D)
}
compute_distance_jaccard <- function(mat_bin) {
  M <- as.matrix(mat_bin); M[is.na(M)] <- 0; M <- (M > 0) + 0L
  rs <- rowSums(M); inter <- M %*% t(M)
  union <- outer(rs, rs, "+") - inter
  sim <- inter / pmax(union, 1); sim[union == 0] <- 1
  D <- 1 - sim; diag(D) <- 0; as.matrix(D)
}
compute_distance_single <- function(mat_bin, method = c("binary","jaccard")) {
  method <- match.arg(method)
  if (method == "binary") compute_distance_binary(mat_bin) else compute_distance_jaccard(mat_bin)
}
upper_tri_vec <- function(W) {
  ut <- upper.tri(W, diag = FALSE)
  as.numeric(W[ut])
}
node_strength_abs <- function(W) rowSums(abs(W), na.rm = TRUE)

assign_items_by_strength <- function(models_list) {
  # models_list: lista de objetos estimateNetwork (uno por subgrupo)
  S <- sapply(models_list, function(m) node_strength_abs(m$graph))
  if (is.null(dim(S))) S <- matrix(S, ncol = 1)
  rn <- rownames(models_list[[1]]$graph)
  if (is.null(rn)) rn <- paste0("V", seq_len(nrow(models_list[[1]]$graph)))
  asign <- max.col(S, ties.method = "first") # argmax por fila
  tibble(Variable = rn, group = asign)
}

# Drop columns with zero variance (all same value or all NA)
drop_zero_variance <- function(df) {
  zs <- vapply(df, function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))
  list(df = df[ , !zs, drop = FALSE], dropped = names(df)[zs])
}

# Compute node strengths for a graph, filling missing items with 0
strength_from_graph <- function(W, all_items) {
  s <- rowSums(abs(W), na.rm = TRUE)
  out <- setNames(rep(0, length(all_items)), all_items)
  nm <- names(s); if (is.null(nm)) nm <- rownames(W)
  out[nm] <- s
  out
}

entropy_normalized <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log(p)) / log(length(counts))
}

save_csv <- function(df, name) {
  readr::write_csv(df, file.path(out_dir, name))
}

# ----------------------- Carga y trimming ------------------------------
stopifnot(file.exists(input_csv))
base <- read.csv(input_csv, stringsAsFactors = FALSE)
stopifnot(all(c(descriptivo_grupal, "puntaje_total") %in% names(base)))

# Recorte por percentiles
q_lo <- as.numeric(quantile(base$puntaje_total, probs = trim_q[1], na.rm = TRUE))
q_hi <- as.numeric(quantile(base$puntaje_total, probs = trim_q[2], na.rm = TRUE))
# Asegurar enteros y rango válido para ventanas deslizantes
q_lo <- ceiling(q_lo)
q_hi <- floor(q_hi)
base_trim <- base %>% filter(puntaje_total >= q_lo, puntaje_total <= q_hi)

lo <- min(base_trim$puntaje_total, na.rm = TRUE)
hi <- max(base_trim$puntaje_total, na.rm = TRUE)
message(sprintf("Rango IGI recortado: [%d, %d] (P10-P90)", lo, hi))

# --------------------- Estabilidad por (w, k) --------------------------
stability_for_k <- function(k) {
  # Ventanas deslizantes (en serie, como pediste)
  start_vals <- seq(from = lo, to = hi - w + 1, by = 1)
  
  # Contenedores agregados (por (w,k), todas las ventanas)
  agg_varstab <- list()
  agg_partstab <- list()
  agg_coassign <- list()
  agg_edges_all <- list()
  
  # Bucle SECUENCIAL por ventana
  for (sv in start_vals) {
    ev <- sv + w - 1
    window_id <- sprintf("[%d,%d]", sv, ev)
    
    sub <- base_trim %>% dplyr::filter(puntaje_total >= sv, puntaje_total <= ev)
    if (nrow(sub) < min_cluster_size) {
      message(sprintf("[w=%d, k=%d, %s] Saltada: n<min_cluster_size", w, k, window_id))
      next
    }
    
    # --- Cluster de personas (distancia fija) ---
    M  <- as.matrix(sub[, descriptivo_grupal, drop = FALSE])
    D  <- compute_distance_single(M, analysis_distance)
    clh <- hclust(as.dist(D), method = "ward.D2")
    sub$sub_cluster <- cutree(clh, k = k)
    sizes <- table(sub$sub_cluster)
    if (any(sizes < min_cluster_size)) {
      message(sprintf("[w=%d, k=%d, %s] Saltada: cluster pequeño", w, k, window_id))
      next
    }
    clusters <- sort(unique(sub$sub_cluster))
    
    # --- Redes de referencia por subgrupo ---
    models_ref <- lapply(clusters, function(cc){
      dat_cc <- sub %>% dplyr::filter(sub_cluster == cc) %>% dplyr::select(dplyr::all_of(descriptivo_grupal))
      bootnet::estimateNetwork(dat_cc, default = "IsingFit", tuning = 0.25)
    })
    
    # --- Bootstrap por réplica (FORK) ---
    items <- colnames(sub[, descriptivo_grupal, drop = FALSE])
    k_sub  <- length(clusters)
    
    cl <- makeCluster(n_cores, type = cluster_type)
    registerDoParallel(cl)
    time_init <- Sys.time()
    
    rep_list <- foreach(b = 1:B,
                        .packages = c("bootnet","dplyr"),
                        .export   = c("sub","clusters","descriptivo_grupal",
                                      "drop_zero_variance","strength_from_graph")) %dopar% {
                                        W_list  <- vector("list", k_sub)
                                        S_vecs  <- vector("list", k_sub)
                                        for (i in seq_along(clusters)) {
                                          cc <- clusters[i]
                                          sub_cc <- dplyr::filter(sub, sub_cluster == cc)
                                          idx <- sample.int(nrow(sub_cc), size = nrow(sub_cc), replace = TRUE)
                                          dat_b <- sub_cc[idx, descriptivo_grupal, drop = FALSE]
                                          dv <- drop_zero_variance(dat_b)
                                          dat_est <- dv$df
                                          fit_i <- bootnet::estimateNetwork(dat_est, default = "IsingFit", tuning = 0.25)
                                          W_list[[i]] <- fit_i$graph
                                          S_vecs[[i]] <- strength_from_graph(fit_i$graph, all_items = colnames(dat_b))
                                        }
                                        list(W_list = W_list, S_list = S_vecs)
                                      }
    
    time_fin <- Sys.time()
    message(sprintf("[w=%d, k=%d, %s] Réplicas estimadas (B=%d) en %s",
                    w, k, window_id, B, as.character(difftime(time_fin, time_init))))
    try(parallel::stopCluster(cl), silent = TRUE)
    
    # --- Fuerzas (B x p) y grafos por subred ---
    S_list <- vector("list", length = length(clusters))
    W_boot <- vector("list", length = length(clusters))
    names(S_list) <- names(W_boot) <- paste0("subred_", seq_along(clusters))
    
    for (i in seq_along(clusters)) {
      S_mat <- t(vapply(rep_list, function(rr) rr$S_list[[i]][items], numeric(length(items))))
      colnames(S_mat) <- items
      W_boot[[i]] <- lapply(rep_list, function(rr) rr$W_list[[i]])
      S_list[[i]] <- S_mat
    }
    
    # --- Asignación de referencia (argmax de fuerza) ---
    assign_ref <- assign_items_by_strength(models_ref)
    ref_labels <- assign_ref$group
    names(ref_labels) <- assign_ref$Variable
    items_ord <- assign_ref$Variable
    p <- length(items_ord)
    
    # --- Ensamblaje de asignaciones (FORK, por réplica) ---
    cl <- makeCluster(n_cores, type = cluster_type)
    registerDoParallel(cl)
    time_init <- Sys.time()
    
    assignments_mat <- foreach(b = 1:B, .combine = "cbind",
                               .packages = c("dplyr"),
                               .export   = c("S_list","items_ord","k")) %dopar% {
                                 Sb <- vapply(S_list, function(S) S[b, ], numeric(length(items_ord)))
                                 as.integer(max.col(Sb, ties.method = "first"))
                               }
    
    time_fin <- Sys.time()
    message(sprintf("[w=%d, k=%d, %s] Ensamblaje de asignaciones (B=%d) en %s",
                    w, k, window_id, B, as.character(difftime(time_fin, time_init))))
    try(parallel::stopCluster(cl), silent = TRUE)
    colnames(assignments_mat) <- paste0("b", seq_len(B))
    rownames(assignments_mat) <- items_ord
    
    # --- Métricas por ventana ---
    p_hat <- rowMeans(assignments_mat == ref_labels)
    ent <- apply(assignments_mat, 1, function(v) {
      tabulate(v, nbins = k) |> entropy_normalized()
    })
    n_groups_var <- apply(assignments_mat, 1, function(v) length(unique(v)))
    
    coassign <- matrix(0, nrow = p, ncol = p, dimnames = list(items_ord, items_ord))
    for (b in seq_len(B)) {
      g <- assignments_mat[, b]
      coassign <- coassign + outer(g, g, FUN = "==")
    }
    coassign <- coassign / B
    
    ARI_vec <- sapply(seq_len(B), function(b) adjustedRandIndex(ref_labels, assignments_mat[, b]))
    ARI_mean <- mean(ARI_vec)
    ARI_ci   <- as.numeric(quantile(ARI_vec, probs = c(0.025, 0.975), names = FALSE))
    
    # --- Acumular CSV (por ventana) ---
    agg_varstab[[window_id]] <- tibble(
      window_start = sv, window_end = ev, window_id = window_id,
      Variable     = items_ord,
      ref_group    = ref_labels,
      p_hat        = p_hat,
      entropy      = ent,
      n_groups     = n_groups_var
    )
    
    agg_partstab[[window_id]] <- tibble(
      w = w, k = k, B = B,
      window_start = sv, window_end = ev, window_id = window_id,
      ARI_mean = ARI_mean,
      ARI_q025 = ARI_ci[1],
      ARI_q975 = ARI_ci[2]
    )
    
    coassign_df <- as.data.frame(coassign) %>%
      mutate(window_start = sv, window_end = ev, window_id = window_id,
             .before = 1) %>%
      mutate(Variable = rownames(coassign), .before = 1)
    agg_coassign[[window_id]] <- coassign_df
    
    # --- Estabilidad de aristas (top-25 por subred) para el GRAN PDF ---
    edge_panel_list <- list()
    for (i in seq_along(clusters)) {
      W_ref <- models_ref[[i]]$graph
      rn <- rownames(W_ref); if (is.null(rn)) rn <- paste0("V", seq_len(nrow(W_ref)))
      ut <- upper.tri(W_ref, diag = FALSE)
      idx_pairs <- which(ut, arr.ind = TRUE)
      edges_ref <- tibble(
        i = idx_pairs[,1],
        j = idx_pairs[,2],
        n1 = rn[idx_pairs[,1]],
        n2 = rn[idx_pairs[,2]],
        w_ref = W_ref[ut]
      ) %>%
        mutate(n1s = pmin(n1, n2), n2s = pmax(n1, n2),
               w_ref_abs = abs(w_ref)) %>%
        arrange(desc(w_ref_abs)) %>%
        slice(1:min(25, n()))  # top-25 por |peso_ref|
      
      Ws_i <- W_boot[[i]]
      # Cuantiles ignorando réplicas donde la arista no existe
      q_na <- function(x, prob) {
        if (all(is.na(x))) return(NA_real_)
        as.numeric(stats::quantile(x, prob, na.rm = TRUE, names = FALSE))
      }
      get_quant_by_names <- function(n1, n2) {
        w_b <- vapply(Ws_i, function(Wb) {
          rn <- rownames(Wb); cn <- colnames(Wb)
          if (is.null(rn)) rn <- colnames(Wb)
          if (is.null(cn)) cn <- rownames(Wb)
          if (!(n1 %in% rn && n2 %in% cn)) return(NA_real_)
          as.numeric(Wb[n1, n2])
        }, numeric(1))
        c(lo = q_na(w_b, 0.025),
          hi = q_na(w_b, 0.975))
      }
      Q <- t(mapply(get_quant_by_names, edges_ref$n1, edges_ref$n2))
      edges_plot <- dplyr::bind_cols(edges_ref, as.data.frame(Q))
      edges_plot$subred <- paste0("Subred ", i)
      edges_plot$window_id <- window_id
      edges_plot$window_start <- sv
      edges_plot$window_end   <- ev
      edge_panel_list[[i]] <- edges_plot
    }
    edges_all <- bind_rows(edge_panel_list)
    agg_edges_all[[window_id]] <- edges_all
  } # fin loop ventanas
  
  # ------------------ Guardar CSV agregados por (w,k) ------------------
  if (length(agg_varstab)) {
    varstab_all <- bind_rows(agg_varstab)
    save_csv(varstab_all, sprintf("variable_stability_w%d_k%d.csv", w, k))
  }
  if (length(agg_partstab)) {
    partstab_all <- bind_rows(agg_partstab)
    save_csv(partstab_all, sprintf("partition_stability_w%d_k%d.csv", w, k))
  }
  if (length(agg_coassign)) {
    coassign_all <- bind_rows(agg_coassign)
    save_csv(coassign_all, sprintf("coassign_w%d_k%d.csv", w, k))
  }
  
  # ------------------ GRAN PDF: todas las ventanas por (w,k) -----------
  if (length(agg_edges_all)) {
    edges_all_k <- bind_rows(agg_edges_all)
    
    # Ordenar ventanas en el eje-x por window_start
    win_levels <- edges_all_k %>%
      distinct(window_id, window_start) %>%
      arrange(window_start) %>%
      pull(window_id)
    edges_all_k$window_id <- factor(edges_all_k$window_id, levels = win_levels)
    
    # Paginación automática si hay muchas ventanas
    ncols_per_page <- 12
    pages <- split(win_levels, ceiling(seq_along(win_levels) / ncols_per_page))
    
    pdf_path <- file.path(out_dir, sprintf("edge_stability_w%d_k%d.pdf", w, k))
    grDevices::pdf(pdf_path, width = max(10, 2.6 * min(ncols_per_page, length(win_levels))), height = 10)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    
    for (pg in seq_along(pages)) {
      cols <- pages[[pg]]
      dat_pg <- edges_all_k %>% filter(window_id %in% cols)
      
      p_edges <- ggplot(dat_pg, aes(x = w_ref, y = reorder(paste(n1, "—", n2), w_ref_abs))) +
        geom_segment(aes(x = lo, xend = hi, yend = reorder(paste(n1, "—", n2), w_ref_abs)),
                     alpha = 0.6, na.rm = TRUE) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_point(color = "red", size = 0.8) +
        facet_grid(rows = vars(subred), cols = vars(window_id), scales = "free_y") +
        labs(x = "Peso de arista (ref y CI 2.5–97.5% bootstrap)",
             y = "Aristas (top-25 por |peso_ref|)",
             title = sprintf("Estabilidad de aristas — w=%d, k=%d (ventanas en [%d,%d], trim %d–%d%%)",
                             w, k, lo, hi, round(trim_q[1]*100), round(trim_q[2]*100))) +
        theme_minimal(base_size = 9) +
        theme(strip.text.y = element_text(face = "bold"))
      
      print(p_edges)
    }
    grDevices::dev.off()
    message(sprintf("[w=%d, k=%d] PDF generado: %s", w, k, normalizePath(pdf_path)))
  }
  
  invisible(TRUE)
}

# ------------------------- Ejecutar -------------------------------------------
time_init_all <- Sys.time()
for (k in k_values) {
  try(stability_for_k(k), silent = FALSE)
}
time_fin_all <- Sys.time()
message(sprintf("Finalizado en %s", as.character(difftime(time_fin_all, time_init_all))))
