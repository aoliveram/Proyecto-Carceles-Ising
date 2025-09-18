# 04-veredicto-k-y-distancia.R
# Analiza resultados de k-grid y emite veredicto sobre k y distancia (binary vs jaccard)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

# Parámetros
files_dir <- "psy_net_files"
cor_threshold <- 0.7    # Debe coincidir con el usado en los heatmaps

# Carga de archivos (usa CSV si existen; puedes cambiar a RDS si prefieres)
kgrid_path      <- file.path(files_dir, "optimal_k_grid_by_window_1.csv")
opt_by_win_path <- file.path(files_dir, "optimal_k_by_window_1.csv")
opt_sum_path    <- file.path(files_dir, "optimal_k_global_summary_1.csv")

stopifnot(file.exists(kgrid_path))
kgrid <- suppressMessages(readr::read_csv(kgrid_path, show_col_types = FALSE)) %>%
  filter(!is.na(composite_score), !is.na(silhouette_score), !is.na(mean_pears_cor)) %>%
  mutate(distance = as.character(distance))

if (file.exists(opt_by_win_path)) {
  opt_by_win <- suppressMessages(readr::read_csv(opt_by_win_path, show_col_types = FALSE)) %>%
    mutate(distance = as.character(distance))
} else {
  # Reconstruye k óptimo por ventana si no existe el CSV
  opt_by_win <- kgrid %>%
    group_by(distance, start_val, end_val) %>%
    slice_min(composite_score, n = 1, with_ties = TRUE) %>%
    mutate(k_opt = k) %>%
    ungroup() %>%
    select(distance, start_val, end_val, k_opt, silhouette_score, mean_pears_cor, composite_score, cluster_sizes)
}

# 1) Resumen simple por distancia: k global candidatos
k_global_mode <- opt_by_win %>%
  count(distance, k_opt, name = "freq") %>%
  group_by(distance) %>%
  arrange(distance, desc(freq), k_opt) %>%
  slice_head(n = 1) %>%
  rename(k_mode = k_opt) %>%
  select(distance, k_mode, freq)

k_global_min_mean <- kgrid %>%
  group_by(distance, k) %>%
  summarise(mean_composite = mean(composite_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(distance) %>%
  arrange(distance, mean_composite, k) %>%
  slice_head(n = 1) %>%
  rename(k_minComp = k)

candidates <- k_global_mode %>%
  left_join(k_global_min_mean, by = "distance")

# 2) Decidir k por distancia (regla: si k_mode != k_minComp y la ganancia de composite > eps, elige k_minComp)
eps <- 0.03  # umbral de “ventaja clara” en composite
decide_k_for_distance <- function(dist_name) {
  cand <- candidates %>% filter(distance == dist_name)
  if (nrow(cand) == 0) return(NULL)
  kmode    <- cand$k_mode[1]
  kmincomp <- cand$k_minComp[1]
  # composites medios en el k de interés
  comp_kmode <- kgrid %>% filter(distance == dist_name, k == kmode) %>%
    summarise(mc = mean(composite_score, na.rm = TRUE), .groups = "drop") %>% pull(mc)
  comp_kmin  <- kgrid %>% filter(distance == dist_name, k == kmincomp) %>%
    summarise(mc = mean(composite_score, na.rm = TRUE), .groups = "drop") %>% pull(mc)
  chosen_k <- if (!is.na(comp_kmode) && !is.na(comp_kmin) && (comp_kmode - comp_kmin) > eps) kmincomp else kmode
  tibble(distance = dist_name,
         k_mode = kmode, k_minComp = kmincomp,
         mean_comp_at_k_mode = comp_kmode,
         mean_comp_at_k_min = comp_kmin,
         chosen_k = chosen_k)
}

by_distance <- bind_rows(lapply(unique(kgrid$distance), decide_k_for_distance))

# 3) Métricas de respaldo al k escogido por distancia
support_metrics <- kgrid %>%
  inner_join(by_distance %>% select(distance, chosen_k), by = c("distance" = "distance", "k" = "chosen_k")) %>%
  group_by(distance) %>%
  summarise(
    n_windows                = n(),
    mean_composite_chosen_k  = mean(composite_score, na.rm = TRUE),
    sd_composite_chosen_k    = sd(composite_score, na.rm = TRUE),
    mean_silhouette_chosen_k = mean(silhouette_score, na.rm = TRUE),
    frac_high_mean_cor       = mean(mean_pears_cor > cor_threshold, na.rm = TRUE),  # menor = mejor
    .groups = "drop"
  )

# 4) Acuerdo / desacuerdo de k_opt por ventana entre distancias
agreement <- opt_by_win %>%
  select(distance, start_val, k_opt) %>%
  distinct() %>%
  pivot_wider(names_from = distance, values_from = k_opt) %>%
  mutate(agree = ifelse(!is.na(binary) & !is.na(jaccard) & binary == jaccard, TRUE, FALSE))

agree_rate <- agreement %>%
  summarise(
    n = n(),
    agree_n = sum(agree, na.rm = TRUE),
    agree_rate = mean(agree, na.rm = TRUE)
  )

# 5) Veredicto final de distancia:
#    criterio principal: menor mean_composite_chosen_k; desempate: mayor silhouette y menor frac_high_mean_cor
ranked_dist <- support_metrics %>%
  arrange(mean_composite_chosen_k, desc(mean_silhouette_chosen_k), frac_high_mean_cor)

recommended <- ranked_dist %>% slice(1)
recommended_distance <- recommended$distance[1]
recommended_k <- by_distance %>% filter(distance == recommended_distance) %>% pull(chosen_k)

# 6) Salida a consola
cat("\n===== Veredicto por distancia =====\n")
print(by_distance)
cat("\n===== Métricas de respaldo (en k elegido por distancia) =====\n")
print(support_metrics)
cat("\nAcuerdo binary vs jaccard (k_opt por ventana):\n")
print(agree_rate)

cat("\n===== Veredicto FINAL =====\n")
cat(sprintf("Distancia recomendada: %s\n", recommended_distance))
cat(sprintf("k recomendado: %s\n", recommended_k))
cat(sprintf("mean_composite: %.4f | silhouette: %.4f | frac_high_mean_cor(>%.2f): %.3f\n",
            recommended$mean_composite_chosen_k,
            recommended$mean_silhouette_chosen_k,
            cor_threshold,
            recommended$frac_high_mean_cor))

# (Opcional) pequeño vistazo a cómo varía composite con k por distancia
# Descomenta para un resumen extra en consola:
# kgrid %>%
#   group_by(distance, k) %>%
#   summarise(mean_comp = mean(composite_score, na.rm = TRUE),
#             mean_sil  = mean(silhouette_score, na.rm = TRUE),
#             frac_hi   = mean(mean_pears_cor > cor_threshold, na.rm = TRUE),
#             .groups = "drop") %>%
#   arrange(distance, mean_comp, k) %>%
#   print(n = 100)