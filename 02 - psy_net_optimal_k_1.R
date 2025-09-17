# Vemos la correlación de Pearson para una separación de la red en k=3 clusters.
# 1) Analiza ventanas deslizantes de 10 puntos del puntaje IGI con k=3 subclusters.
# 2) Para cada ventana: clusteriza sujetos (distancia binaria con stats::dist y hclust Ward.D2),
#    estima redes de Ising por subcluster (IsingFit) y compara redes entre subclusters con
#    correlación de Pearson de las matrices de pesos (menor = más distintas).
# 3) Calcula silhouette promedio por ventana y un puntaje compuesto (1 - silhouette) + media de
#    correlaciones entre redes (usado como referencia diagnóstica).
# 4) Evalúa asociación entre subclusters y variables externas (SEXO, salud mental, edad, estado civil)
#    mediante tablas de contingencia y pruebas de chi-cuadrado.
# 5) Exporta resultados a CSV (cor_results_sliding_*.csv, chi_results_sliding_*.csv) y genera gráficos
#    de evolución y heatmaps en psy_net_plots/.
# 6) Repite el mismo flujo para dos esquemas de binarización: {3,2,1},{0}->{0},{1} y {3,2},{0,1}->{0},{1}.
# 7) Notas: (a) "binary" no es Jaccard; si se requiere Jaccard usar vegan::vegdist(..., method="jaccard", binary=TRUE).
#           (b) El bloque final de análisis por rangos de riesgo depende de objetos definidos en otros scripts.

# Archivos que genera (# = 1,2):
# - psy_net_files/cor_results_sliding_#.csv
# - psy_net_files/chi_results_sliding_#.csv
# - psy_net_plots/cor_evolucion_#.png

library(bootnet)
library(qgraph)
library(dplyr)
library(tidyr)
#install.packages("EGAnet") # Problems on CHPC
#library(EGAnet)

# Este código debe tener una parte para estimar el óptimo en una parte del loop.
# Ese código es optimal_k_estimate.R
# Todo lo que está antes es simplemente para cargar los datos y procesarlos.
# Esto ya no sería necesario con:

# Cargamos datos procesados

#base_igi <- read.csv("psy_net_files/base_igi.csv")

base_igi_bin_1 <- read.csv("psy_net_files/base_igi_bin_1.csv")
base_igi_bin_2 <- read.csv("psy_net_files/base_igi_bin_2.csv")

base_igi_bin_ising_1 <- read.csv("psy_net_files/base_igi_bin_ising_1.csv")
base_igi_bin_ising_2 <- read.csv("psy_net_files/base_igi_bin_ising_2.csv")

verificar_binario <- function(data) {
  # Vemos que solo hayan 0's y 1's
  columnas_invalidas <- data %>%
    summarise(across(everything(), ~ list(setdiff(unique(.), c(0, 1))))) %>%
    pivot_longer(everything(), names_to = "columna", values_to = "valores_invalidos") %>%
    filter(lengths(valores_invalidos) > 0)  # Filtrar columnas con valores fuera de {0,1}
  
  if (nrow(columnas_invalidas) > 0) {
    # Vemos si tiene valores inválidos y cuáles son
    columnas_invalidas %>%
      rowwise() %>%
      mutate(valores_invalidos = paste(unlist(valores_invalidos), collapse = ", ")) %>%
      print()
  } else {
    print("Solo hay 0's y 1's en las columnas seleccionadas.")
  }
}

descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR24","PAR25","PAR26","PAR27",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

#verificar_binario(base_igi[,descriptivo_grupal])

verificar_binario(base_igi_bin_1[,descriptivo_grupal])
verificar_binario(base_igi_bin_2[,descriptivo_grupal])

verificar_binario(base_igi_bin_ising_1)
verificar_binario(base_igi_bin_ising_2)

# ------------------------------------------------------------------------------
# Bootnet - Ventana de 10 puntos  
# ------------------------------------------------------------------------------

# vamos al centro de nuestro objetivo: 
# veamos si para riesgos similares se producen diferencias sustantivas en las redes generadas. 
# Para ello definiremos ventanas de riesgos de 5, 10, 15 puntos y los análisis los haremos dentro de cada ventana.

# Las ventanas serán móviles. Es decir, si la ventana es de 10 puntos, se calcula para
# ventana que parte en 0 hasta 9, luego de 1 hasta 10, y así sucesivamente hasta 34 a 43.

# ---------- Formato {3,2,1},{0} -> {0},{1} ----------------------

# 1. Definimos contenedores para resultados ---
cor_results_sliding_1 <- data.frame()
chi_results_sliding_1 <- data.frame()
cross_tables_sliding_1 <- list()

# Definimos las variables externas a explorar después en las regresiones
variables_externas <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# 2. Definimos parámetros de la ventana ---
window_size <- 10 # Modificar la ventana para otros ejercicios
max_puntaje_observado <- max(base_igi_bin_1$puntaje_total)  
inicios_ventana <- 0:(max_puntaje_observado - window_size + 1)

# 3. Bucle principal por ventanas "traslapadas" ---
time_init <- Sys.time()
for (start_val in inicios_ventana) {
  
  # (a) Determinamos el rango de la ventana
  end_val <- start_val + window_size - 1  #  
  rango_puntaje <- start_val:end_val
  
  # (b) Filtramos la base por rango_puntaje  
  sub_bloque <- base_igi_bin_1 %>%
    dplyr::filter(puntaje_total %in% rango_puntaje)
  
  # Revisamos que sub_bloque tenga suficientes observaciones para clusterizar
  if (nrow(sub_bloque) < 50) {
    next # Evitemos problemas si el grupo es muy pequeño. .
  }
  
  # (c) Seleccionamos las variables de interés
  sub_bloque_items <- sub_bloque[, descriptivo_grupal]
  
  # (d) Hacemos el Cluster analysis (distancia Jaccard y Ward.D2)
  matriz_binaria <- as.matrix(sub_bloque_items)
  distancias <- dist(matriz_binaria, method = "binary") # d= 1 - shared_1's/tatal_1's
  clust_hier <- hclust(distancias, method = "ward.D2") # squared Euclidean distances. To call hclust()
  
  # (e) Cortamos en k clusters.
  # Aquí usamos un cluster jerarquico. Hay que darle una vuelta a si usamos algún 
  # mecanimos diferente para la generación de los clusters que emerga de los proios 
  # datos y no impuestos por nosotros como en este caso ¿por qué 2, 3, 4 clusters?
  sub_bloque$sub_cluster <- cutree(clust_hier, k = 3) # cambiar el número de clusters en pruebas de robustez
  
  # Revisamos los tamaños de cada cluster
  n_c1 <- sum(sub_bloque$sub_cluster == 1)
  n_c2 <- sum(sub_bloque$sub_cluster == 2)
  n_c3 <- sum(sub_bloque$sub_cluster == 3)
  
  # (f) Estimamos redes de Ising por sub_cluster (si hay suficiente n en cada cluster)
  # Nota para después en caso que haya errores para ventanas muy pequeñas: 
  # Para evitar errores, podríamos considerar después exigir un mínimo de observaciones 
  # en cada cluster. Por ejemplo agregar algo como: if (n_c1 < 10 | n_c2 < 10 | n_c3 < 10) next. 
  
  # Corremos los modelos Ising
  
  model_sub_1 <- estimateNetwork(
    data    = sub_bloque_items[sub_bloque$sub_cluster == 1, ],
    default = "IsingFit"
  )
  model_sub_2 <- estimateNetwork(
    data    = sub_bloque_items[sub_bloque$sub_cluster == 2, ],
    default = "IsingFit"
  )
  model_sub_3 <- estimateNetwork(
    data    = sub_bloque_items[sub_bloque$sub_cluster == 3, ],
    default = "IsingFit"
  )
  
  edges_1 <- model_sub_1$graph
  edges_2 <- model_sub_2$graph
  edges_3 <- model_sub_3$graph
  
  correlacion1_2 <- cor(as.vector(edges_1), as.vector(edges_2))
  correlacion1_3 <- cor(as.vector(edges_1), as.vector(edges_3))
  correlacion2_3 <- cor(as.vector(edges_2), as.vector(edges_3))
  
  # Evaluate cluster quality (silhouette width: 
  # close 1: obsv. well clusteres
  # around 0: obsv. lies between two clusters
  # negative: placed in the wrong cluster)
  silhouette_score <- mean(cluster::silhouette(
    sub_bloque$sub_cluster, distancias)[, "sil_width"])
  
  # Calculate composite score (lower is better)
  composite_score <- (1 - silhouette_score) + mean(c(correlacion1_2,correlacion1_3,correlacion2_3))
  
  # (g) Guardamos las correlaciones y tamaño de clusters en cor_results_sliding_1
  cor_results_sliding_1 <- rbind(
    cor_results_sliding_1,
    data.frame(
      min_puntaje   = start_val,
      max_puntaje   = end_val,
      cor_1_2       = round(correlacion1_2, 3),
      cor_1_3       = round(correlacion1_3, 3),
      cor_2_3       = round(correlacion2_3, 3),
      n_subcluster1 = n_c1,
      n_subcluster2 = n_c2,
      n_subcluster3 = n_c3
    )
  )
  
  # (h) Tablas de contingencia y test de chi-cuadrado para variables_externas
  for (var_externa in variables_externas) {
    
    # Construimos tabla de contingencia
    tabla <- table(sub_bloque[[var_externa]], sub_bloque$sub_cluster)
    
    # Almacenamos tabla en la lista
    nombre_tabla <- paste(var_externa, "_", start_val, "_", end_val, sep = "")
    cross_tables_sliding_1[[nombre_tabla]] <- tabla
    
    # Ejecutamos test de chi-cuadrado 
    test_result <- suppressWarnings(chisq.test(tabla))
    
    chi_results_sliding_1 <- rbind(
      chi_results_sliding_1,
      data.frame(
        min_puntaje = start_val,
        max_puntaje = end_val,
        variable    = var_externa,
        X2          = round(test_result$statistic, 2),
        df          = test_result$parameter,
        p_value     = test_result$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}
time_fin <- Sys.time()
time_total <- time_fin - time_init

write.csv(cor_results_sliding_1, "psy_net_files/cor_results_sliding_1.csv", row.names = FALSE)
write.csv(chi_results_sliding_1, "psy_net_files/chi_results_sliding_1.csv", row.names = FALSE)

cor_results_sliding_1 <- read.csv("psy_net_files/cor_results_sliding_1.csv")
chi_results_sliding_1 <- read.csv("psy_net_files/chi_results_sliding_1.csv")

# 4. Al terminar el bucle, revisamos los resultados ---

# Correlaciones entre matrices, con el rango de puntajes y n de cada subcluster
cor_results_sliding_1

# Resultados de chi cuadrado
chi_results_sliding_1

# Tablas de contingencia almacenadas
names(cross_tables_sliding_1)  # para ver qué tablas hay

# Creamos objeto para plot
cor_long_1 <- cor_results_sliding_1 %>%
  select(min_puntaje, max_puntaje, cor_1_2, cor_1_3, cor_2_3) %>%
  pivot_longer(
    cols = c("cor_1_2", "cor_1_3", "cor_2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )

# EL MISMO PLOT QUE EN PPT
png("psy_net_plots/cor_evolucion_1.png", width = 800*0.9, height = 600*0.9)
ggplot(cor_long_1, aes(x = min_puntaje, y = correlacion)) +
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
dev.off()

# Convertimos en factores
cor_long_1$clusters_comparados_factor <- factor(
  cor_long_1$clusters_comparados, 
  levels = c("cor_1_2", "cor_1_3", "cor_2_3")
)

ggplot(cor_long_1, aes(x = min_puntaje, y = clusters_comparados_factor, fill = correlacion)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "red", high = "blue", mid = "white",
    midpoint = 0.5, limit = c(0,1), space = "Lab",
    name = "Correlación"
  ) +
  theme_minimal() +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Par de subclusters",
    title = "Heatmap de correlaciones de redes entre subclusters"
  )

# ---------- Formato {3,2},{0} -> {0,1},{1} ----------------------

# 1. Definimos contenedores para resultados ---
cor_results_sliding_2 <- data.frame()
chi_results_sliding_2 <- data.frame()
cross_tables_sliding_2 <- list()

# Definimos las variables externas a explorar después en las regresiones
variables_externas <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# 2. Definimos parámetros de la ventana ---
window_size <- 10 # Modificar la ventana para otros ejercicios
max_puntaje_observado <- max(base_igi_bin_2$puntaje_total)  
inicios_ventana <- 0:(max_puntaje_observado - window_size + 1)

# 3. Bucle principal por ventanas "traslapadas" ---
time_init <- Sys.time()
for (start_val in inicios_ventana) {
  
  # (a) Determinamos el rango de la ventana
  end_val <- start_val + window_size - 1  #  
  rango_puntaje <- start_val:end_val
  
  # (b) Filtramos la base por rango_puntaje  
  sub_bloque <- base_igi_bin_2 %>%
    dplyr::filter(puntaje_total %in% rango_puntaje)
  
  # Revisamos que sub_bloque tenga suficientes observaciones para clusterizar
  if (nrow(sub_bloque) < 50) {
    next # Evitemos problemas si el grupo es muy pequeño. .
  }
  
  # (c) Seleccionamos las variables de interés
  sub_bloque_items <- sub_bloque[, descriptivo_grupal]
  
  # (d) Hacemos el Cluster analysis (distancia Jaccard y Ward.D2)
  matriz_binaria <- as.matrix(sub_bloque_items)
  distancias <- dist(matriz_binaria, method = "binary") # d= 1 - shared_1's/tatal_1's
  clust_hier <- hclust(distancias, method = "ward.D2") # squared Euclidean distances. To call hclust()
  
  # (e) Cortamos en k clusters.
  # Aquí usamos un cluster jerarquico. Hay que darle una vuelta a si usamos algún 
  # mecanimos diferente para la generación de los clusters que emerga de los proios 
  # datos y no impuestos por nosotros como en este caso ¿por qué 2, 3, 4 clusters?
  sub_bloque$sub_cluster <- cutree(clust_hier, k = 3) # cambiar el número de clusters en pruebas de robustez
  
  # Revisamos los tamaños de cada cluster
  n_c1 <- sum(sub_bloque$sub_cluster == 1)
  n_c2 <- sum(sub_bloque$sub_cluster == 2)
  n_c3 <- sum(sub_bloque$sub_cluster == 3)
  
  # (f) Estimamos redes de Ising por sub_cluster (si hay suficiente n en cada cluster)
  # Nota para después en caso que haya errores para ventanas muy pequeñas: 
  # Para evitar errores, podríamos considerar después exigir un mínimo de observaciones 
  # en cada cluster. Por ejemplo agregar algo como: if (n_c1 < 10 | n_c2 < 10 | n_c3 < 10) next. 
  
  # Corremos los modelos Ising
  
  model_sub_1 <- estimateNetwork(
    data    = sub_bloque_items[sub_bloque$sub_cluster == 1, ],
    default = "IsingFit"
  )
  model_sub_2 <- estimateNetwork(
    data    = sub_bloque_items[sub_bloque$sub_cluster == 2, ],
    default = "IsingFit"
  )
  model_sub_3 <- estimateNetwork(
    data    = sub_bloque_items[sub_bloque$sub_cluster == 3, ],
    default = "IsingFit"
  )
  
  edges_1 <- model_sub_1$graph
  edges_2 <- model_sub_2$graph
  edges_3 <- model_sub_3$graph
  
  correlacion1_2 <- cor(as.vector(edges_1), as.vector(edges_2))
  correlacion1_3 <- cor(as.vector(edges_1), as.vector(edges_3))
  correlacion2_3 <- cor(as.vector(edges_2), as.vector(edges_3))
  
  # Evaluate cluster quality (silhouette width: 
  # close 1: obsv. well clusteres
  # around 0: obsv. lies between two clusters
  # negative: placed in the wrong cluster)
  silhouette_score <- mean(cluster::silhouette(
    sub_bloque$sub_cluster, distancias)[, "sil_width"])
  
  # Calculate composite score (lower is better)
  composite_score <- (1 - silhouette_score) + mean(c(correlacion1_2,correlacion1_3,correlacion2_3))
  
  # (g) Guardamos las correlaciones y tamaño de clusters en cor_results_sliding_2
  cor_results_sliding_2 <- rbind(
    cor_results_sliding_2,
    data.frame(
      min_puntaje   = start_val,
      max_puntaje   = end_val,
      cor_1_2       = round(correlacion1_2, 3),
      cor_1_3       = round(correlacion1_3, 3),
      cor_2_3       = round(correlacion2_3, 3),
      n_subcluster1 = n_c1,
      n_subcluster2 = n_c2,
      n_subcluster3 = n_c3
    )
  )
  
  # (h) Tablas de contingencia y test de chi-cuadrado para variables_externas
  for (var_externa in variables_externas) {
    
    # Construimos tabla de contingencia
    tabla <- table(sub_bloque[[var_externa]], sub_bloque$sub_cluster)
    
    # Almacenamos tabla en la lista
    nombre_tabla <- paste(var_externa, "_", start_val, "_", end_val, sep = "")
    cross_tables_sliding_2[[nombre_tabla]] <- tabla
    
    # Ejecutamos test de chi-cuadrado 
    test_result <- suppressWarnings(chisq.test(tabla))
    
    chi_results_sliding_2 <- rbind(
      chi_results_sliding_2,
      data.frame(
        min_puntaje = start_val,
        max_puntaje = end_val,
        variable    = var_externa,
        X2          = round(test_result$statistic, 2),
        df          = test_result$parameter,
        p_value     = test_result$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}
time_fin <- Sys.time()
time_total <- time_fin - time_init

write.csv(cor_results_sliding_2, "psy_net_files/cor_results_sliding_2.csv", row.names = FALSE)
write.csv(chi_results_sliding_2, "psy_net_files/chi_results_sliding_2.csv", row.names = FALSE)

cor_results_sliding_2 <- read.csv("psy_net_files/cor_results_sliding_2.csv")
chi_results_sliding_2 <- read.csv("psy_net_files/chi_results_sliding_2.csv")

# 4. Al terminar el bucle, revisamos los resultados ---

cor_results_sliding_2

chi_results_sliding_2

names(cross_tables_sliding_2)  # para ver qué tablas hay

cor_long_2 <- cor_results_sliding_2 %>%
  select(min_puntaje, max_puntaje, cor_1_2, cor_1_3, cor_2_3) %>%
  pivot_longer(
    cols = c("cor_1_2", "cor_1_3", "cor_2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )

# EL MISMO PLOT QUE EN PPT
png("psy_net_plots/cor_evolucion_2.png", width = 800*0.9, height = 600*0.9)
ggplot(cor_long_2, aes(x = min_puntaje, y = correlacion)) +
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
dev.off()

# Convertimos en factores
cor_long_2$clusters_comparados_factor <- factor(
  cor_long_2$clusters_comparados, 
  levels = c("cor_1_2", "cor_1_3", "cor_2_3")
)

ggplot(cor_long_2, aes(x = min_puntaje, y = clusters_comparados_factor, fill = correlacion)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "red", high = "blue", mid = "white",
    midpoint = 0.5, limit = c(0,1), space = "Lab",
    name = "Correlación"
  ) +
  theme_minimal() +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Par de subclusters",
    title = "Heatmap de correlaciones de redes entre subclusters"
  )








# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# TO_DO: Implementar el análisis de las redes con los parámetros de 02 - psy_net_optimal_k_2







# ------------------------------------------------------------------------------
# Análisis psicométricos diferenciados por Rangos de Riesgo  
# ------------------------------------------------------------------------------

base_ising_2 <- base_ising_2 %>%
  dplyr::mutate(risk_category = cut(puntaje_total, 
                                    breaks = c(-Inf, 19, 29, Inf),  # Nuevos breaks
                                    labels = c("Low", "Intermediate", "High"),  # Nombres de los grupos
                                    include.lowest = TRUE))

red_bajo_2 <- estimateNetwork(
  base_ising_2 %>% 
    filter(risk_category == "Low") %>% 
    select(-risk_category),  # Eliminar la columna de clasificación
  default = "IsingFit"
)
red_intermedio_2 <- estimateNetwork(
  base_ising_2 %>% 
    filter(risk_category == "Intermediate") %>% 
    select(-risk_category),  # Eliminar la columna de clasificación
  default = "IsingFit"
)
red_alto_2 <- estimateNetwork(
  base_ising_2 %>% 
    filter(risk_category == "High") %>% 
    select(-risk_category),  # Eliminar la columna de clasificación
  default = "IsingFit"
)

matriz_pesos_bajo_2 <- red_bajo_2$graph
matriz_pesos_medio_2 <- red_intermedio_2$graph
matriz_pesos_alto_2 <- red_alto_2$graph

png("psy_net_recidivism_plots/red_criteriofijo_bajo_2.png", width = 1000, height = 1000)
#red_criteriofijo_bajo <- qgraph(matriz_pesos_bajo, layout = "spring", title = "Riesgo Bajo")
qgraph(matriz_pesos_bajo_2, layout = layout_eganet, title = "Riesgo Bajo")
dev.off()
#layout_eganet <- red_criteriofijo_bajo$layout

png("psy_net_recidivism_plots/red_criteriofijo_medio_2.png", width = 1000, height = 1000)
qgraph(matriz_pesos_medio_2, layout = layout_eganet, title = "Riesgo Medio")
dev.off()

png("psy_net_recidivism_plots/red_criteriofijo_alto_2.png", width = 1000, height = 1000)
qgraph(matriz_pesos_alto_2, layout = layout_eganet, title = "Riesgo Alto")
dev.off()

# Veamos ahora si se detectan comunidades diferentes en estas tres redes obtenidas.

# EGA para riesgo bajo
datos_ajustados_bajo <- base_ising %>%
  filter(risk_category == "Low") %>%
  select(-risk_category) %>%
  mutate(across(everything(), ~ . + 1))  # Mapea de -1:1 a 0:2

ega_bajo <- EGA(datos_ajustados_bajo)

# EGA para riesgo intermedio
datos_ajustados_bajo <- base_ising %>%
  filter(risk_category == "Intermediate") %>%
  select(-risk_category) %>%
  mutate(across(everything(), ~ . + 1))  # Mapea de -1:1 a 0:2

#ega_intermedio <- EGA(datos_ajustados_bajo)
EGA(datos_ajustados_bajo)

# EGA para riesgo alto
datos_ajustados_alto <- base_ising %>%
  filter(risk_category == "High") %>%
  select(-risk_category) %>%
  mutate(across(everything(), ~ . + 1))  # Mapea de -1:1 a 0:2

#ega_intermedio <- EGA(datos_ajustados_alto)
EGA(datos_ajustados_alto)

# Nota: Las visualizaciones sugieren que hay espacio para una exploración más detallada.
