library(bootnet)
library(qgraph)
library(igraph)

# ------------------------------------------------------------------------------
# Modelo de Ising - correlación condicionando el resto de variables
# ------------------------------------------------------------------------------

# Cargamos datos procesados

base_igi <- read.csv("psy_net_recidivism_files/base_igi_ising.csv")
base_ptje <- read.csv("psy_net_recidivism_files/base_igi_ising.csv")

# NOTA: Aquí hicimos un Ising apegado al estándar de valores 1 y -1. 
# No obstante, en etapas posteriores esto generaba errores para algunos cálculos que 
# requerían sólo valores positivos, por lo que finalmente lo dejamos con valores 0s y 1s.
# Dejamos igual el código para generar una versión con 1s y -1s por si eventualmente cambianos nuevamente. 

# Convertimos los valores 0 en valore -1 en todas las variables. 
convert_zeros_to_neg1 <- function(df, var_names) { # Iteramos por cada nombre de variable 
  for (col in var_names) { # Verificamos que la columna exista en el data frame 
    if (col %in% names(df)) { # Reemplazamos los valores 0 por -1 
      df[[col]][df[[col]] == 0] <- 0 } 
    else { warning(paste("La columna", col, "no existe en el data frame.")) } 
  } 
  return(df) 
}

# Cambiamos nombre de variables
descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR24","PAR25","PAR26","PAR27",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

# Seleccionamos variables de interés (LO HARÉ CON BASE_PTJE)
base_igi_ising_2 <-  base_ptje[,colnames(base_ptje) %in% descriptivo_grupal]
base_ising_2 <- convert_zeros_to_neg1(base_igi_ising_2, descriptivo_grupal)

# TO_DO: implementar con IsingFit directamente.
red_dicotomica_2 <- bootnet::estimateNetwork(base_ising_2, 
                                             default = "IsingFit",
                                             #principalDirection =T,
                                             tuning = 0.25,
                                             labels = descriptivo_grupal)

png("psy_net_recidivism_plots/red_base_igi.png", width = 1000, height = 1000)
qgraph(red_dicotomica_2$graph, layout = "spring", labels = colnames(red_dicotomica_2$graph))
dev.off()

# ------------------------------------------------------------------------------
# Estabilidad - ver psy_net_stability.R
# ------------------------------------------------------------------------------

# \\~\\

# ------------------------------------------------------------------------------
# Análisis psicométricos por Rangos de Riesgo  (Analizar después con k adecuado)
# ------------------------------------------------------------------------------

# \\~\\

# ------------------------------------------------------------------------------
# Bootnet - Ventana de 10 puntos  
# ------------------------------------------------------------------------------

# vamos al centro de nuestro objetivo: 
# veamos si para riesgos similares se producen diferencias sustantivas en las redes generadas. 
# Para ello definiremos ventanas de riesgos de 5, 10, 15 puntos y los análisis los haremos dentro de cada ventana.

# Las ventanas serán móviles. Es decir, si la ventana es de 10 puntos, se calcula para
# ventana que parte en 0 hasta 9, luego de 1 hasta 10, y así sucesivamente hasta 34 a 43.

# 1. Definimos contenedores para resultados ---
cor_results_sliding <- data.frame()
chi_results_sliding <- data.frame()
cross_tables_sliding <- list()

# Definimos las variables externas a explorar después en las regresiones
variables_externas <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# 2. Definimos parámetros de la ventana ---
window_size <- 10 # Modificar la ventana para otros ejercicios
max_puntaje_observado <- max(base_ptje$puntaje_total)  #-->> usamos ¿¿ base_igi ?? base_ptje
inicios_ventana <- 0:(max_puntaje_observado - window_size + 1)

# 3. Bucle principal por ventanas "traslapadas" ---
for (start_val in inicios_ventana) {
  
  # (a) Determinamos el rango de la ventana
  end_val <- start_val + window_size - 1  #  --->>> ¿por qué -5 y no -1? CAMBIÉ A -1
  rango_puntaje <- start_val:end_val
  
  # (b) Filtramos la base por rango_puntaje  -->> usamos ¿¿ base_igi ?? base_ptje
  sub_bloque <- base_ptje %>%
    dplyr::filter(puntaje_total %in% rango_puntaje)
  
  # Revisamos que sub_bloque tenga suficientes observaciones para clusterizar
  if (nrow(sub_bloque) < 50) {
    next # Evitemos problemas si el grupo es muy pequeño. .
  }
  
  # (c) Seleccionamos las variables de interés
  sub_bloque_items <- sub_bloque[, descriptivo_grupal]
  
  # (d)
  
  
  
  
  
  
  
  
  
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
  
  # (g) Guardamos las correlaciones y tamaño de clusters en cor_results_sliding
  cor_results_sliding <- rbind(
    cor_results_sliding,
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
    cross_tables_sliding[[nombre_tabla]] <- tabla
    
    # Ejecutamos test de chi-cuadrado 
    test_result <- suppressWarnings(chisq.test(tabla))
    
    chi_results_sliding <- rbind(
      chi_results_sliding,
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

# --- Revisión objeto clust_hier (rango puntaje 32 33 34 35 36 37 - 'el último clust_hier generado')

# Inspecting the Dendrogram
png("psy_net_recidivism_plots/method_dendrogram.png", width = 800, height = 600)
plot(clust_hier)
#heights <- sort(clust_hier$height, decreasing = TRUE)
#largest_gap <- max(diff(heights))
#threshold <- heights[which(diff(heights) == largest_gap)]
abline(h = 1.75, col = "red")
dev.off() 

# Elbow Method (Within-Cluster Sum of Squares)
png("psy_net_recidivism_plots/method_elbow.png", width = 800, height = 600)
wss <- sapply(1:10, function(k) { 
  sum(cutree(clust_hier, k = k)^2)
})
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "WSS") # “elbow” in the curve
dev.off()

# Silhouette Method
library(cluster)
png("psy_net_recidivism_plots/method_silhouette.png", width = 800, height = 600)
silhouette_scores <- sapply(2:10, function(k) {
  mean(silhouette(cutree(clust_hier, k = k), distancias)[, 3])
})
plot(2:10, silhouette_scores, type = "b", xlab = "Number of Clusters", ylab = "Average Silhouette Width")
optimal_k <- which.max(silhouette_scores)
dev.off()

# Gap Statistic (To_Do : Lento - optimizar !!!! )
library(cluster)
png("psy_net_recidivism_plots/method_gap.png", width = 800, height = 600)
make_cluster <- function(x, k) {
  list(cluster = cutree(hclust(dist(x), method = "ward.D2"), k = k))
}
gap_stat <- clusGap(
  as.matrix(sub_bloque_items), 
  FUN = make_cluster, 
  K.max = 10, 
  B = 100, # number of Monte Carlo (“bootstrap”) samples
  spaceH0 = "scaledPCA" # second choice from Tibshirani 2001, pag 414.
)
plot(gap_stat) # (the “1-SE rule”) optimal is chosen as the smallest k  where `Gap_k`
grid()         # is within one standard error of the maximum value. k=5 in this case.
# The optimal number of clusters is the smallest k such that: 
# Gap(k) ≥ Gap(k+1) - s_{k+1}. 
# Where s_{k+1} is the standard error of Gap(k+1)
dev.off()


# --- Fin revisión


# 4. Al terminar el bucle, revisamos los resultados ---

# Correlaciones entre matrices, con el rango de puntajes y n de cada subcluster
cor_results_sliding

# Resultados de chi cuadrado
chi_results_sliding

# Tablas de contingencia almacenadas
names(cross_tables_sliding)  # para ver qué tablas hay

#library(tidyverse)
library(tidyr)

cor_long <- cor_results_sliding %>%
  select(min_puntaje, max_puntaje, cor_1_2, cor_1_3, cor_2_3) %>%
  pivot_longer(
    cols = c("cor_1_2", "cor_1_3", "cor_2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )

cor_long_b <- cor_long[1:96,]

# EL MISMO PLOT QUE EN PPT
ggplot(cor_long_b, aes(x = min_puntaje, y = correlacion)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ clusters_comparados, ncol = 1) +
  theme_minimal() +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Correlación de la red",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  )

# Convertimos en factores
cor_long$clusters_comparados_factor <- factor(
  cor_long$clusters_comparados, 
  levels = c("cor_1_2", "cor_1_3", "cor_2_3")
)
cor_long_b$clusters_comparados_factor <- factor(
  cor_long_b$clusters_comparados, 
  levels = c("cor_1_2", "cor_1_3", "cor_2_3")
)

ggplot(cor_long, aes(x = min_puntaje, y = clusters_comparados_factor, fill = correlacion)) +
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
ggplot(cor_long_b, aes(x = min_puntaje, y = clusters_comparados_factor, fill = correlacion)) +
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

library(EGAnet)

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
