# 1) Usamos codificación {3,2,1},{0} --> {0},{1}
# 2) Quitar PAR24 y HD6 --> muy poca varianza
# 3) Nos quedamos con los vantanas de puntaje inicial 10 - 20
# 4) Limpiaremos los casos para aislar solo aquellos que corresponden a delitos de
#   4.1 Robo
#   4.2 Droga
#   4.3 Delitos sexuales
#   4.4 Delitos económicos
#   4.5 Violencia intrafamiliar
# 5) Seleccionaremos a los reos con casos puros: sin imputaciones de más de 1 delito.

library(dplyr)
library(tidyr)
library(ggplot2)
library(bootnet)
library(qgraph)
library(cluster) # For silhouette calculation
library(viridis) # For color scales
library(igraph)

# ------------------------------------------------------------------------------
# Creación objeto para 03 - intermediate_values_split.R
# ------------------------------------------------------------------------------

base_igi_bin_1 <- read.csv("psy_net_files/base_igi_bin_1.csv")
reinc_por_reo <- readRDS("documentos carceles - ising/Datos reos - sensible/reinc_por_reo.rds")
datos_reos_reinc <- readRDS("documentos carceles - ising/Datos reos - sensible/datos_reos_reinc.rds")
delitos_unicos_clasificados <- readRDS("psy_net_files/delitos_unicos_clasificados.rds") 

length(unique(base_igi_bin_1$COD_PERS)) # 17276
length(unique(reinc_por_reo$COD_PERS)) # 97123
length(datos_reos_reinc$COD_PERS) # 203130 --> hay reos que han reincidido varias veces 

# Como hay reos que han reincidido varias veces, hay que tomar una decisión: 
#   1) tratamos cada reincidencia como un caso nuevo ('nuevo reo'), o
#   2) tratamos cada reo como único, y consideramos solo su primera reincidencia
# Acá adoptamos la opción 1).

# Seleccionamos "COD_PERS" , "SEXO", "COMUNA_DOMICILIO", "FECHA_NACIMIENTO", "DELITOS"
reos_filtrados <- datos_reos_reinc %>%
  select(COD_PERS, SEXO, COMUNA_DOMICILIO, FECHA_NACIMIENTO, EGRESO, DELITOS) #203130

# Filtramos las columnas que tienen delitos compuestos. Solo un delito (#171755) 
reos_filtrados <- reos_filtrados %>%
  mutate(
    num_delitos = sapply(strsplit(DELITOS, "; "), length), # Columna temporal
    DELITOS = gsub(" \\(\\d+\\)", "", DELITOS)
  ) %>%
  filter(num_delitos == 1) %>%
  select(-num_delitos) # Eliminamos columna temporal

# Porcentaje de duplicados
(1 - n_distinct(reos_filtrados$COD_PERS)/nrow(reos_filtrados)) * 100

# -------------------- Creamos base_igi_bin_1_filtrado (solo códigos únicos) -----

# all(unique(base_igi_bin_1$COD_PERS) %in% unique(reos_filtrados$COD_PERS))
# 
# # 7006 códigos están en ambas bases de datos, un 40.553 % de las disponibles en base_igi
# # (para el resto de casos, con num_delitos > 1, hay otro 23.385 %)
# cod_pers_baseigi_inter_reinc <- intersect(unique(base_igi_bin_1$COD_PERS), unique(reos_filtrados$COD_PERS))
# (length(cod_pers_baseigi_inter_reinc) / length(unique(base_igi_bin_1$COD_PERS))) * 100

base_igi_bin_1_filtrado <- base_igi_bin_1 %>%
  filter(COD_PERS %in% reos_filtrados$COD_PERS)

# -------------------- ... -----------------------------------------------------
# -------------------- ... -----------------------------------------------------

# base_igi_bin_1 SOLO TIENE CÓDIGOS ÚNICOS, ASÍ QUE NO LO PUEDO USAR, PORQUE HAY AMBIGUEDAD 
# DE CUÁL ES LA OCASIÓN DE INGRESO DE CADA REO. DEBO VOLVER A 01 - psy_net_recidivism.R y arreglarlo.
#
# ACTUALIZACIÓN: base_igi SOLO CONTIENE CÓDIGOS ÚNICOS. HAY QUE SELECCIONAR AQUELLA FILA 
# QUE CONTIENE LA VENTANA TEMPORAL CUANDO SE HIZO EL IGI.
# length(unique(base_igi_bin_1$COD_PERS)) == length(base_igi_bin_1$COD_PERS) # TRUE
# length(unique(base_igi$COD_PERS)) == length(base_igi$COD_PERS) # TRUE

# OTRA OPCIÓN ES VER CUÁNTO VARÍA LA CATEGORÍA DEL DELITO DE AQUELLOS QUE TIENEN 1 ACUSACIÓN

# -------------------- ... -----------------------------------------------------

# # Unir categorías al dataframe principal
# reos_filtrados <- reos_filtrados %>% 
#   left_join(delitos_unicos_clasificados %>% select(DELITOS, CATEGORIA_DELITO), 
#             by = "DELITOS") %>% 
#   mutate(CATEGORIA_DELITO = ifelse(is.na(CATEGORIA_DELITO), 
#                                    "Falla categoría", 
#                                    CATEGORIA_DELITO))
# 
# # 0 en 'Falla categoría' !
# nrow(reos_filtrados %>% filter(CATEGORIA_DELITO == "Falla categoría")) 
# 
# # 70.6 % de los reos con solo 1 delito han reingresado más de 1 vez !
# reos_filtrados_multiReinc <- reos_filtrados %>%
#   group_by(COD_PERS) %>%
#   filter(n() > 1) %>%  # Selecciona códigos con más de 1 registro
#   ungroup()
# 
# analisis_variacion <- reos_filtrados_multiReinc %>%
#   group_by(COD_PERS) %>%
#   summarise(
#     Total_Delitos = n(),
#     Categorias_Unicas = n_distinct(CATEGORIA_DELITO),
#     Lista_Categorias = toString(unique(CATEGORIA_DELITO))
#   ) %>%
#   arrange(desc(Categorias_Unicas))
# 
# # Distribución de categorías únicas
# distribucion <- analisis_variacion %>%
#   count(Categorias_Unicas, name = "N_Reos") %>%
#   mutate(Porcentaje = round(N_Reos/nrow(analisis_variacion)*100, 2))
# 
# distribucion
# output: Parte importante tiene 2 categorías ---> cambian de categoría.
#  Categorias_Unicas | N_Reos | Porcentaje
# 1                  | ~17k  | 45.3%
# 2                  | ~15k  | 40.3%
# 3                  | ~4k   | 11.1%
# 4                  | ~1k   | 2.6%
# 5                  | 211   | 0.5%

# -------------------- ... -----------------------------------------------------

# SÍ, VARÍA HARTO, ASÍ QUE NO PUEDO OCUPAR CUALQUIERA, DEBO USAR EL QUE CORRESPONDE
# AL MOMENTO EN QUE SE HIZO LA ENTREVISTRA

# min(base_igi_bin_1_filtrado$FECHA_ENTREVISTA)
# max(base_igi_bin_1_filtrado$FECHA_ENTREVISTA)

# Veamos el porcentaje de entrevistas realizadas antes del egreso de los reos

# # 1. Unir ambas bases por COD_PERS
# datos_comparacion <- base_igi_bin_1_filtrado %>% 
#   inner_join(datos_reos_reinc %>% select(COD_PERS, EGRESO), by = "COD_PERS")
# 
# # 2. Convertir fechas a formato Date
# datos_comparacion <- datos_comparacion %>%
#   mutate(
#     FECHA_ENTREVISTA = as.Date(FECHA_ENTREVISTA),
#     EGRESO = as.Date(EGRESO)
#   )
# 
# # 3. Calcular condición y porcentaje 
# datos_comparacion %>%
#   summarise(
#     total_casos = n(),
#     entrevistas_antes = sum(FECHA_ENTREVISTA < EGRESO, na.rm = TRUE),
#     porcentaje = round(entrevistas_antes / total_casos * 100, 2)
#   )
# 
# # 38.67 % de las entrevistas fueron antes de EGRESO. Es decir, podemos ocupar datos de 
# length(datos_reos_reinc$COD_PERS) * 0.3867
# # 78550 reos.

# -------------------- ... -----------------------------------------------------
# -------------------- ... -----------------------------------------------------

# El problema es que un mismo código COD_PERS está una única vez en base_igi_bin_1_filtrado, 
# pero puede estar varias veces en reos_filtrados. Por ahora, cada vez que no exista un 
# único valor de código para hacer un match, se tendrá que elegir aquel que tenga el valor 
# reos_filtrados$EGRESO más posterior, es decir, el egreso que está más en el futuro.

# Unir categorías al dataframe principal
reos_filtrados <- reos_filtrados %>% 
  left_join(delitos_unicos_clasificados %>% select(DELITOS, CATEGORIA_DELITO), 
            by = "DELITOS") %>% 
  mutate(CATEGORIA_DELITO = ifelse(is.na(CATEGORIA_DELITO), 
                                   "Falla categoría", 
                                   CATEGORIA_DELITO))

# 0 en 'Falla categoría' !
nrow(reos_filtrados %>% filter(CATEGORIA_DELITO == "Falla categoría")) 

# Paso 1: Crear dataframe con el último egreso por COD_PERS (88350 obs)
reos_ultimo_egreso <- reos_filtrados %>%
  mutate(EGRESO = as.Date(EGRESO)) %>%
  group_by(COD_PERS) %>%
  arrange(desc(EGRESO)) %>%
  slice(1) %>%  # Selecciona el registro más reciente por COD_PERS
  ungroup() %>%
  select(COD_PERS, DELITOS, CATEGORIA_DELITO)

# Paso 2: Unir con base_igi_bin_1_filtrado
base_igi_bin_1_filtrado <- base_igi_bin_1_filtrado %>%
  left_join(reos_ultimo_egreso, by = "COD_PERS")

# Nombres columnas
colnames(base_igi_bin_1_filtrado)

# ------------------------------------------------------------------------------
# Bootnet - Ventana de 10 puntos  
# ------------------------------------------------------------------------------

descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD7","HD8",#"HD6",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR25","PAR26","PAR27", #"PAR24",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

# Definimos las variables externas a explorar después en las regresiones
variables_externas <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

#### HACER ANÁLISIS SEPARANDO POR TIPOS DE DELITO:
# 4) Limpiaremos los casos para aislar solo aquellos que corresponden a delitos de
#   4.1 Robo
#   4.2 Droga
#   4.3 Delitos sexuales
#   4.4 Delitos económicos
#   4.5 Violencia intrafamiliar

delitos_lista <- unique(base_igi_bin_1_filtrado$CATEGORIA_DELITO)

obs_por_delito <- as.data.frame(table(base_igi_bin_1_filtrado$CATEGORIA_DELITO))
colnames(obs_por_delito) <- c("Categoría", "Casos")
obs_por_delito[order(-obs_por_delito$Casos), ]

# Función principal ------------------------------------------------------------
barrido_param_igi_puntaje <- function(
    base_igi_bin,
    categoria_delito = NULL,
    descriptivo_grupal,
    variables_externas,
    window_size,
    window_starts,
    k_clusters = 3,
    min_window_obs = 50,
    min_cluster_obs = 10,
    distance_method = "binary",
    hclust_method = "ward.D2") {
  
  # Inicializar objetos para resultados
  resultados_cor <- data.frame()
  resultados_chi2 <- data.frame()
  tablas_contingencia <- list()
  
  # 1. Filtrar por categoría delito si se especifica
  if (!is.null(categoria_delito)) {
    if (!"CATEGORIA_DELITO" %in% colnames(base_igi_bin)) {
      stop("Columna 'CATEGORIA_DELITO' no encontrada en el dataframe")
    }
    if (!categoria_delito %in% unique(base_igi_bin$CATEGORIA_DELITO)) {
      stop(paste("Categoría inválida. Opciones válidas:\n", 
                 paste(unique(base_igi_bin$CATEGORIA_DELITO), collapse = "\n")))
    }
    base_igi_bin <- base_igi_bin %>% 
      dplyr::filter(CATEGORIA_DELITO == categoria_delito)
  }
  
  # 2. Validar existencia de columnas requeridas
  columnas_requeridas <- c("puntaje_total", descriptivo_grupal, variables_externas)
  if (!all(columnas_requeridas %in% colnames(base_igi_bin))) {
    faltantes <- setdiff(columnas_requeridas, colnames(base_igi_bin))
    stop(paste("Columnas faltantes:", paste(faltantes, collapse = ", ")))
  }
  
  # 3. Procesar cada ventana
  for (inicio in window_starts) {
    fin <- inicio + window_size - 1
    sub_bloque <- base_igi_bin %>%
      dplyr::filter(puntaje_total >= inicio, puntaje_total <= fin)
    
    # Saltar ventanas con pocos datos
    if (nrow(sub_bloque) < min_window_obs) next
    
    # 4. Clusterización
    matriz_items <- as.matrix(sub_bloque[, descriptivo_grupal])
    distancias <- dist(matriz_items, method = distance_method)
    clust <- hclust(distancias, method = hclust_method)
    sub_bloque$sub_cluster <- cutree(clust, k = k_clusters)
    
    # 5. Validar tamaño de clusters
    conteos_cluster <- table(sub_bloque$sub_cluster)
    if (any(conteos_cluster < min_cluster_obs)) next
    
    # 6. Modelos Ising por cluster
    modelos_ising <- lapply(1:k_clusters, function(k) {
      datos_cluster <- sub_bloque[sub_bloque$sub_cluster == k, descriptivo_grupal]
      if (nrow(datos_cluster) == 0) return(matrix(NA))
      estimateNetwork(datos_cluster, default = "IsingFit")$graph
    })
    
    # 7. Calcular correlaciones entre modelos (versión corregida)
    correlacion1_2 <- cor(
      as.vector(modelos_ising[[1]]),
      as.vector(modelos_ising[[2]]),
      use = "pairwise.complete.obs",
      method = "pearson"  # Especificar método explícitamente
    )
    
    correlacion1_3 <- cor(
      as.vector(modelos_ising[[1]]),
      as.vector(modelos_ising[[3]]),
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    
    correlacion2_3 <- cor(
      as.vector(modelos_ising[[2]]),
      as.vector(modelos_ising[[3]]),
      use = "pairwise.complete.obs",
      method = "pearson"
    )
    
    # 8. Métricas de calidad
    silhouette <- mean(cluster::silhouette(sub_bloque$sub_cluster, distancias)[, "sil_width"])
    puntaje_compuesto <- (1 - silhouette) + mean(c(correlacion1_2,correlacion1_3,correlacion2_3))
    
    # 9. Almacenar resultados correlacionales
    resultados_cor <- rbind(resultados_cor, data.frame(
      inicio_ventana = inicio,
      fin_ventana = fin,
      cor_cluster1_2 = round(correlacion1_2, 3),
      cor_cluster1_3 = round(correlacion1_3, 3),
      cor_cluster2_3 = round(correlacion2_3, 3),
      n_cluster1 = conteos_cluster[1],
      n_cluster2 = conteos_cluster[2],
      n_cluster3 = conteos_cluster[3],
      silhouette = round(silhouette, 3),
      puntaje_compuesto = round(puntaje_compuesto, 3)
    ))
    
    # 10. Pruebas chi-cuadrado para variables externas
    for (var in variables_externas) {
      tabla <- table(sub_bloque[[var]], sub_bloque$sub_cluster)
      test_chi <- suppressWarnings(chisq.test(tabla))
      
      tablas_contingencia[[paste(var, inicio, fin, sep = "_")]] <- tabla
      
      resultados_chi2 <- rbind(resultados_chi2, data.frame(
        inicio_ventana = inicio,
        fin_ventana = fin,
        variable = var,
        estadistico_chi2 = round(test_chi$statistic, 2),
        valor_p = ifelse(
          test_chi$p.value < 0.001, 
          "<0.001", 
          format(round(test_chi$p.value, 3), nsmall = 3)
        )
      ))
    }
  }
  
  # 11. Retornar resultados estructurados
  list(
    resultados_correlaciones = resultados_cor,
    resultados_chi2 = resultados_chi2,
    tablas_contingencia = tablas_contingencia
  )
}

# ------------------------------------------------------------------------------

# Bucle por característica de delito ---

print(delitos_lista[7])

time_init <- Sys.time()
barrido_param_robo <- barrido_param_igi_puntaje(
  base_igi_bin = base_igi_bin_1_filtrado,
  categoria_delito = delitos_lista[7], # "Robo"
  descriptivo_grupal = descriptivo_grupal, 
  variables_externas = variables_externas,
  window_size = 10,
  window_starts = seq(10, 20, by = 1)
)
time_fin <- Sys.time()
time_total <- time_fin - time_init  # 48.11 sec

str(barrido_param_robo)

barrido_param_robo$resultados_correlaciones

saveRDS(barrido_param_robo, "psy_net_files/barrido_param_robo.rds")
barrido_param_robo <- readRDS("psy_net_files/barrido_param_robo.rds")

# Transformar datos a formato largo
cor_long_robo <- barrido_param_robo$resultados_correlaciones %>%
  select(inicio_ventana, fin_ventana, cor_cluster1_2, cor_cluster1_3, cor_cluster2_3) %>%
  pivot_longer(
    cols = c("cor_cluster1_2", "cor_cluster1_3", "cor_cluster2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )

# Crear gráfico facetado
png("psy_net_plots/barrido_param_robo.png", width = 800*0.9, height = 600*0.9)
ggplot(cor_long_robo, aes(x = inicio_ventana, y = correlacion)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ clusters_comparados, ncol = 1) +
  theme_minimal() + 
  ylim(0, 1) +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Correlación de la red",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  ) +
  theme(
    panel.grid.minor = element_line(color = "gray93"),
    panel.grid.major = element_line(color = "gray93"),
    strip.text = element_text(face = "bold")
  )
dev.off()

# ------------------------------------------------------------------------------
# Comparación redes con distintos delitos
# ------------------------------------------------------------------------------

# Función principal ------------------------------------------------------------
barrido_param_igi_puntaje_2 <- function(
    base_igi_bin,
    categoria_delito = NULL,
    descriptivo_grupal,
    variables_externas,
    window_size,
    window_starts,
    k_clusters = 3,
    min_window_obs = 50,
    distance_method = "binary",
    hclust_method = "ward.D2") {
  
  # Inicializar objetos para resultados
  resultados_cor <- data.frame()
  resultados_chi2 <- data.frame()
  
  # 1. Filtrar por categoría delito si se especifica
  if (!is.null(categoria_delito)) {
    if (!"CATEGORIA_DELITO" %in% colnames(base_igi_bin)) {
      stop("Columna 'CATEGORIA_DELITO' no encontrada en el dataframe")
    }
    base_igi_bin1 <- base_igi_bin %>% 
      dplyr::filter(CATEGORIA_DELITO == categoria_delito[1])
    base_igi_bin2 <- base_igi_bin %>% 
      dplyr::filter(CATEGORIA_DELITO == categoria_delito[2])
    base_igi_bin3 <- base_igi_bin %>% 
      dplyr::filter(CATEGORIA_DELITO == categoria_delito[3])
  }
  
  # 2. Validar existencia de columnas requeridas
  
  # 3. Procesar cada ventana
  for (inicio in window_starts) {
    fin <- inicio + window_size - 1
    
    sub_bloque1 <- base_igi_bin1 %>%
      dplyr::filter(puntaje_total >= inicio, puntaje_total <= fin)
    sub_bloque1 <- sub_bloque1[, descriptivo_grupal]
    sub_bloque2 <- base_igi_bin2 %>%
      dplyr::filter(puntaje_total >= inicio, puntaje_total <= fin)
    sub_bloque2 <- sub_bloque2[, descriptivo_grupal]
    sub_bloque3 <- base_igi_bin3 %>%
      dplyr::filter(puntaje_total >= inicio, puntaje_total <= fin)
    sub_bloque3 <- sub_bloque3[, descriptivo_grupal]
    
    # Saltar ventanas con pocos datos
    if (nrow(sub_bloque1) < min_window_obs || nrow(sub_bloque2) < min_window_obs
        || nrow(sub_bloque3) < min_window_obs) {
      message('ventana con muy pocos datos')
      next
      } 
    
    # 4. Clusterización
    # matriz_items <- as.matrix(sub_bloque[, descriptivo_grupal])
    # distancias <- dist(matriz_items, method = distance_method)
    # clust <- hclust(distancias, method = hclust_method)
    # sub_bloque$sub_cluster <- cutree(clust, k = k_clusters)
    
    # 5. Validar tamaño de clusters
    
    # 6. Modelos Ising por cluster
    model_sub_1 <- estimateNetwork(
      data    = sub_bloque1,
      default = "IsingFit"
    )
    model_sub_2 <- estimateNetwork(
      data    = sub_bloque2,
      default = "IsingFit"
    )
    model_sub_3 <- estimateNetwork(
      data    = sub_bloque3,
      default = "IsingFit"
    )
    
    # 7. Calcular correlaciones entre modelos (versión corregida)
    correlacion1_2 <- cor(as.vector(model_sub_1$graph), as.vector(model_sub_2$graph))
    correlacion1_3 <- cor(as.vector(model_sub_1$graph), as.vector(model_sub_3$graph))
    correlacion2_3 <- cor(as.vector(model_sub_2$graph), as.vector(model_sub_3$graph))
    
    # 8. Métricas de calidad
    # silhouette <- mean(cluster::silhouette(sub_bloque$sub_cluster, distancias)[, "sil_width"])
    # puntaje_compuesto <- (1 - silhouette) + mean(c(correlacion1_2,correlacion1_3,correlacion2_3))
    
    # 9. Almacenar resultados correlacionales
    resultados_cor <- rbind(resultados_cor, data.frame(
      inicio_ventana = inicio,
      fin_ventana = fin,
      cor_cluster1_2 = round(correlacion1_2, 3),
      cor_cluster1_3 = round(correlacion1_3, 3),
      cor_cluster2_3 = round(correlacion2_3, 3),
      n_net1 = nrow(sub_bloque1),
      n_net2 = nrow(sub_bloque2),
      n_net3 = nrow(sub_bloque3)#,
      # silhouette = round(silhouette, 3),
      # puntaje_compuesto = round(puntaje_compuesto, 3)
    ))
    
  }
  
  # 11. Retornar resultados estructurados
  list(
    resultados_correlaciones = resultados_cor,
    resultados_chi2 = resultados_chi2
  )
}

# ------------------------------------------------------------------------------

# Bucle por característica de delito ---

# 4) Limpiaremos los casos para aislar solo aquellos que corresponden a delitos de
#   4.1 Robo
#   4.2 Droga
#   4.3 Delitos sexuales
#   4.4 Delitos económicos
#   4.5 Violencia intrafamiliar

obs_por_delito[order(-obs_por_delito$Casos), ]

print(delitos_lista[7]) # "Robo"
print(delitos_lista[4]) # "Delitos de Drogas"
print(delitos_lista[3]) # "Delitos Sexuales"
categorias_delitos_input <- c(delitos_lista[7], delitos_lista[4], delitos_lista[3])

time_init <- Sys.time()
barrido_param_1 <- barrido_param_igi_puntaje_2(
  base_igi_bin = base_igi_bin_1_filtrado,
  categoria_delito = categorias_delitos_input,
  descriptivo_grupal = descriptivo_grupal, 
  variables_externas = variables_externas,
  window_size = 10,
  window_starts = seq(0, 30, by = 1),
  min_window_obs = 10
)
time_fin <- Sys.time()
time_total <- time_fin - time_init  # 1.05 min

str(barrido_param_1)

barrido_param_1$resultados_correlaciones

saveRDS(barrido_param_1, "psy_net_files/barrido_param_1.rds")
barrido_param_1 <- readRDS("psy_net_files/barrido_param_1.rds")

# Transformar datos a formato largo
cor_long_1 <- barrido_param_1$resultados_correlaciones %>%
  select(inicio_ventana, fin_ventana, cor_cluster1_2, cor_cluster1_3, cor_cluster2_3) %>%
  pivot_longer(
    cols = c("cor_cluster1_2", "cor_cluster1_3", "cor_cluster2_3"),
    names_to = "clusters_comparados",
    values_to = "correlacion"
  )

# Crear gráfico facetado
png("psy_net_plots/barrido_param_Rob-Sex-Drog.png", width = 800*0.9, height = 600*0.9)
ggplot(cor_long_1, aes(x = inicio_ventana, y = correlacion)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ clusters_comparados, ncol = 1,
             labeller = labeller(clusters_comparados = c(
               "cor_cluster1_2" = "Robo y Drogas",
               "cor_cluster1_3" = "Robo y Sexuales",
               "cor_cluster2_3" = "Drogas y Sexuales"
             ))) +
  theme_minimal() + 
  ylim(0, 1) +
  labs(
    x = "Puntaje mínimo de la ventana",
    y = "Correlación de la red",
    title = "Evolución de la correlación entre subclusters según ventana de puntajes"
  ) +
  theme(
    panel.grid.minor = element_line(color = "gray93"),
    panel.grid.major = element_line(color = "gray93"),
    strip.text = element_text(face = "bold")
  )
dev.off()
