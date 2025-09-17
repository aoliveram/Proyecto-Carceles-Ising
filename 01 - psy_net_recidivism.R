# 1) Limpieza y binarización del IGI (dos esquemas: {3,2,1}/{0} y {3,2}/{1,0}).
# 2) Cálculo de puntaje_total IGI y descriptivos/plots básicos por riesgo.
# 3) Estimación de redes Ising globales por esquema (IsingFit); layout y métricas/comunidades.
# 4) Redes estratificadas por riesgo (bajo, intermedio, alto) y comparación.
# 5) Barrido por ventanas de puntaje (ancho ~10) con subclusters (k=3): redes por subcluster,
#    correlación entre grafos, silhouette y puntaje compuesto; chi-cuadrado con variables externas.
#
# Archivos que genera (salida principal):
#  - psy_net_files/base_igi.csv
#  - psy_net_files/base_igi_ising.csv
#  - psy_net_files/base_ptje_ising.csv
#  - psy_net_files/base_igi_bin_1.csv
#  - psy_net_files/base_igi_bin_2.csv
#  - psy_net_files/cor_results_sliding_1.csv
#  - psy_net_files/cor_results_sliding_2.csv
#  - psy_net_files/chi_results_sliding_1.csv
#  - psy_net_files/chi_results_sliding_2.csv
#  - psy_net_files/delitos_unicos_clasificados.rds
#  - Gráficos varios en psy_net_plots/

library(dplyr)
library(tidyr)
library(readxl)
library(EGAnet)
library(stringr)
library(car)
library(psychonetrics)
library(cooccur)
library(visNetwork)
library(igraph)
library(ggplot2)

# ------------------------------------------------------------------------------
# Análisis Exploratorio
# ------------------------------------------------------------------------------

# Importamos
base_igi <- read.csv("psy_net_recidivism3_revisar_files/base_psy_net_recidivism.csv")

length(unique(base_igi$COD_PERS)) == length(base_igi$COD_PERS)

# DECRIPTIVOS --- 

sin_edad <- which(base_igi$EDAD == "S.I.")

# Convertimos la variable EDAD a factor para manejar valores como "S.I."
base_igi$EDAD <- as.factor(base_igi$EDAD)

# Creamos la variable 'rangos_edad' agrupando las edades en las categorías deseadas
base_igi$rangos_edad <- dplyr::case_when(
  as.numeric(as.character(base_igi$EDAD)) < 25 ~ "Menores de 25",
  as.numeric(as.character(base_igi$EDAD)) >= 25 & as.numeric(as.character(base_igi$EDAD)) <= 39 ~ "Entre 25 y 39",
  as.numeric(as.character(base_igi$EDAD)) >= 40 ~ "40 o más",
  TRUE ~ "Sin información"
)

# Convertimos la nueva variable en factor para facilitar el análisis posterior
base_igi$rangos_edad <- factor(base_igi$rangos_edad, 
                               levels = c("Menores de 25", "Entre 25 y 39", "40 o más", "Sin información"))

# Estado civil
base_igi$recod_estado_civil <- dplyr::case_when(
  base_igi$ESTADO_CIVIL %in% c("CASADO", "CONVIVIENTE CIVIL") ~ "Con pareja",
  base_igi$ESTADO_CIVIL %in% c("SEPARADO/DIVORCIADO", "VIUDO", "SOLTERO") ~ "Sin pareja",
  TRUE ~ NA_character_
)

base_igi$recod_estado_civil <- factor(base_igi$recod_estado_civil, 
                                      levels = c("Con pareja", "Sin pareja"))

# Tablas descriptivas 
table(base_igi$rangos_edad)
table(base_igi$recod_estado_civil)
table(base_igi$SEXO)
table(base_igi$COD_SALUD_MENTAL)

# LIMPIEZA DATOS ---

base_igi$X.8. <- car::recode(base_igi$X.8.,"2=NA")
base_igi$X.17. <- car::recode(base_igi$X.17.,"4=NA")
base_igi$X.21. <- car::recode(base_igi$X.21.,"2=NA")
base_igi$X.23. <- car::recode(base_igi$X.23.,"4=NA")
base_igi$X.27. <- car::recode(base_igi$X.27.,"4=NA")
base_igi$X.35. <- car::recode(base_igi$X.35.,"2=NA")
base_igi$X.39. <- car::recode(base_igi$X.39.,"2=NA")

# Eliminamos las variables
eliminar <- c("X.42.", "X.43.", "X.44.", "X.45.", "X.47.", "X.48.", "X.49." ,
              "X.51.", "X.52.", "X.53.", "X.54.", "X.55.", "X.56.", "X.57.",
              "X.58.") # se eliminan porque son subitems, por ejemplo
              # 42, 43, 44 y 45 son subítems de 41; 47, 48 y 49 lo son de 46, etc.

base_igi <- base_igi[,!(colnames(base_igi) %in% eliminar)]

# Cambiamos nombre de variables
descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD6","HD7","HD8",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR24","PAR25","PAR26","PAR27",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

colnames(base_igi)[8:50] <- descriptivo_grupal

# RE-CODIFICACION ---

# Nos quedamos solo con los reos con todas ítems con respuesta válida
base_igi <- base_igi %>%
  filter(if_all(all_of(descriptivo_grupal), ~ !is.na(.)))

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

verificar_binario(base_igi[,descriptivo_grupal])

# Binarizamos las variables que tienen cuatro valores. 
# El principio rector será 
# "0 = no necesidad de intervención"
# "1 = Necesidad de intervención"

recode_variables_1 <- function(x) {
  # Variables a recodificar
  vars_to_recode <- c("EDU15", "EDU16", "EDU17", "FAM18", "FAM19", "FAM20", 
                      "UTL23", "PAR25", "PAR27", "CAD30", "CAD31", "PRO36", "PRO37")
  # Verificamos que las variables existan en el data frame
  missing_vars <- setdiff(vars_to_recode, names(x))
  if (length(missing_vars) > 0) {
    warning("Las siguientes variables no existen en el data frame: ",
            paste(missing_vars, collapse = ", "))
  }
  # Re-codificamos
  for (col in vars_to_recode) {
    if (col %in% names(x)) {
      x[[col]] <- ifelse(x[[col]] %in% c(3, 2, 1), 0,  # {3,2,1} -> {0} Bajo riesgo
                         ifelse(x[[col]] %in% c(0), 1,  # {0} -> {1} Necesita intervención
                                NA))  # Para cualquier otro valor
    }
  }
  return(x)
}

base_igi_bin_1 <- recode_variables_1(base_igi)      # -->> creación base_igi_bin  !!!

recode_variables_2 <- function(x) {
  # Variables a recodificar
  vars_to_recode <- c("EDU15", "EDU16", "EDU17", "FAM18", "FAM19", "FAM20", 
                      "UTL23", "PAR25", "PAR27", "CAD30", "CAD31", "PRO36", "PRO37")
  # Verificamos que las variables existan en el data frame
  missing_vars <- setdiff(vars_to_recode, names(x))
  if (length(missing_vars) > 0) {
    warning("Las siguientes variables no existen en el data frame: ",
            paste(missing_vars, collapse = ", "))
  }
  # Re-codificamos
  for (col in vars_to_recode) {
    if (col %in% names(x)) {
      x[[col]] <- ifelse(x[[col]] %in% c(3, 2), 0,  # {3,2} -> {0} Bajo riesgo
                         ifelse(x[[col]] %in% c(0,1), 1,  # {0,1} -> {1} Necesita intervención
                                NA))  # Para cualquier otro valor
    }
  }
  return(x)
}

base_igi_bin_2 <- recode_variables_2(base_igi)      # -->> creación base_igi_bin  !!!

# Verificamos
verificar_binario(base_igi_bin_1[,descriptivo_grupal])
verificar_binario(base_igi_bin_2[,descriptivo_grupal])

# PUNTAJES DE LA MUESTRA ---

# Función para sumar los valores por fila de las variables en descriptivo_grupal
sumar_por_fila <- function(base, variables) {
  
  # Verificar si las variables están en la base de datos
  variables_presentes <- variables %in% colnames(base)
  
  if (!all(variables_presentes)) {
    stop("Algunas variables en 'descriptivo_grupal' no están presentes en la base de datos.")
  }
  
  # Seleccionamos las columnas correspondientes y sumamos por fila
  fila_sumas <- rowSums(base[ , variables], na.rm = TRUE)
  
  return(fila_sumas)
}

puntaje_total_1 <- sumar_por_fila(base_igi_bin_1, descriptivo_grupal)
puntaje_total_2 <- sumar_por_fila(base_igi_bin_2, descriptivo_grupal)

base_igi_bin_1$puntaje_total <- puntaje_total_1
base_igi_bin_2$puntaje_total <- puntaje_total_2

png("psy_net_plots/hist_puntajes_igi_1.png", width = 800, height = 600)
hist(puntaje_total_1,
     main = "Distribución de puntajes IGI {3,2,1},{0} -> {0},{1}",
     ylab = "Número de casos",
     xlab = "Puntage total IGI",
     col = "gold")
dev.off()
png("psy_net_plots/hist_puntajes_igi_2.png", width = 800, height = 600)
hist(puntaje_total_2,
     main = "Distribución de puntajes IGI {3,2},{0,1} -> {0},{1}",
     ylab = "Número de casos",
     xlab = "Puntage total IGI",
     col = "gold")
dev.off()

# MAPA DE CALOR ---

# Solo filas con puntajes
base_igi_bin_1 <- base_igi_bin_1[!is.na(base_igi_bin_1$puntaje_total), ]
base_igi_bin_2 <- base_igi_bin_2[!is.na(base_igi_bin_2$puntaje_total), ]

valores_puntaje_1 <- sort(unique(base_igi_bin_1$puntaje_total))
valores_puntaje_2 <- sort(unique(base_igi_bin_2$puntaje_total))

# Creamos una base para guardar la proporción de reos con cada variable como 1's.
base_proporciones_1 <- matrix(0, nrow = length(descriptivo_grupal), ncol = length(valores_puntaje_1))
colnames(base_proporciones_1) <- valores_puntaje_1
row.names(base_proporciones_1) <- descriptivo_grupal

# Filtramos por Puntaje (j) y por Variable de interés (i)
for (j in seq_along(valores_puntaje_1)) {
  casos_filtrados <- base_igi_bin_1$puntaje_total == valores_puntaje_1[j]
  total_reos_ptje_j <- sum(casos_filtrados)
  
  if (total_reos_ptje_j > 0) {
    for (i in seq_along(descriptivo_grupal)) {
      variable_i <- base_igi_bin_1[casos_filtrados, descriptivo_grupal[i]]
      
      # Proporción de reos con ptje j que tienen activo el ítem i
      z <- sum(variable_i, na.rm = TRUE) / total_reos_ptje_j
      
      if (!is.na(z)) {
        base_proporciones_1[i, j] <- z
      }
    }
  }
}

base_proporciones_1 <- as.data.frame(base_proporciones_1)

# Media por fila
base_proporciones_1$MediaProporciones <- rowMeans(base_proporciones_1, na.rm = TRUE)

# Creamos una columna con la categoría (HD, EDU, PAT, PAR, CAD...)
base_proporciones_1$Categoria <- gsub("[^A-Za-z]+", "", rownames(base_proporciones_1))

# Ordenamos dentro de cada categoría por la media de proporciones
base_proporciones_1 <- base_proporciones_1[order(base_proporciones_1$Categoria, -base_proporciones_1$MediaProporciones), ]

# Identificar los cambios de categoría para segmentar el plot
category_breaks_1 <- cumsum(table(base_proporciones_1$Categoria))

# Transformar matriz para ggplot2 [antes (43 obs 105 var) ahora (1763 obs 3 var)]
base_proporciones_1 <- as.data.frame(as.table(as.matrix(
  base_proporciones_1[, !(colnames(base_proporciones_1) %in% c("MediaProporciones", "Categoria"))]
)))
names(base_proporciones_1) <- c("Variable", "PuntajeTotal", "Frecuencia")

# Creamos el mapa de calor con líneas grises entre categorías
png("psy_net_plots/heatmap_prop_1.png", width = 800, height = 600)
ggplot(base_proporciones_1, aes(x = PuntajeTotal, y = Variable, fill = Frecuencia)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Proporción de activación de ítems segun puntaje total {3,2,1},{0} -> {0},{1}", 
       x = "Puntaje Total", 
       y = "Variable", 
       fill = "Proporción") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = category_breaks_1 + 0.5, color = "gold", size = 0.5, linetype = "dashed")  # Líneas grises separando categorías
dev.off()



# Creamos una base para guardar la proporción de reos con cada variable como 1's.
base_proporciones_2 <- matrix(0, nrow = length(descriptivo_grupal), ncol = length(valores_puntaje_2))
colnames(base_proporciones_2) <- valores_puntaje_2
row.names(base_proporciones_2) <- descriptivo_grupal

# Filtramos por Puntaje (j) y por Variable de interés (i)
for (j in seq_along(valores_puntaje_2)) {
  casos_filtrados <- base_igi_bin_2$puntaje_total == valores_puntaje_2[j]
  total_reos_ptje_j <- sum(casos_filtrados)
  
  if (total_reos_ptje_j > 0) {
    for (i in seq_along(descriptivo_grupal)) {
      variable_i <- base_igi_bin_2[casos_filtrados, descriptivo_grupal[i]]
      
      # Proporción de reos con ptje j que tienen activo el ítem i
      z <- sum(variable_i, na.rm = TRUE) / total_reos_ptje_j
      
      if (!is.na(z)) {
        base_proporciones_2[i, j] <- z
      }
    }
  }
}

base_proporciones_2 <- as.data.frame(base_proporciones_2)

base_proporciones_2$MediaProporciones <- rowMeans(base_proporciones_2, na.rm = TRUE)

# Creamos una columna con la categoría (HD, EDU, PAT, PAR, CAD...)
base_proporciones_2$Categoria <- gsub("[^A-Za-z]+", "", rownames(base_proporciones_2))

# Ordenamos dentro de cada categoría por la media de proporciones
base_proporciones_2 <- base_proporciones_2[order(base_proporciones_2$Categoria, -base_proporciones_2$MediaProporciones), ]

# Identificar los cambios de categoría para segmentar el plot
category_breaks_2 <- cumsum(table(base_proporciones_2$Categoria))

# Transformar matriz para ggplot2 [antes (43 obs 105 var) ahora (1763 obs 3 var)]
base_proporciones_2 <- as.data.frame(as.table(as.matrix(
  base_proporciones_2[, !(colnames(base_proporciones_2) %in% c("MediaProporciones", "Categoria"))]
)))
names(base_proporciones_2) <- c("Variable", "PuntajeTotal", "Frecuencia")

# Creamos el mapa de calor con líneas grises entre categorías
png("psy_net_plots/heatmap_prop_2.png", width = 800, height = 600)
ggplot(base_proporciones_2, aes(x = PuntajeTotal, y = Variable, fill = Frecuencia)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Proporción de activación de ítems segun puntaje total {3,2},{0,1} -> {0},{1}", 
       x = "Puntaje Total", 
       y = "Variable", 
       fill = "Proporción") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = category_breaks_2 + 0.5, color = "gold", size = 0.5, linetype = "dashed")  # Líneas grises separando categorías
dev.off()



 # PLOT REINCIDENTES VS NO-REINCIDENTES ---

base_igi$reincidencia <- ifelse(is.na(base_igi$FECHA_REINGRESO)==T, 0, 1)

# Creamos la nueva variable categórica en base a los tramos definidos
base_igi$riesgo_cat <- cut(
  base_igi_bin_1$puntaje_total,
  breaks = c(-Inf, 4, 10, 19, 29, Inf),  # Límites de los tramos
  labels = c("Muy bajo", "Bajo", "Medio", "Alto", "Muy alto"),  # Categorías
  right = TRUE  # Incluir el límite superior en cada intervalo
)

base_igi$riesgo <- cut(
  base_igi_bin_1$puntaje_total,
  breaks = c(-Inf, 4, 10, 19, 29, Inf),  # Límites de los tramos
  labels = c(1,2,3,4,5),  # Categorías
  right = TRUE  # Incluir el límite superior en cada intervalo
)
base_igi$riesgo <- as.numeric(base_igi$riesgo)

tabla_proporciones <- table(base_igi$riesgo, base_igi$reincidencia)
tabla_proporciones <- prop.table(tabla_proporciones, margin = 2)
tabla_proporciones <- as.data.frame(tabla_proporciones)
tabla_proporciones$Var1 <- as.numeric(as.character(tabla_proporciones$Var1))

# Generamos el gráfico base para Var2 == 0
png("psy_net_plots/riesgo_reinc_vs_no_reinc.png", width = 800, height = 600)
plot(
  tabla_proporciones$Var1[tabla_proporciones$Var2 == 0], 
  tabla_proporciones$Freq[tabla_proporciones$Var2 == 0], 
  type = "o", col = "blue", xlab = "Riesgo", ylab = "Proporción", 
  ylim = c(0, max(tabla_proporciones$Freq)), 
  main = "Proporciones por Riesgo y Reincidencia {3,2,1},{0} -> {0},{1}", xaxt = "n"
)
# Añadimos el eje X con etiquetas personalizadas
axis(1, at = tabla_proporciones$Var1[tabla_proporciones$Var2 == 0], 
     labels = c("Muy bajo", "Bajo", "Medio", "Alto", "Muy alto"))
# Añadimos la línea para Var2 == 1
lines(
  tabla_proporciones$Var1[tabla_proporciones$Var2 == 1], 
  tabla_proporciones$Freq[tabla_proporciones$Var2 == 1], 
  type = "o", col = "red"
)
# Agregamos una leyenda
legend("topleft", legend = c("Reincidencia = 0", "Reincidencia = 1"),
       col = c("blue", "red"), lty = 1, pch = 1)
# Agregamos una cuadrícula
grid()
dev.off()

round(cor(base_igi_bin_1$puntaje_total,base_igi$reincidencia),2)
round(cor(base_igi_bin_2$puntaje_total,base_igi$reincidencia),2)

# MAPA DE CALOR ---

# Ordenamos todas las variables por la media de proporciones (sin agrupar por categorías)
risk_category_breaks <- c(-Inf, 4, 10, 19, 29, Inf)

# Creamos el mapa de calor con líneas verticales en los breaks
ggplot(base_proporciones_1, aes(x = as.numeric(as.character(PuntajeTotal)), y = Variable, fill = Frecuencia)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Mapa de calor con categorías de riesgo", 
       x = "Puntaje Total", 
       y = "Variable", 
       fill = "Proporción") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = risk_category_breaks, color = "gold", size = 0.5, linetype = "dashed")  # Líneas verticales

# GUARDAMOS DATOS ---

# Seleccionamos variables de interés
base_igi_bin_ising_1 <-  base_igi_bin_1[,colnames(base_igi_bin_1) %in% descriptivo_grupal]
base_igi_bin_ising_1$puntaje_total <-  base_igi_bin_1$puntaje_total
base_igi_bin_ising_2 <-  base_igi_bin_2[,colnames(base_igi_bin_2) %in% descriptivo_grupal]
base_igi_bin_ising_2$puntaje_total <-  base_igi_bin_2$puntaje_total

write.csv(base_igi, "psy_net_files/base_igi.csv", row.names = FALSE)

write.csv(base_igi_bin_1, "psy_net_files/base_igi_bin_1.csv", row.names = FALSE)
write.csv(base_igi_bin_ising_1, "psy_net_files/base_igi_bin_ising_1.csv", row.names = FALSE)

write.csv(base_igi_bin_2, "psy_net_files/base_igi_bin_2.csv", row.names = FALSE)
write.csv(base_igi_bin_ising_2, "psy_net_files/base_igi_bin_ising_2.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Modelo de Ising - correlación condicionando el resto de variables
# ------------------------------------------------------------------------------

# Cargamos datos procesados

base_igi <- read.csv("psy_net_files/base_igi.csv")

base_igi_bin_1 <- read.csv("psy_net_files/base_igi_bin_1.csv")
base_igi_bin_2 <- read.csv("psy_net_files/base_igi_bin_2.csv")

base_igi_bin_ising_1 <- read.csv("psy_net_files/base_igi_bin_ising_1.csv")
base_igi_bin_ising_2 <- read.csv("psy_net_files/base_igi_bin_ising_2.csv")

# NOTA: Aquí hicimos un Ising apegado al estándar de valores 1 y -1. 
# No obstante, en etapas posteriores esto generaba errores para algunos cálculos que 
# requerían sólo valores positivos, por lo que finalmente lo dejamos con valores 0s y 1s.
# Dejamos igual el código para generar una versión con 1s y -1s por si eventualmente cambianos nuevamente. 

# Convertimos los valores 0 en valore -1 en todas las variables. 
# convert_zeros_to_neg1 <- function(df, var_names) { # Iteramos por cada nombre de variable 
#   for (col in var_names) { # Verificamos que la columna exista en el data frame 
#     if (col %in% names(df)) { # Reemplazamos los valores 0 por -1 
#       df[[col]][df[[col]] == 0] <- 0 } 
#     else { warning(paste("La columna", col, "no existe en el data frame.")) } 
#   } 
#   return(df) 
# }

#base_igi_bin_ising <- convert_zeros_to_neg1(base_igi_bin_ising, descriptivo_grupal)

#require(IsingFit)
#require(bootnet)
#require(qgraph)
#require(psychonetrics)
#install.packages("bootnet")
library(bootnet)
library(qgraph)
library(igraph)

verificar_binario(base_igi_bin_ising_1[,descriptivo_grupal]) # Check
verificar_binario(base_igi_bin_ising_2[,descriptivo_grupal]) 

# Estimamos la red
red_dicotomica_1 <- bootnet::estimateNetwork(base_igi_bin_ising_1, 
                                           default = "IsingFit",
                                           tuning = 0.25,
                                           labels = descriptivo_grupal)

# Generar el layout con qgraph y guardarlo
png("psy_net_plots/plots_redes/red_base_igi_bin_1.png", width = 1000, height = 1000)
qg <- qgraph(red_dicotomica_1$graph, layout = "spring", labels = colnames(red_dicotomica_1$graph))
dev.off()

# Extraer y guardar el layout generado
layout_guardado <- qg$layout
write.csv(layout_guardado, "psy_net_plots/plots_redes/layout_red_base_igi_bin.csv", row.names = FALSE)


# Estimamos la red
red_dicotomica_2 <- bootnet::estimateNetwork(base_igi_bin_ising_2, 
                                           default = "IsingFit",
                                           tuning = 0.25,
                                           labels = descriptivo_grupal)

# Generar el layout con qgraph y guardarlo
png("psy_net_plots/plots_redes/red_base_igi_bin_2.png", width = 1000, height = 1000)
layout_red <- as.matrix(read.csv("psy_net_plots/plots_redes/layout_red_base_igi_bin.csv"))
qg <- qgraph(red_dicotomica_2$graph, layout = layout_red, labels = colnames(red_dicotomica_2$graph))
dev.off()


# Convertimos - base_igi_bin_ising_1 - a igraph
grafo <- graph_from_adjacency_matrix(red_dicotomica_1$graph, mode = "undirected", weighted = TRUE)
E(grafo)$weight <- abs(E(grafo)$weight)

# Métricas comunes
centralidad_grado <- degree(grafo)
centralidad_intermediación <- betweenness(grafo)
centralidad_cercanía <- closeness(grafo)
cbind(centralidad_grado, centralidad_intermediación,centralidad_cercanía)

# Comunidades
comunidades <- cluster_walktrap(grafo)
plot(comunidades, grafo)

comunidades_tabla <- cbind(comunidades$names,comunidades$membership)
comunidades_tabla <- as.data.frame(comunidades_tabla)
comunidades_tabla <- comunidades_tabla[order(comunidades_tabla$V2),]
comunidades_tabla

# ------------------------------------------------------------------------------
# Estabilidad - ver psy_net_stability.R
# ------------------------------------------------------------------------------

# \\~\\

# ------------------------------------------------------------------------------
# Análisis psicométricos diferenciados por Rangos de Riesgo  
# ------------------------------------------------------------------------------

# Ahora analicemos cuáles atributos aparecen correlacionados positiva y negativamente entre sí.

# Creamos categorías de puntajes de riesgo
base_igi_bin_ising_1 <- base_igi_bin_ising_1 %>%
  dplyr::mutate(risk_category = cut(puntaje_total, 
                                    breaks = c(-Inf, 19, 29, Inf),  # Nuevos breaks
                                    labels = c("Low", "Intermediate", "High"),  # Nombres de los grupos
                                    include.lowest = TRUE))

base_igi_bin_ising_2 <- base_igi_bin_ising_2 %>%
  dplyr::mutate(risk_category = cut(puntaje_total, 
                                    breaks = c(-Inf, 19, 29, Inf),  # Nuevos breaks
                                    labels = c("Low", "Intermediate", "High"),  # Nombres de los grupos
                                    include.lowest = TRUE))

# Analizamos correlaciones por categoría de riesgo
cor_low_1 <- cor(base_igi_bin_ising_1 %>% filter(risk_category == "Low") %>% select(-risk_category))
cor_intermediate_1 <- cor(base_igi_bin_ising_1 %>% filter(risk_category == "Intermediate") %>% select(-risk_category))
cor_high_1 <- cor(base_igi_bin_ising_1 %>% filter(risk_category == "High") %>% select(-risk_category)) # HD6 no tiene varianza

cor_low_2 <- cor(base_igi_bin_ising_2 %>% filter(risk_category == "Low") %>% select(-risk_category))
cor_intermediate_2 <- cor(base_igi_bin_ising_2 %>% filter(risk_category == "Intermediate") %>% select(-risk_category))
cor_high_2 <- cor(base_igi_bin_ising_2 %>% filter(risk_category == "High") %>% select(-risk_category))


library(corrplot)
corrplot(cor_low_1, main = "Correlaciones (Puntajes Bajos)")
corrplot(cor_intermediate_1, main = "Correlaciones (Puntajes Intermedios)")
corrplot(cor_high_1, main = "Correlaciones (Puntajes Altos)")

# Para un primer acercamiento, dividamos la base en tres grupos de riesgo: bajo,intermedio y alto.

red_bajo_1 <- bootnet::estimateNetwork(
  base_igi_bin_ising_1 %>% filter(risk_category == "Low") %>% 
    select(-c(risk_category, puntaje_total)),  # Eliminar la columna de clasificación
  default = "IsingFit"
)
red_intermedio_1 <- bootnet::estimateNetwork(
  base_igi_bin_ising_1 %>% filter(risk_category == "Intermediate") %>% 
    select(-c(risk_category, puntaje_total)),  # Eliminar la columna de clasificación
  default = "IsingFit"
)
red_alto_1 <- bootnet::estimateNetwork(
  base_igi_bin_ising_1 %>% filter(risk_category == "High") %>% 
    select(-c(risk_category, puntaje_total)),  # Eliminar la columna de clasificación
  default = "IsingFit"
)

matriz_pesos_bajo_1 <- red_bajo_1$graph
matriz_pesos_medio_1 <- red_intermedio_1$graph
matriz_pesos_alto_1 <- red_alto_1$graph

red_bajo_2 <- bootnet::estimateNetwork(
  base_igi_bin_ising_2 %>% filter(risk_category == "Low") %>% 
    select(-c(risk_category, puntaje_total)),  # Eliminar la columna de clasificación
  default = "IsingFit"
)
red_intermedio_2 <- bootnet::estimateNetwork(
  base_igi_bin_ising_2 %>% filter(risk_category == "Intermediate") %>% 
    select(-c(risk_category, puntaje_total)),  # Eliminar la columna de clasificación
  default = "IsingFit"
)
red_alto_2 <- bootnet::estimateNetwork(
  base_igi_bin_ising_2 %>% filter(risk_category == "High") %>% 
    select(-c(risk_category, puntaje_total)),  # Eliminar la columna de clasificación
  default = "IsingFit"
)

matriz_pesos_bajo_2 <- red_bajo_2$graph
matriz_pesos_medio_2 <- red_intermedio_2$graph
matriz_pesos_alto_2 <- red_alto_2$graph

# Cargar el layout guardado
layout_red <- as.matrix(read.csv("psy_net_plots/plots_redes/layout_red_base_igi_bin.csv"))

# Generarmos los grafos con qgraph para inspección visual de las diferencias
png("psy_net_plots/plots_redes/red_ptje_1_bin_1_2.png", width = 1000, height = 1000)
qgraph(matriz_pesos_bajo_1, layout = layout_red, title = "Riesgo Bajo {3,2,1},{0} -> {0},{1}")
dev.off()
png("psy_net_plots/plots_redes/red_ptje_2_bin_1.png", width = 1000, height = 1000)
qgraph(matriz_pesos_medio_1, layout = layout_red, title = "Riesgo Medio {3,2,1},{0} -> {0},{1}")
dev.off()
png("psy_net_plots/plots_redes/red_ptje_3_bin_1.png", width = 1000, height = 1000)
qgraph(matriz_pesos_alto_1, layout = layout_red, title = "Riesgo Alto {3,2,1},{0} -> {0},{1}")
dev.off()

png("psy_net_plots/plots_redes/red_ptje_1_bin_2.png", width = 1000, height = 1000)
qgraph(matriz_pesos_bajo_2, layout = layout_red, title = "Riesgo Bajo {3,2},{0,1} -> {0},{1}")
dev.off()
png("psy_net_plots/plots_redes/red_ptje_2_bin_2.png", width = 1000, height = 1000)
qgraph(matriz_pesos_medio_2, layout = layout_red, title = "Riesgo Medio {3,2},{0,1} -> {0},{1}")
dev.off()
png("psy_net_plots/plots_redes/red_ptje_3_bin_2.png", width = 1000, height = 1000)
qgraph(matriz_pesos_alto_2, layout = layout_red, title = "Riesgo Alto {3,2},{0,1} -> {0},{1}")
dev.off()


# Veamos ahora si se detectan comunidades diferentes en estas tres redes obtenidas.

# # EGA para riesgo bajo
# datos_ajustados_bajo <- base_igi_bin_ising %>%
#   filter(risk_category == "Low") %>%
#   select(-risk_category) %>%
#   mutate(across(everything(), ~ . + 1))  # Mapea de -1:1 a 0:2
# 
# ega_bajo <- EGA(datos_ajustados_bajo)
# 
# # EGA para riesgo intermedio
# datos_ajustados_bajo <- base_igi_bin_ising %>%
#   filter(risk_category == "Intermediate") %>%
#   select(-risk_category) %>%
#   mutate(across(everything(), ~ . + 1))  # Mapea de -1:1 a 0:2
# 
# #ega_intermedio <- EGA(datos_ajustados_bajo)
# EGA(datos_ajustados_bajo)
# 
# # EGA para riesgo alto
# datos_ajustados_alto <- base_igi_bin_ising %>%
#   filter(risk_category == "High") %>%
#   select(-risk_category) %>%
#   mutate(across(everything(), ~ . + 1))  # Mapea de -1:1 a 0:2
# 
# #ega_intermedio <- EGA(datos_ajustados_alto)
# EGA(datos_ajustados_alto)

# Nota: Las visualizaciones sugieren que hay espacio para una exploración más detallada.


# ------------------------------------------------------------------------------
# Bootnet - Ventana de 10 puntos  
# ------------------------------------------------------------------------------

# Ahora que sabemos que parece haber diferencias en las redes para distintos 
# niveles de riesgo (como era lógico esperar), vamos al centro de nuestro objetivo: 
# veamos si para riesgos similares se producen diferencias sustantivas en las redes generadas. 
# Para ello definiremos ventanas de riesgos de 5, 10, 15 puntos y los análisis los haremos dentro de cada ventana.

# Las ventanas serán móviles. Es decir, si la ventana es de 10 puntos, se calcula para
# ventana que parte en 0 hasta 9, luego de 1 hasta 10, y así sucesivamente hasta 34 a 43.

# ---------- Formato {3,2,1},{0} -> {0},{1} ----------------------

# 1. Definimos contenedores para resultados ---
cor_results_sliding <- data.frame()
chi_results_sliding <- data.frame()
cross_tables_sliding <- list()

# Definimos las variables externas a explorar después en las regresiones
variables_externas <- c("SEXO", "COD_SALUD_MENTAL", "rangos_edad", "recod_estado_civil")

# 2. Definimos parámetros de la ventana ---
window_size <- 10 # Modificar la ventana para otros ejercicios
max_puntaje_observado <- max(base_igi_bin_1$puntaje_total)  
inicios_ventana <- 0:(max_puntaje_observado - window_size + 1)

# 3. Bucle principal por ventanas "traslapadas" ---
for (start_val in inicios_ventana) {
  
  # (a) Determinamos el rango de la ventana
  end_val <- start_val + window_size - 5  #  --->>> ¿por qué -5 y no -1?
  rango_puntaje <- start_val:end_val
  
  # (b) Filtramos la base por puntaje         -->> usamos ¿¿ base_igi ??
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
png("psy_net_plots/method_dendrogram.png", width = 800, height = 600)
plot(clust_hier)
#heights <- sort(clust_hier$height, decreasing = TRUE)
#largest_gap <- max(diff(heights))
#threshold <- heights[which(diff(heights) == largest_gap)]
abline(h = 1.75, col = "red")
dev.off() 

# Elbow Method (Within-Cluster Sum of Squares)
png("psy_net_plots/method_elbow.png", width = 800, height = 600)
wss <- sapply(1:10, function(k) { 
  sum(cutree(clust_hier, k = k)^2)
})
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "WSS") # “elbow” in the curve
dev.off()

# Silhouette Method
library(cluster)
png("psy_net_plots/method_silhouette_2.png", width = 800, height = 600)
silhouette_scores <- sapply(2:10, function(k) {
  mean(silhouette(cutree(clust_hier, k = k), distancias)[, 3])
})
plot(2:10, silhouette_scores, type = "b", xlab = "Number of Clusters", ylab = "Average Silhouette Width")
optimal_k <- which.max(silhouette_scores)
dev.off()

# Gap Statistic (To_Do : Lento - optimizar !!!! )
library(cluster)
png("psy_net_plots/method_gap.png", width = 800, height = 600)
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

# TO_DO: seguir con 00 - psy_net_recidivism3_revisar_original