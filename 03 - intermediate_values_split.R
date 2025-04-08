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
library(cluster) # For silhouette calculation
library(viridis) # For color scales
library(qgraph)
library(igraph)
library(doParallel)

# ------------------------------------------------------------------------------
# Creación objeto para 03 - intermediate_values_split.R
# ------------------------------------------------------------------------------

base_igi_bin <- read.csv("psy_net_files/base_igi_bin_1.csv")
reinc_por_reo <- readRDS("documentos carceles - ising/Datos reos - sensible/reinc_por_reo.rds")

length(unique(base_igi_bin$COD_PERS)) # 17276
length(unique(reinc_por_reo$COD_PERS)) # 97123
length(datos_reos_reinc$COD_PERS) # 203130 --> hay reos que han reincidido varias veces 

# Como hay reos que han reincidido varias veces, hay que tomar una decisión: 
#   1) tratamos cada reincidencia como un caso nuevo ('nuevo reo'), o
#   2) tratamos cada reo como único, y consideramos solo su primera reincidencia
# Acá adoptamos la opción 1).

# Seleccionamos "COD_PERS" , "SEXO", "COMUNA_DOMICILIO", "FECHA_NACIMIENTO", "DELITOS"
reos_filtrados <- datos_reos_reinc %>%
  select(COD_PERS, SEXO, COMUNA_DOMICILIO, FECHA_NACIMIENTO, DELITOS) #203130

# Filtramos las columnas que tienen delitos compuestos: más de un delito (#171755) 
reos_filtrados <- reos_filtrados %>%
  mutate(
    num_delitos = sapply(strsplit(DELITOS, "; "), length), # Columna temporal
    DELITOS = gsub(" \\(\\d+\\)", "", DELITOS)
  ) %>%
  filter(num_delitos == 1) %>%
  select(-num_delitos) # Eliminamos columna temporal

# Porcentaje de duplicados
(1 - n_distinct(reos_filtrados$COD_PERS)/nrow(reos_filtrados)) * 100

# -------------------- Creamos base_igi_bin_filtrado (solo códigos únicos) -----

all(unique(base_igi_bin$COD_PERS) %in% unique(reos_filtrados$COD_PERS))

# 7006 códigos están en ambas bases de datos, un 40.553 % de las disponibles en base_igi
# (para el resto de casos, con num_delitos > 1, hay otro 23.385 %)
cod_pers_baseigi_inter_reinc <- intersect(unique(base_igi_bin$COD_PERS), unique(reos_filtrados$COD_PERS))
(length(cod_pers_baseigi_inter_reinc) / length(unique(base_igi_bin$COD_PERS))) * 100

base_igi_bin_filtrado <- base_igi_bin %>%
  filter(COD_PERS %in% reos_filtrados$COD_PERS)

# base_igi_bin SOLO TIENE CÓDIGOS ÚNICOS, ASÍ QUE NO LO PUEDO USAR, PORQUE HAY AMBIGUEDAD 
# DE CUÁL ES LA OCASIÓN DE INGRESO DE CADA REO. DEBO VOLVER A 01 - psy_net_recidivism.R y arreglarlo.
#
# ACTUALIZACIÓN: base_igi SOLO CONTIENE CÓDIGOS ÚNICOS. HAY QUE SELECCIONAR AQUELLA FILA 
# QUE CONTIENE LA VENTANA TEMPORAL CUANDO SE HIZO EL IGI.
length(unique(base_igi_bin$COD_PERS)) == length(base_igi_bin$COD_PERS) # TRUE
length(unique(base_igi$COD_PERS)) == length(base_igi$COD_PERS) # TRUE

# OTRA OPCIÓN ES VER CUÁNTO VARÍA LA CATEGORÍA DEL DELITO DE AQUELLOS QUE TIENEN 1 ACUSACIÓN

# -------------------- ... -----------------------------------------------------

delitos_unicos_clasificados <- readRDS("psy_net_files/delitos_unicos_clasificados.rds") 

# Unir categorías al dataframe principal
reos_filtrados <- reos_filtrados %>% 
  left_join(delitos_unicos_clasificados %>% select(DELITOS, CATEGORIA_DELITO), 
            by = "DELITOS") %>% 
  mutate(CATEGORIA_DELITO = ifelse(is.na(CATEGORIA_DELITO), 
                                   "Falla categoría", 
                                   CATEGORIA_DELITO))

# 0 en 'Falla categoría' !
nrow(reos_filtrados %>% filter(CATEGORIA_DELITO == "Falla categoría")) 

# 70.6 % de los reos con solo 1 delito han reingresado más de 1 vez !
reos_filtrados_multiReinc <- reos_filtrados %>%
  group_by(COD_PERS) %>%
  filter(n() > 1) %>%  # Selecciona códigos con más de 1 registro
  ungroup()

analisis_variacion <- reos_filtrados_multiReinc %>%
  group_by(COD_PERS) %>%
  summarise(
    Total_Delitos = n(),
    Categorias_Unicas = n_distinct(CATEGORIA_DELITO),
    Lista_Categorias = toString(unique(CATEGORIA_DELITO))
  ) %>%
  arrange(desc(Categorias_Unicas))

# Distribución de categorías únicas
distribucion <- analisis_variacion %>%
  count(Categorias_Unicas, name = "N_Reos") %>%
  mutate(Porcentaje = round(N_Reos/nrow(analisis_variacion)*100, 2))

distribucion
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

min(base_igi_bin_filtrado$FECHA_ENTREVISTA)
max(base_igi_bin_filtrado$FECHA_ENTREVISTA)

# Veamos el porcentaje de entrevistas realizadas antes del egreso de los reos

# 1. Unir ambas bases por COD_PERS
datos_comparacion <- base_igi_bin_filtrado %>% 
  inner_join(datos_reos_reinc %>% select(COD_PERS, EGRESO), by = "COD_PERS")

# 2. Convertir fechas a formato Date
datos_comparacion <- datos_comparacion %>%
  mutate(
    FECHA_ENTREVISTA = as.Date(FECHA_ENTREVISTA),
    EGRESO = as.Date(EGRESO)
  )

# 3. Calcular condición y porcentaje 
datos_comparacion %>%
  summarise(
    total_casos = n(),
    entrevistas_antes = sum(FECHA_ENTREVISTA < EGRESO, na.rm = TRUE),
    porcentaje = round(entrevistas_antes / total_casos * 100, 2)
  )

# 38.67 % de las entrevistas fueron antes de EGRESO. Es decir, podemos ocupar datos de 
length(datos_reos_reinc$COD_PERS) * 0.3867
# 78550 reos.



# -------------------- RESTO DEL CÓDIGO... ---------------------------


descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD7","HD8",#"HD6",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR25","PAR26","PAR27", #"PAR24",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")