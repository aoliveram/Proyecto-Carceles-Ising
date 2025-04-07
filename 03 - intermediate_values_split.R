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
  filter(num_delitos > 1) %>%
  select(-num_delitos) # Eliminamos columna temporal


# -------------------- Formato {3,2,1},{0} -> {0},{1} --------------------------


all(unique(base_igi_bin$COD_PERS) %in% unique(reos_filtrados$COD_PERS))

# 8731 códigos están en ambas bases de datos, un 40.553 % de las disponibles en base_igi
# (para el resto de casos, con num_delitos > 1, hay otro 23.385 %)
cod_pers_baseigi_inter_reinc <- intersect(unique(base_igi_bin$COD_PERS), unique(reos_filtrados$COD_PERS))
(length(cod_pers_baseigi_inter_reinc) / length(unique(base_igi_bin$COD_PERS))) * 100

reos_filtrados$






descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD7","HD8",#"HD6",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR25","PAR26","PAR27", #"PAR24",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")