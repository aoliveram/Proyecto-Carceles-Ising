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
#    Esto se hace en el script 'datos reos - exploracion.R'

library(dplyr)
library(tidyr)
library(ggplot2)
library(bootnet)
library(cluster) # For silhouette calculation
library(viridis) # For color scales
library(qgraph)
library(igraph)
library(doParallel)

descriptivo_grupal <- c("HD1","HD2","HD3","HD4","HD5","HD7","HD8",#"HD6",
                        "EDU9","EDU10","EDU11","EDU12","EDU13","EDU14","EDU15","EDU16","EDU17",
                        "FAM18","FAM19","FAM20","FAM21",
                        "UTL22","UTL23",
                        "PAR25","PAR26","PAR27", #"PAR24",
                        "CAD28","CAD29","CAD30","CAD31","CAD32","CAD33","CAD34","CAD35",
                        "PRO36","PRO37","PRO38","PRO39",
                        "PAT40","PAT41","PAT42","PAT43")

# -------------------- Formato {3,2,1},{0} -> {0},{1} --------------------------

base_igi_bin <- read.csv("psy_net_files/base_igi_bin_1.csv")
reinc_por_reo <- readRDS("documentos carceles - ising/Datos reos - sensible/reinc_por_reo.rds")

length(unique(base_igi_bin$COD_PERS))
length(unique(reinc_por_reo$COD_PERS))

all(unique(base_igi_bin$COD_PERS) %in% unique(reinc_por_reo$COD_PERS))

# 8731 códigos están en ambas bases de datos, un 50.538 % de las disponibles en base_igi
cod_pers_baseigi_inter_reinc <- intersect(unique(base_igi_bin$COD_PERS), unique(reinc_por_reo$COD_PERS))
(length(cod_pers_baseigi_inter_reinc) / length(unique(base_igi_bin$COD_PERS))) * 100

reinc_por_reo$