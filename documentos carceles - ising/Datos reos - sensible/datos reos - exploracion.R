# Algunos datos:

# 0) 712626 observaciones
# 1) 346244 son CONDENADOS
# 2) 546721 son CONDENADOS, con múltiple aparicion de COD_PERS
# 3) 338786 son CONDENADOS, con múltiple aparicion de COD_PERS, diferentes de TRASLADO
# 4) 203130 son CONDENADOS, con múltiple aparicion de COD_PERS, diferentes de TRASLADO, INGRESO > EGRESO
# 5) 97123 son CONDENADOS, con múltiple aparicion de COD_PERS, diferentes de TRASLADO, INGRESO > EGRESO, y CÓDIGO ÚNICO

# 1868 DELITOS ÚNICOS al descomponer (AUN MUY SUCIO) EJ: 
# FAL
# FALS
# FALSIF
# FALSIFICA
# FALSIFICAC
# FALSIFICACION
# FALSIFICACION (1)
# FALSIFICACION (6)
# FALSIFICACIÓN DE BILLETES ART. 64 LEY ORGÂNICA

# ------------------------------------------------------------------------------

library(readxl)
library(tidyr)
library(dplyr)

# Leer la primera hoja del archivo
datos_reos <- read_excel("documentos carceles - ising/Datos reos - sensible/BD1 INGRESOS.xlsx")

str(datos_reos)

# filtrar por códigos únicos
# reos_codigos <- unique(datos_reos$COD_PERS) # 302215 codigos unicos 

#filtrar por CONDENADOS
# datos_reos$`CALIDAD_PROCESAL AL INGRESO` == "CONDENADO"

# usar esos filtros para obtener los delitos ya filtrados
# datos_reos$DELITOS
# length(unique(datos_reos$DELITOS)) # 40658 entradas únicas
# length(datos_reos$DELITOS)-length(unique(datos_reos$DELITOS))

# datos_reos$CONDENA # prescindible

unique(datos_reos$...15) # ¿qué es esto? es como motivo por el que ya no está??
# sum(datos_reos$...15 == "FUGA (009)[C]", rm.na=TRUE)

# 1) --------Filtramos por CONDENADO--------------------------------------------

datos_reos_reinc <- datos_reos %>%
  filter(`CALIDAD_PROCESAL AL INGRESO` == "CONDENADO")

# 2) --------Filtramos por al menos dos apariciones del mismo COD_PERS----------

cod_frecuencia <- datos_reos %>%
  group_by(COD_PERS) %>%
  tally() %>%
  filter(n > 1)  # Filtrar solo los que aparecen más de una vez

datos_reos_reinc <- datos_reos %>%
  filter(COD_PERS %in% cod_frecuencia$COD_PERS)

length(unique(datos_reos_reinc$COD_PERS)) # 136309 (vs 302215): muchos no son condenados

# # usar esos filtros para obtener los delitos ya filtrados
# datos_reos_reinc$DELITOS
# length(unique(datos_reos_reinc$DELITOS)) # 22914 (vs 40658): solo algunas entradas son de condenados
# length(datos_reos_reinc$DELITOS)-length(unique(datos_reos_reinc$DELITOS))
# 
# # Hay casi 1/4 de número de Delitos en comparación a los reos con ¿reincidencia??.
# length(unique(datos_reos_reinc$DELITOS)) / length(unique(datos_reos_reinc$COD_PERS))

# # Aseguremos si es que HAY o no REINCIDENCIA
# head(datos_reos_reinc$`F_INGRESO AL EP`)

# 3) --------Filtramos por TRASLADO---------------------------------------------

# Filtrar datos eliminando los traslados. 163656 original vs 82514 nuevo
datos_reos_reinc <- datos_reos_reinc %>%
  filter(`TIPO DE INGRESO` != "TRASLADO")

# 4) --------Filtramos por reos con INGRESOS posteriores a EGRESOS--------------

# Convertir formatos de FECHAS ingreso / egreso
datos_reos_reinc <- datos_reos_reinc %>%
  mutate(`F_INGRESO AL EP` = as.Date(`F_INGRESO AL EP`),
         EGRESO = as.Date(as.numeric(EGRESO), origin = "1899-12-30"))

head(datos_reos_reinc$`F_INGRESO AL EP`)
head(datos_reos_reinc$EGRESO)

#Ordenar por código de reo y fecha de ingreso
datos_reos_reinc <- datos_reos_reinc %>%
  arrange(COD_PERS, `F_INGRESO AL EP`)

# Identificar reincidencias. 261029 original vs 163656 nuevo - 203130
datos_reos_reinc <- datos_reos_reinc %>%
  group_by(COD_PERS) %>%
  mutate(reincide = `F_INGRESO AL EP` > lag(EGRESO, default = as.Date(NA))) %>%
  filter(reincide == TRUE) %>%  # Mantener solo reincidentes
  ungroup()

# Guardar en formato RDS
saveRDS(datos_reos_reinc, "documentos carceles - ising/Datos reos - sensible/datos_reos_reinc.rds")
datos_reos_reinc <- readRDS("documentos carceles - ising/Datos reos - sensible/datos_reos_reinc.rds")

# --------Vemos los DELITOS UNICOS----------------------------------------------

head(datos_reos_reinc$DELITOS, 20)

# # Separar delitos en múltiples filas --> 2426 delitos únicos
# delitos_unicos <- datos_reos_reinc %>%
#   separate_rows(DELITOS, sep = "; ") %>%  # Divide las entradas compuestas
#   #mutate(DELITOS = gsub(" \\(\\d+\\)", "", DELITOS)) %>%  # Elimina los números entre paréntesis
#   distinct(DELITOS)  # Obtenemos solo delitos únicos
# delitos_unicos <- delitos_unicos %>% arrange(DELITOS)
# 
# head(delitos_unicos)
# 
# # Guardar en CSV
# #write.csv(delitos_unicos, "documentos carceles - ising/Datos reos - sensible/delitos_unicos.csv", row.names = FALSE)
# # Guardar en formato RDS
# saveRDS(delitos_unicos, "documentos carceles - ising/Datos reos - sensible/delitos_unicos.rds")
# delitos_unicos <- readRDS("documentos carceles - ising/Datos reos - sensible/delitos_unicos.rds")

# Sin paréntesis
delitos_unicos <- datos_reos_reinc %>%
  separate_rows(DELITOS, sep = "; ") %>%  # Divide las entradas compuestas
  mutate(DELITOS = gsub(" \\(\\d+\\)", "", DELITOS)) %>%  # Elimina los números entre paréntesis
  distinct(DELITOS)  # Obtenemos solo delitos únicos
delitos_unicos <- delitos_unicos %>% arrange(DELITOS)

# Guardar en CSV
write.csv(delitos_unicos, "documentos carceles - ising/Datos reos - sensible/delitos_unicos.csv", row.names = FALSE)
# Guardar en formato RDS
saveRDS(delitos_unicos, "documentos carceles - ising/Datos reos - sensible/delitos_unicos.rds")
delitos_unicos <- readRDS("documentos carceles - ising/Datos reos - sensible/delitos_unicos.rds")


# 5) --------Nuevo objeto con UNA FILA POR REO-------------------------------------

length(unique(datos_reos_reinc$COD_PERS))

# Contar número de ingresos por reo
reinc_por_reo <- datos_reos_reinc %>%
  group_by(COD_PERS) %>%
  summarise(n_ingresos = n(),  # Número de veces en prisión
            delitos = list(DELITOS))  # Lista con todos sus delitos

# Convertir lista de delitos a columnas separadas
reinc_por_reo <- reinc_por_reo %>%
  unnest_wider(delitos, names_sep = "_")

# Guardar en CSV
#write.csv(reinc_por_reo, "documentos carceles - ising/Datos reos - sensible/reinc_por_reo.csv", row.names = FALSE)
# Guardar en formato RDS
saveRDS(reinc_por_reo, "documentos carceles - ising/Datos reos - sensible/reinc_por_reo.rds")
reinc_por_reo <- readRDS("documentos carceles - ising/Datos reos - sensible/reinc_por_reo.rds")
