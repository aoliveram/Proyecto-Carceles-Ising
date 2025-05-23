
library(dplyr)
library(stringr)
library(openxlsx)

# cargamos datos limpiados de "Datos reos - exploracion.R"
delitos_unicos <- readRDS("documentos carceles - ising/Datos reos - sensible/delitos_unicos.rds")

# Cream,os nuevo objeto con columna de clasificación
delitos_unicos_clasificados <- delitos_unicos %>%
  mutate(
    # Convertimos a mayúscula, quitamos espacios inicio/final
    DELITO_LIMPIO = str_squish(toupper(DELITOS)),
    
    # Nueva columna con la categoría
    # El orden importa: las categorías más graves deberían ir primero
    # El " ^ " indica que la pista está al comienzo. El " $ ", al final.
    
    CATEGORIA_DELITO = case_when(
      
      # Homicidio y relacionados
      str_detect(DELITO_LIMPIO, "HOMICIDIO|PARRICIDIO|FEMICIDIO|INFANTICIDIO|HOM|MUER|MATAR") ~ "Homicidio y delitos relacionados",
      
      # Secuestro y sustraccion menores
      str_detect(DELITO_LIMPIO, "SECUESTRO|SUSTRACCION M|TRATA") ~ "Secuestro y Trata",
      
      #Atentados
      str_detect(DELITO_LIMPIO, "ATENTAD|EXPLOSIV|QUÍM") ~ "Atentados",
      
      # Robo (detecta cualquier forma de robo)
      str_detect(DELITO_LIMPIO, "ROB") ~ "Robo", 
      
      # Hurto
      str_detect(DELITO_LIMPIO, "HURT|SAQUEO") ~ "Hurto",
      
      # Lesiones
      str_detect(DELITO_LIMPIO, "LESION|LESIONES|LESI$") ~ "Lesiones a Terceros",
      
      # Amenazas
      str_detect(DELITO_LIMPIO, "AMENA|AMEN$") ~ "Amenazas",
      
      # Delitos de Drogas
      str_detect(DELITO_LIMPIO, "DRO|VEGET|TRAFICO|MICROTR|LEY 20000|LEY 20\\.000|ESTUPE") ~ "Delitos de Drogas",
      
      # Delitos Sexuales
      str_detect(DELITO_LIMPIO, "VIOLAC|ABUSO S|DELITO SEX|SEXU|PORN|PROSTI") ~ "Delitos Sexuales",
      
      # Porte / Tenencia de Armas
      str_detect(DELITO_LIMPIO, "POSESIÓN, T|PORTE ILEGAL|TENENCIA ILEGAL|DE ARMA") ~ "Porte/Tenencia de Armas",
      
      # Daños
      str_detect(DELITO_LIMPIO, "DAÑO|DAÑOS") ~ "Daños a la Propiedad",
      
      # Delitos de Tránsito (principalmente conducción ebriedad)
      str_detect(DELITO_LIMPIO, "CONDUC. E|ESTADO D|CONDUCCION EBRIEDAD|ALCO|BAJO INFLUEN|TRANSIT") ~ "Delitos de Tránsito",
      
      # Estafas / Delitos económicos
      str_detect(DELITO_LIMPIO, "ESTAF|APROPIACION INDEBIDA|FRAUD|ADULTERAC|FALSIF|LAVAD") ~ "Delitos Económ/Estafas/Falsific",
      
      # Receptación
      str_detect(DELITO_LIMPIO, "RECEPTA") ~ "Receptación",
      
      # Desacato / Oposición
      str_detect(DELITO_LIMPIO, "DESACATO|OPOSICION|RESISTENCIA|RIÑA|BUENAS COS|ORDEN P|DESORD") ~ "Oposición a la Autoridad y Riña",
      
      # Violencia Intrafamiliar
      str_detect(DELITO_LIMPIO, "VIOLENCIA IN|INTRAF|INTRA F|ABANDONO DE NI|MALTRATO INF") ~ "Violencia Intrafamiliar (VIF)",
      
      # Cohecho
      str_detect(DELITO_LIMPIO, "COHEC") ~ "Cohecho",
      
      # Multa - Pena sustitutiva
      str_detect(DELITO_LIMPIO, "^MULT") ~ "Multa -Pena sustitutiva",
      
      # Si no hay coincidencias, "Otros"
      TRUE ~ "Otros Delitos / No Especificado"
    )
  )

#resumen de cuántos casos caen en cada categoría
print(count(delitos_unicos_clasificados, CATEGORIA_DELITO, sort = TRUE))

# Guardar en Excel 
write.csv(delitos_unicos_clasificados, "documentos carceles - ising/Datos reos - sensible/delitos_unicos_clasificados.csv", row.names = FALSE)
write.xlsx(delitos_unicos_clasificados, "documentos carceles - ising/Datos reos - sensible/delitos_unicos_clasificados.xlsx", overwrite = TRUE)
# Guardar en formato RDS ('Datos reos - sensible' y 'psy_net_files')
saveRDS(delitos_unicos_clasificados, "documentos carceles - ising/Datos reos - sensible/delitos_unicos_clasificados.rds")
saveRDS(delitos_unicos_clasificados, "psy_net_files/delitos_unicos_clasificados.rds")

# Cargar
delitos_unicos_clasificados <- readRDS("documentos carceles - ising/Datos reos - sensible/delitos_unicos_clasificados.rds")

# 
table(delitos_unicos_clasificados$CATEGORIA_DELITO)
