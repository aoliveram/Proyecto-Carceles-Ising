######################################################
# Esto es lento. No volver a correr hasta version final
######################################################

#vamos a mirar la estabilidad de las estimaciones hechas de los valores de los links  para robustez.
estabilidad <- bootnet::bootnet(red_dicotomica, nBoots = 1000, type = "nonparametric", nCores = 5)
plot(estabilidad)

# Extraemos la información de estabilidad de edges
edges_stability <- as.data.frame(estabilidad$bootTable)

edges_filtered <- dplyr::filter(edges_stability, type == "edge")

edges_summary <- edges_filtered %>%
  dplyr::select(name, node1, node2, value, rank_avg, rank_min, rank_max) %>%
  dplyr::arrange(desc(rank_avg))

edges_summary <- edges_summary %>%
  dplyr::arrange(rank_avg)

# Agregamos una columna con la variabilidad en los rangos
edges_summary <- edges_summary %>%
  dplyr::mutate(variability = rank_max - rank_min)

ggplot(edges_summary, aes(x = rank_avg, y = variability)) +
  geom_point() +
  labs(title = "Estabilidad de Conexiones",
       x = "Rango Promedio (rank_avg)",
       y = "Variabilidad (rank_max - rank_min)") +
  theme_minimal()

edges_summary <- edges_summary %>%
  dplyr::mutate(
    classification = dplyr::case_when(
      rank_avg < 100 & variability < 50 ~ "Estable",       # Baja variabilidad y buen ranking
      variability > 100 ~ "Inestable",                     # Alta variabilidad
      rank_avg > 700 ~ "Poco relevante",                   # Alto rank_avg
      TRUE ~ "Sin clasificar"                              # Cualquier otro caso (opcional)
    )
  )

write.csv(edges_summary, paste0(aqui,"/pape/ising/estabilidad_edges_modelo_general.csv"))

conexiones_estables <- edges_summary %>%
  filter(rank_avg < 100, variability < 50)

conexiones_fuertes <- conexiones_estables %>%
  arrange(value) %>%
  head(10)

conexiones_positivas <- edges_summary %>% filter(value > 0)

ggplot(conexiones_estables, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Distribución de Valores de Conexión", x = "Valor de Conexión", y = "Frecuencia") +
  theme_minimal()
# dado que la estabilidad de las conexiones positivas es muy baja, hay que volver a probar con un segundo método de estimación el bootstrap. Como esto es lento, dejar comentado para hacerlo correr de noche
red_tmfg <- estimateNetwork(datos, default = "TMFG")
estabilidad_tmfg <- bootnet(red_tmfg, nBoots = 1000, type = "nonparametric")


#------------------

# Nota: por etapas posteriores dejamos esta exploración comentada. No es necesario volver a correr. 

# Vemaos la estabilidad de las redes/conexiones generadas. 

# Estabilidad para la red de riesgo bajo
# estabilidad_bajo <- bootnet(red_bajo, nBoots = 1000, type = "nonparametric")
# 
# # Estabilidad para la red de riesgo intermedio
# estabilidad_intermedio <- bootnet(red_intermedio, nBoots = 1000, type = "nonparametric")
# 
# # Estabilidad para la red de riesgo alto
# estabilidad_alto <- bootnet(red_alto, nBoots = 1000, type = "nonparametric")
# 
# 
# plot(estabilidad_bajo)
# plot(estabilidad_intermedio)
# plot(estabilidad_alto)