
#-----------------------------------------------------------------------------


# --- REVISIÓN objeto clust_hier (rango puntaje 32 33 34 35 36 37 - 'el último clust_hier generado')

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

# --- Fin REVISIÓN



# ----------------------------------------------------------------------------


install.packages("DiscreteGapStatistic")
library(DiscreteGapStatistic)
library(cluster)

?clusGapDiscr0

# Discrete application of clusGap

# cG_obj <- clusGapDiscr0 (
#     x = sub_bloque_items,
#     FUNcluster = cluster::pam,
#     K.max = 8,
#     B = 100, # default B=nrow(x)
#     value.range = "DS", # vector with all categories. DS: Data Support option
#     verbose = interactive(), # progress output
#     distName = "hamming", # recomendado para binario
#     useLog = TRUE, # Use log function after estimating `W.k`
#     Input2Alg = 'distMatr' # input for FUNcluster
# )

cG_obj <- clusGapDiscr(
  x = sub_bloque_items,
  clusterFUN = 'pam',
  K.max = 8,
  B = 100, # default B=nrow(x)
  value.range = "DS", # vector with all categories. DS: Data Support option
  verbose = interactive(), # progress output
  distName = "hamming", # recomendado para binario
  useLog = TRUE # Use log function after estimating `W.k`
)


# Criteria to determine number of clusters k

outDF <- data.frame(cG_obj['Tab'])
outDF[, 'nClus'] <- 1:nrow(outDF)
colnames(outDF) <- c('logW', 'E.logW', 'gap', 'SE.sim', 'nClus')

myTab <- outDF %>% subset(is.finite(logW) & !is.na(gap) & !is.na(SE.sim))

selK <- with(myTab, cluster::maxSE(f = gap, SE.f = SE.sim, method = meth) )

#A numerical value from 1 to K.max, contained in the input `cG_obj` object.
myTab[selK, 'nClus']


# findK <- function (cG_obj, meth = "Tibs2001SEmax"){
#   logW <- gap <- SE.sim <- NULL
#   if (!meth %in% c("minSE", "minGap", "maxChange")) {
#     
#     ## Tibs2001 SEmax criterion
#     # myTab <- data.frame(cG_obj$Tab) %>%
#     #    mutate(nClus = 1:nrow(cG_obj$Tab) ) %>%
#     #    subset(is.finite(logW) & !is.na(gap) & !is.na(SE.sim))
#     
#     outDF <- data.frame(cG_obj['Tab'])
#     outDF[, 'nClus'] <- 1:nrow(outDF)
#     colnames(outDF) <- c('logW', 'E.logW', 'gap', 'SE.sim', 'nClus')
#     
#     myTab <- outDF %>%
#       subset(is.finite(logW) & !is.na(gap) & !is.na(SE.sim))
#     selK <- with(myTab,
#                  cluster::maxSE(f = gap,
#                                 SE.f = SE.sim,
#                                 method = meth) )
#     myTab[selK, 'nClus']
#     
#   }
#   else if (meth == "minSE") {
#     which.min(cG_obj$Tab[, "SE.sim"])
#   }
#   else if (meth == "minGap") {
#     which.min(cG_obj$Tab[, "gap"])
#   }
#   else if (meth == "maxChange") {
#     which.max(abs(cG_obj$Tab[, "gap"])) + 1
#   }
# }