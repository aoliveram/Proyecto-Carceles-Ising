# Discrete application of clusGap

cG_obj <- clusGapDiscr0 (
    x = sub_bloque_items,
    FUNcluster = cluster::pam,
    K.max = 8,
    B = 100, # default B=nrow(x)
    value.range = "DS", # vector with all categories. DS: Data Support option
    verbose = interactive(), # progress output
    distName = "hamming", # recomendado para binario
    useLog = TRUE, # Use log function after estimating `W.k`
    Input2Alg = 'distMatr' # input for FUNcluster
)

# cG_obj <- clusGapDiscr(x = sub_bloque_items, 
#                        clusterFUN = 'pam', 
#                        K.max = 8, 
#                        B = 100, # default B=nrow(x)
#                        value.range = "DS", # vector with all categories. DS: Data Support option
#                        verbose = interactive(), # progress output
#                        distName = "hamming", # recomendado para binario
#                        useLog = TRUE # Use log function after estimating `W.k`
#                        )


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