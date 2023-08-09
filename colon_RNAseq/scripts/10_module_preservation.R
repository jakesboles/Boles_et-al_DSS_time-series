#This script conducts the actual module preservation analysis

library(WGCNA)

load("colon_RNAseq/data_objects/module_preservation_input.RData")

#Distal colon vs GSE131032----
setLabels <- c("Distal", "GSE131032")
multiExpr <- list(Distal = list(data = dcolon.1), GSE131032 = list(data = dfs$gse1))
multiColor <- list(Distal = colors$gse1, GSE131032 = colors$gse1)
nSets <- 2

system.time({
  mp.d1 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Distal colon vs GSE168053 (distal)----
setLabels <- c("Distal", "GSE168053_dist")
multiExpr <- list(Distal = list(data = dcolon.2), GSE131032 = list(data = dfs$gse2_distal))
multiColor <- list(Distal = colors$gse2_distal, GSE131032 = colors$gse2_distal)
nSets <- 2

system.time({
  mp.d2 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Distal colon vs GSE210405----
setLabels <- c("Distal", "GSE210405")
multiExpr <- list(Distal = list(data = dcolon.3), GSE131032 = list(data = dfs$gse3))
multiColor <- list(Distal = colors$gse3, GSE131032 = colors$gse3)
nSets <- 2

system.time({
  mp.d3 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Proximal colon vs GSE131032----
setLabels <- c("Proximal", "GSE131032")
multiExpr <- list(Proximal = list(data = pcolon.1), GSE131032 = list(data = dfs$gse1))
multiColor <- list(Proximal = colors$gse1, GSE131032 = colors$gse1)
nSets <- 2

system.time({
  mp.p1 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Proximal colon vs GSE168053 (proximal)----
setLabels <- c("Proximal", "GSE168053_dist")
multiExpr <- list(Proximal = list(data = pcolon.2), GSE131032 = list(data = dfs$gse2_proximal))
multiColor <- list(Proximal = colors$gse2_proximal, GSE131032 = colors$gse2_proximal)
nSets <- 2

system.time({
  mp.p2 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Proximal colon vs GSE210405----
setLabels <- c("Proximal", "GSE210405")
multiExpr <- list(Proximal = list(data = pcolon.3), GSE131032 = list(data = dfs$gse3))
multiColor <- list(Proximal = colors$gse3, GSE131032 = colors$gse3)
nSets <- 2

system.time({
  mp.p3 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})

save(mp.d1, mp.d2, mp.d3,
     mp.p1, mp.p2, mp.p3,
     file = "colon_RNAseq/data_objects/module_preservation.RData")
