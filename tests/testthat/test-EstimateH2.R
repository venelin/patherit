library(patherit)
library(ggplot2)
library(data.table)

N <- 10000
tr <- ape::rtree(N)
z <- POUMM::rVNodesGivenTreePOUMM(tr, z0 = 3, alpha = 2, theta = 3, sigma = 1, sigmae = .4)[1:N]

# Enable parallel execution of the POUMM and PMM chains


cluster <- parallel::makeCluster(parallel::detectCores(logical = FALSE), outfile="./mcmc.log")
doParallel::registerDoParallel(cluster)

anH2 <- estimateH2(
  z, tr, 
  methods = list(PP=TRUE, 
                 POUMM=list(nSamplesMCMC=1e6, 
                            parallelMCMC = TRUE, 
                            verbose = TRUE), 
                 PMM=FALSE), verbose = TRUE)

parallel::stopCluster(cluster)



anH2$fits$POUMM$dof
summary(anH2)

plot(anH2$summaries$PP, facets=NULL) + 
  facet_grid(facets=tauQuantileType~., scales = "free")

corrTable <- corrProfile(anH2, verbose = TRUE)




tau <- seq(0, 2, by = .01)

corr <- covFunPOUMM(anH2$fits$POUMM, corr=TRUE)(tau)

corrHPDFun <- covHPDFunPOUMM(anH2$fits$POUMM, corr = TRUE, startMCMC = 5e5)  

corrHPD <- corrHPDFun(tau)

corrTablePOUMM <- data.table(tau=tau, corr = corr, 
                             corrLower = corrHPD[, 1], corrUpper = corrHPD[, 2])

corrTableSizes <- corrTable[, list(N=max(N), K = max(K)), by = tauQuantileType]

palette <- c("#000000", "#999999", "#0072B2", "#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442")

ggplot(corrTable[tauQuantileType %in% c("V", "D", "Q", "M", "A")], 
       aes(x = tauMean, y = est, 
           ymin = CI.lower, ymax = CI.upper, col = simulationMethod)) + 
  geom_pointrange(position = position_dodge(.05), pch = 20, size = .3) +
  scale_x_continuous(limits = c(0, 2)) + 
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Phylogenetic distance") +
  ylab("Correlation") +
  guides(color = "none") + 
  theme_bw() +
  stat_function(fun = covFunPOUMM(anH2$fits$POUMM, corr=TRUE)) +
  #geom_line(data = corrTablePOUMM, aes(x = tau, y = corr), inherit.aes = FALSE) +
  #geom_ribbon(data = corrTablePOUMM, 
  #            aes(x = tau, ymin = corrLower, ymax = corrUpper, alpha=0.6), 
  #            inherit.aes = FALSE) + 
  facet_grid(tauQuantileType~.) +
  geom_text(data = corrTableSizes[tauQuantileType %in% c("V", "D", "Q", "M", "A")], 
            inherit.aes = FALSE,
            aes(label=paste0(K, " PPs per bin"), 
                x = 0, y = 1, hjust = "left", vjust = "top"), 
            color = "black", size = 3)  +
  scale_color_manual(values = palette)

