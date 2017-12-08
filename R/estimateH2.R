# Copyright 2017 Venelin Mitov
#
# This file is part of patherit.
# 
# patherit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# patherit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#' Estimate the heritability of a continuous trait in a phylogenetically 
#' linked population. 
#' 
#' @description This function estimates the broad-sense heritability of a trait
#' from trait measurements at the tips of a phylogenetic tree connecting the carriers
#' of this trait. By default the function performs a phylogenetic pair (PP) analysis,
#' a phylogenetic mixed model fit (PMM) and a phylogenetic Ornstein-Uhlenbeck 
#' mixed model fit (POUMM) and returns an object of class H2Analysis. Use the 
#' generic functions summary and plot to see the heritability estimates.
#' 
#' @param z A numeric vector of length equal to the number of tips in tree matching
#' their order.
#' @param tree A phylo object. 
#' @param methods A named list of logical values and/or lists, specifying each 
#' estimating method. The names of the elements should start with PP, POUMM,
#' or PMM. For example, specifying \code{methods = list(PP = TRUE, POUMM_1 = TRUE, 
#' POUMM_2 = list(nSamplesMCMC = 1e6), PMM = TRUE)} will cause the methods PP and 
#' PMM to be executed with default settings, while the method POUMM will be 
#' executed with a million iterations in each MCMC chain. 
#' The values in the list methods should be either logical or lists.
#' To disable the fit of a method, simply exclude it from the list or set it to FALSE.
#' 
#' @param verbose Logical: should messages be written on the screen at runtime.
#' @details The order of elements in z should match the order of the tips in tree. 
#' By default the function will run six MCMC chains of length 100,000. Without 
#' parallelization of the MCMC chains and the likelihood calculation, this may 
#' take more than 20 minutes to run on a tree of 1000 tips.
#' 
#' @return an S3 object of class H2Analysis. Use the generic functions summary
#' and plot to display the results.
#'
#' @import POUMM
#' 
#' @seealso \code{\link{PP}}, \code{\link{POUMM::POUMM}}, \code{\link{POUMM::specifyPOUMM}}.
#'  
#' @examples 
#' \dontrun{
#' # Please, read the package vignette for more detailed examples.
#' N <- 200
#' tr <- ape::rtree(N)
#' z <- POUMM::rVNodesGivenTreePOUMM(tr, 0, 2, 3, 1, 1)[1:N]
#' 
#' # Enable parallel execution of the POUMM and PMM chains
#' 
#' cluster <- parallel::makeCluster(parallel::detectCores(logical = FALSE))
#' doParallel::registerDoParallel(cluster)
#' 
#' H2Analysis <- estimateH2(z, tree)
#' 
#' parallel::stopCluster(cluster)
#' 
#' summary(H2Analysis)
#' }
#' 
#'@references 
#'  Mitov, V. and Stadler, T. (2016). The heritability of pathogen traits - 
#'  definitions and estimators. bioRxiv, 058503
#'  doi: https://doi.org/10.1101/058503
#'  
#'  Mitov, V., and Stadler, T. (2017). Fast and Robust Inference of Phylogenetic 
#'  Ornstein-Uhlenbeck Models Using Parallel Likelihood Calculation. bioRxiv, 115089. 
#'  https://doi.org/10.1101/115089
#'  
#' @export
estimateH2 <- function(z, tree, 
                       methods = list(PP = TRUE, PMM = TRUE, POUMM = TRUE), 
                       verbose = FALSE) {
  
  if(!(is.numeric(z) & class(tree) == 'phylo' & 
       length(z) == length(tree$tip.label)) ) {
    stop("z should be a numeric vector of length equal to the number of tips in tree; tree should be a phylo object. The order of elements in z should match the order of the tips in tree.")
  } else 
    
    if(any(is.na(z))) {
      stop("There are NA or NaN elements in z.")
    }
  
  zVar <- var(z)
  zSD <- sqrt(zVar)
  zMin <- min(z)
  zMean <- mean(z)
  zMax <- max(z)
  tMin <- min(nodeTimes(tree, tipsOnly = TRUE))
  tMax <- max(nodeTimes(tree, tipsOnly = TRUE))
  tMean <- mean(nodeTimes(tree, tipsOnly = TRUE))
  
  res <- list(z = z, tree = tree, fits = list(), summaries = list(), plots = list())
  
  for(i in 1:length(methods)) {
    method <- methods[[i]]
    
    if( ! (is.list(method) | identical(method, TRUE)))  {
      next
    }
    if(verbose) {
      cat("Performing ", names(methods)[i], '(', toString(method), ")\n")
    }
    
    if(startsWith(names(methods)[i], "PP")) {
      specPP <- methods[[i]]
      specPP2 <- list(z = z, tree = tree)
      
      if( is.list(specPP) ) {
        for(name in names(specPP)) {
          if(! name %in% 'verbose')
            specPP2[[name]] <- specPP[[name]]
        }
      }
      if(verbose) {
        cat("Performning", names(methods)[i], "analysis...\n")
      }
      
      fit <- do.call(PP, specPP2)
      
      res$fits[[names(methods)[i]]] <- fit
      res$summaries[[names(methods)[i]]] <- summary(fit)
      
      
      if(verbose) {
        cat("Summary from", names(methods)[i], "analysis...\n")
        print(res$summaries[[names(methods)[i]]])
      }
      
    } else if(startsWith(names(methods)[i], "POUMM")) {
      specPOUMM <- methods[[i]]
      specPOUMM2 <- list(
        z = z, tree = tree, zMin = zMin, zMean = zMean, zMax = zMax, 
        zVar = zVar, zSD = zSD, tMin = tMin, tMean = tMean, tMax = tMax, 
        parallelMCMC = TRUE)
      
      verbosePOUMM <- FALSE
      if( is.list(specPOUMM) ) {
        verbosePOUMM <- if( is.logical(specPOUMM$verbose) ) specPOUMM$verbose else FALSE
        
        for(name in names(specPOUMM)) {
          if(! name %in% 'verbose')
            specPOUMM2[[name]] <- specPOUMM[[name]]
        }
      }
      
      if(verbose) {
        cat("Validating", names(methods)[i], "specification\n.")
      }
      specPOUMM2 <- do.call(specifyPOUMM_ATH2tMeanSeG0, specPOUMM2)
      
      if(verbose) {
        cat("Performing", names(methods)[i], "analysis...\n")
      }
      
      fit <- POUMM(z, tree, spec = specPOUMM2, verbose = verbosePOUMM)
      
      res$fits[[names(methods)[i]]] <- fit
      
      res$summaries[[names(methods)[i]]] <- list(short = summary(fit), expert = summary(fit, mode="expert"))
      
      if(verbose) {
        cat("Summary from", names(methods)[i], "analysis...\n")
        print(res$summaries[[names(methods)[i]]]$short)
      }
      
    } else if(startsWith(names(methods)[i], "PMM")) {
      specPMM <- methods[[i]]
      specPMM2 <- list(
        z = z, tree = tree, zMin = zMin, zMean = zMean, zMax = zMax, 
        zVar = zVar, zSD = zSD, tMin = tMin, tMean = tMean, tMax = tMax, 
        parallelMCMC = TRUE)
      
      verbosePMM <- FALSE
      if( is.list(specPMM) ) {
        verbosePMM <- if( is.logical(specPMM$verbose) ) specPMM$verbose else FALSE
        
        for(name in names(specPMM)) {
          if(! name %in% 'verbose') 
            specPMM2[[name]] <- specPMM[[name]]
        }
      }
      
      if(verbose) {
        cat("Validating", names(methods)[i], "specification\n.")
      }
      specPMM2 <- do.call(specifyPMM_H2tMeanSeG0, specPMM2)
      
      if(verbose) {
        cat("Performing PMM analysis...\n")
      }
      
      fit <- POUMM(z, tree, spec = specPMM2, verbose = verbosePMM)
      
      res$fits[[names(methods)[i]]] <- fit
      
      res$summaries[[names(methods)[i]]] <- 
        list(short = summary(fit), expert = summary(fit, mode="expert"))
      
      if(verbose) {
        cat("Summary from", names(methods)[i], "analysis...\n")
        print(res$summaries[[names(methods)[i]]]$short)
      }
    }
  }
  
  class(res) <- "H2Analysis"
  res
}

#' @export
summary.H2Analysis <- function(object, useBootstrapsForPP = FALSE, ...) {
  
  if("H2Analysis" %in% class(object)) {
    summaryAll <- NULL
    for(methodName in names(object$fits)) {
      if(startsWith(methodName, "PP")) {
        smm <- 
          summary(object$fits[[methodName]])[
            stat == "rA" & grepl("\\[0%", x = tauQuantile),
            list(method = methodName, tauQuantile, 
                 N, K, stat, MLE = est, PostMean = NA,
                 CI.lower = if(useBootstrapsForPP) bCI.lower else CI.lower, 
                 CI.upper = if(useBootstrapsForPP) bCI.upper else CI.upper, 
                 tauMean, filter)]
        
      } else if(any(startsWith(methodName, c("POUMM", "PMM")))) {
        smm <- summary(object$fits[[methodName]], ...)
        if("PostMean" %in% names(smm)) {
          smm <- smm[
            stat %in% c("H2tMean", "H2e", "H2tInf", "AIC", "AICc", 
                        "alpha", "theta", "sigma", "sigmae", "g0"),
            list(method = methodName,tauQuantile = NA, N, K = NA, stat, MLE,
                 PostMean, 
                 CI.lower = sapply(HPD, function(.) .[1]),
                 CI.upper = sapply(HPD, function(.) .[2]),
                 tauMean = NA, filter = NA)]
        } else {
          smm <- smm[
            stat %in% c("H2tMean", "H2e", "H2tInf", "AIC", "AICc",  
                        "alpha", "theta", "sigma", "sigmae", "g0"),
            list(method = methodName,tauQuantile = NA, N, K = NA, stat, MLE,
                 PostMean = NA, 
                 CI.lower = NA,
                 CI.upper = NA,
                 tauMean = NA, filter = NA)]
        }
      }
      summaryAll <- rbind(summaryAll, smm)
    }
    summaryAll
  } else {
    stop("summary.H2Analysis called on a non H2Analysis-object.")
  }
}

#' Generate a correlation profile for a given tree and data at the tips
#' 
#' The correlation profile for a given tree and trait values at its tips 
#' represents the values of the phylogenetic pair (PP) correlation among PPs 
#' stratified patristic distance. This function allows to compare the 
#' correlation profile emerging from the original data with correlation profiles
#' emerging from data simulated along the tree according to maximum likelihood
#' fits of the POUMM and PMM methods (all ML-fits found in object$fits).
#' 
#' @param object an S3 object of class H2Analysis
#' @param nSim an Integer indicating the number of simulations per ML-fit
#' @param seed optional a seed from the random number generator.
#' @param verbose Logical. Should informative messages be output during the simulations.
#' @param ... Additional parameters passed to summary.PP, e.g. tauQuantileType = "D"
#' 
#' @importFrom data.table copy
#' 
#' @export
corrProfile <- function(object, nSim = 100, seed = NA, verbose = FALSE, ...) {
  if("H2Analysis" %in% class(object)) {
    if(is.null(object$fits$PP)) {
      stop("No PP-analysis included. Possibly, you disabled the PP-analysis in the call to estimateH2?")
    }
    
    pp <- copy(object$fits$PP$pp)
    
    if('z' %in% names(pp)) {
      pp[, z:=NULL]
      pp[, deltaz:=NULL]
    }
    
    corrTable <- summary(object$fits$PP, ...)[stat=="rA"]
    corrTable[, simulationMethod:="Original"]
    corrTable[, MLE:=list(list(NA))]
    corrTable[, simNo:=NA]
    
    corrTable <- corrTable[
      ,
      list(tauQuantile, stat, N, K, 
           est, filter, 
           tauMean, tauMedian, 
           tMean, tMedian,
           tauQuantileType, 
           CI.lower, 
           CI.upper,
           Data = simulationMethod,
           method = "PP",
           MLE)
      ]
    
    if(!is.na(seed)) {
      set.seed(seed)
    }
    
    simulatedData <- list()
    
    for(methodName in names(object$fits)) {
      if(any(startsWith(methodName, c('POUMM', 'PMM')))) {
        
        simulatedDataMethod <- list()
        
        fit <- object$fits[[methodName]]
        par <- fit$spec$parMapping(coef(fit))
        if(is.na(par['g0'])) {
          par['g0'] <- mean(object$z)
        }
        if(verbose) {
          cat("Performing", nSim, methodName, "simulations, MLE:", 
              toString(round(par,2)), "\n")
        }
        
        corrTableMethod <- NULL
        
        for(i in 1:nSim) {
          zSim <- rVNodesGivenTreePOUMM(
            object$tree, par['g0'], par['alpha'], par['theta'], par['sigma'],
            par['sigmae'])[1:length(object$tree$tip.label)]
          
          simulatedDataMethod[[i]] <- zSim
          
          PPSim <- PP(zSim, tree = NULL, dists = NULL, pp = pp, tauQuantiles = object$fits$PP$tauQuantiles)
          smm <- summary(PPSim, ...)[stat=='rA']
          smm[, simulationMethod:=methodName]
          smm[, MLE:=list(list(par))]
          smm[, simNo:=i]
          corrTableMethod <- rbindlist(list(corrTableMethod, smm))
          if(verbose) {
            cat('.')
          }
        }
        cat("done\n")
        
        simulatedData[[methodName]] <- simulatedDataMethod
        
        corrTableMethod <- corrTableMethod[
          ,
          list(stat = unique(stat), N = unique(N), K = unique(K), 
               est = mean(est), filter = unique(filter), 
               tauMean = unique(tauMean), tauMedian = unique(tauMedian), 
               tMean = unique(tMean), tMedian = unique(tMedian),
               tauQuantileType = unique(tauQuantileType), 
               CI.lower = quantile(est, probs = .025, na.rm = TRUE), 
               CI.upper = quantile(est, probs = 0.975, na.rm = TRUE),
               Data = paste0(unique(simulationMethod), " simulation"),
               method = unique(simulationMethod), 
               MLE = MLE[1]),
          by = tauQuantile
          ]
        
        corrTable <- rbind(corrTable, corrTableMethod)
      }
    }
    list(corrTable=corrTable, simulatedData = simulatedData)
  } else {
    stop("object should be of S3 class 'H2Analysis'.")
  }
}

#' Plot a correlation profile
#' 
#' @param object an S3 object of class "H2Analysis"
#' @param corrTable a data.table returned by corrProfile
#' @param tauQuantileTypes a character vector of at least 1 element denoting the 
#' types of quantiles of tau (stratifications) to plot the correlation profile
#' for. All these types must already be available in the corrTable object.
#' 
#' @importFrom ggplot2 ggplot geom_pointrange geom_label scale_x_continuous aes xlab ylab stat_function  facet_grid position_dodge scale_color_manual scale_y_continuous coord_cartesian
#' @export
plotCorrProfile <- function(object, 
                            corrProfile = corrProfile(object, tauQuantileTypes = tauQuantileTypes), 
                            tauQuantileTypes = c("A", "M", "Q", "D", "V"), 
                            palette = c("#999999", "#0072B2", "#CC79A7", 
                                        "#E69F00", "#D55E00", "#56B4E9", 
                                        "#009E73", "#F0E442"), 
                            xlim = NULL, ylim = NULL, 
                            showLabel = TRUE,
                            dodgeRatio = 1/25,
                            labelSize = 2.5, labelColor = "black") {
  if("H2Analysis" %in% class(object)) {
    if(is.null(object$fits$PP)) {
      stop("No PP-analysis included. Possibly, you disabled the PP-analysis in the call to estimateH2?")
    }
    
    force(corrProfile)
    
    corrTable <- corrProfile$corrTable
    
    fitPP <- object$fits$PP
    
    corrTable[, tauQuantileType:=factor(tauQuantileType, 
                                        levels = rev(levels(tauQuantileType)))]
    
    corrTable <- corrTable[tauQuantileType %in% tauQuantileTypes]
    
    maxTau <- corrTable[, max(tauMean)]
    minTau <- corrTable[, min(tauMean)]
    
    if(is.null(xlim)) {
      xlim <- c(min(0, minTau - maxTau*dodgeRatio), maxTau + maxTau*dodgeRatio)  
    }
    if(is.null(ylim)) {
      ylim <- c(0, round(min(1, corrTable[, max(CI.upper)]) + .05, 1)  )
    }
    
    corrTableSizes <- corrTable[Data == "Original", 
                                list(N=max(N), Nbins = length(K), 
                                     K = K[1], K2 = round(mean(K))), 
                                by = tauQuantileType]
    
    if(is.null(names(palette))) {
      Datas <- corrTable[, unique(Data)]
      names(palette)[1:length(Datas)] <- Datas
    }
    
    pl <- ggplot() + 
      geom_pointrange(
        data = corrTable, 
        mapping = aes(x = tauMean, y = est, 
                       ymin = CI.lower, ymax = CI.upper, col = Data),
        pch = 20, size = .25, 
        position = position_dodge((xlim[2] - xlim[1])*dodgeRatio)) +
      #scale_x_continuous(limits = xlim) +
      coord_cartesian(xlim = xlim, ylim = ylim) + 
      xlab(expression(paste("Phylogenetic distance (", tau, ")"))) +
      scale_color_manual(values = palette) + 
      ylab("Correlation")
    
    # for(methodName in names(object$fits)) {
    #   if(any(startsWith(methodName, c('POUMM', 'PMM')))) {
    #     fit <- object$fits[[methodName]]
    #     smmfit <- summary(fit)
    #     methodDataName <- paste0()
    #     pl <- pl + stat_function(
    #       data = corrTable, 
    #       mapping = aes(x = tauMean, y = est, col = Data),
    #       fun = covFunPOUMM(fit, corr=TRUE), 
    #       col=palette[grep(methodName, names(palette))[1]])
    #   }
    # }
    
    if(showLabel) {
      pl <- pl +
        geom_label(
          data = corrTableSizes[tauQuantileType %in% tauQuantileTypes], 
          inherit.aes = FALSE,
          aes(label = paste0(ifelse(Nbins == 1, 
                                    paste0("All ", K, " PPs"),
                                    paste0(K, " CPPs, avgerage ", 
                                           K2, " PPs per stratum"))), 
              x = 0, y = ylim[2], hjust = "left", vjust = "top"), 
          color = labelColor, size = labelSize) 
    }
    
    if(length(tauQuantileTypes) > 1) {
      pl + facet_grid(tauQuantileType~.)
    } else {
      pl
    }
  } else {
    stop("object should be of S3 class 'H2Analysis'.")
  }
}

#' Regression slope of recipient on donor values in an epidemic
#' @param sampledOnly logical, indicating if only sampled individuals should be 
#' included in the heritability calculation
#' @param tMin,tMax numeric, time interval for which to measure the heritability. 
#' if sampledOnly is TRUE this would be the times of sampling for the 
#' individuals; otherwise, this would be the times of infection. 
#' @param corr logical, should donor-recipient correlation be returned instead 
#' of regression slope
#' 
#' @import data.table
#' 
#' @export
b <- function(epidemic=NULL, data=NULL, GEValues, sampledOnly=TRUE, 
              activeOnly=FALSE, tMin=0, tMax=Inf, 
              firstN=Inf, lastN=Inf, atInfection=FALSE, report=FALSE, 
              corr=FALSE) {
  if(is.null(epidemic)&is.null(data)) {
    warning('One of the parameters epidemic or data should be specified, but both were NULL. Returning NA.')
    NA
  } else if(!is.null(data)) {
    couples <- data
  } else {
    couples <- extractDRCouples(epidemic=epidemic, 
                                sampledOnly=sampledOnly, activeOnly=activeOnly, 
                                tMin=tMin, tMax=tMax, firstN=firstN, lastN=lastN)
  }
  
  if(nrow(couples) > 0) {
    couples[, GED0:=GEValues[cbind(envd, gD0)]]
    couples[, GER0:=GEValues[cbind(env, gD)]]
    couples[, zD0:=calcValue(envd, gD0, eD0, GEValues)] 
    couples[, zR0:=calcValue(env, gR0, eR0, GEValues)]
    
    couples[, GED:=GEValues[cbind(envd, gD)]]
    couples[, GER:=GEValues[cbind(env, gR)]]
    couples[, zD:=calcValue(envd, gD, eD, GEValues)]
    couples[, zR:=calcValue(env, gR, eR, GEValues)]
    
    beta0=couples[, cov(zD0, zR0)/var(zD0)]
    beta=couples[, cov(zD, zR)/var(zD)]  
    
    cor0=couples[, cor(zD0, zR0)]
    cor=couples[, cor(zD, zR)]
    
    if(report) {
      list(beta0, beta=beta, cor=cor, cor0=cor0, couples=couples)
    } else {
      if(atInfection) {
        ifelse(corr, cor0, beta0)
      } else {
        ifelse(corr, cor, beta)
      }
    }
  } else {
    NA
  }           
}


#' Estimating intraclass correlation (ICC) through ANOVA
#' @param epidemic a list of objects returned from simulateEpidemic
#' @param data NULL or a data.table such as the element gen in an epidemic list
#' @param GEValues genotype-environment trait values for calculating the z-values; if NULL it is assumed that the 
#' z-values are already in data, or epidemic$gen.
#' @param sampledOnly,activeOnly,tMin,tMax,firstN,lastN parameters passed to extractPop
#' @param by a character string which can be evaluated as expression in the by clause of data.table
#' the times of sampling for the individuals; otherwise, this would be the times of infection.
#' @param report a logical indicating what result should be returned. If FALSE, only the ICC value is returned, 
#' otherwise a list of statistics from the ANOVA calculation
#' @return See report parameter
#' 
#' @import data.table
#' 
#' @export
rA <- function(epidemic=NULL, data=NULL, GEValues=NULL, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf,
               firstN=Inf, lastN=Inf, by=list('gene'), 
               report=FALSE, includePop=TRUE) {
  if(is.null(epidemic)&is.null(data)) {
    warning('One of the parameters epidemic or data should be specified, but both were NULL. Returning NA.')
    NA
  } else if(!is.null(data)) {
    pop <- data
  } else {
    pop <- extractPop(epidemic=epidemic, 
                      sampledOnly=sampledOnly, activeOnly=activeOnly, 
                      tMin=tMin, tMax=tMax, firstN=firstN, lastN=lastN)
  }
  if(nrow(pop)>0) {
    if(!is.null(GEValues)) {
      pop[, z:=calcValue(env, gene, e, GEValues)]
    }
    
    # numbers of indivs in groups
    nums <- pop[, list(ni=length(z)), by=eval(parse(text=by))]
    
    # total number of individuals
    N <- nums[, sum(ni)]
    
    #number of groups
    K <- nrow(nums)
    
    # weighted mean number of individuals in each group
    n0 <- nums[, (N-sum(ni^2/N))/(K-1)]
    
    # grand mean
    zBar <- pop[, mean(z)]
    
    # Total sum of squares
    SSt <- pop[, sum((z-zBar)^2)]
    
    # gruop means assigned to each row to facilitate sums of squares
    pop[, zBari:=mean(z), by=eval(parse(text=by))]
    
    # sum of squares between groups
    SSb <- pop[, sum((zBari-zBar)^2)]
    
    # sum of squares within groups
    SSe <- pop[, sum((z-zBari)^2)]
    
    # mean square between
    MSb <- SSb/(K-1)
    # mean square within
    MSe <- SSe/(N-K)
    
    sigmaG2 <- (MSb-MSe)/n0
    sigmae2 <- MSe
    rA <- sigmaG2/(sigmaG2+sigmae2)
    
    # F-statistics and 95% CI
    F <- MSb/MSe
    Fu <- qf(0.975, df1=K-1, df2=N-K)
    Fl <- 1/qf(0.975, df1=N-K, df2=K-1)
    CI95upper <- (F/Fl-1)/(F/Fl+n0-1)
    CI95lower <- (F/Fu-1)/(F/Fu+n0-1)
    
    # Standard error: use an approximate formula from Lynch (eq. 18.21, p. 562, ch. 18)
    SE = sqrt(2*(1-rA)^2*(1+(n0-1)*rA)^2/(N*(n0-1)))
    
    if(report) {
      res <- list(rA=rA, H2aov=rA, N=N, K=K, nums=nums, 
                  SSt=SSt, SSb=SSb, SSe=SSe, MSb=MSb, MSe=MSe, n0=n0, SE = SE,
                  F=F, Fu=Fu, Fl=Fl, CI95upper=CI95upper, CI95lower=CI95lower, 
                  sigmaG2=sigmaG2, sigmae2=sigmae2)
      if(includePop) {
        res$pop <- pop
      }
      res
    } else {
      rA
    }
  } else {
    NA
  }
  
}


#' Broad-sense heritability of a pathogen trait
#' @param epidemic a list of objects returned from simulateEpidemic
#' @param data NULL or a data.table such as the element gen in an epidemic list
#' @param GEValues genotype-environment trait values for calculating the z-values; if NULL it is assumed that the 
#' z-values are already in data, or epidemic$gen.
#' @param sampledOnly,activeOnly,tMin,tMax,firstN,lastN parameters passed to extractPop
#' @return numeric indicating the estimated heritability for infected individuals in the population. 
#'  
#' @import data.table
#' 
#' @export
R2adj <- function(epidemic=NULL, data=NULL, GEValues=NULL, 
                  sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, 
                  firstN=Inf, lastN=Inf) {
  if(is.null(epidemic)&is.null(data)) {
    warning('One of the parameters epidemic or data should be specified, but both were NULL. Returning NA.')
    NA
  } else if(!is.null(data)) {
    pop <- data
  } else {
    pop <- extractPop(epidemic=epidemic, 
                      sampledOnly=sampledOnly, activeOnly=activeOnly, 
                      tMin=tMin, tMax=tMax, firstN=firstN, lastN=lastN)
  }
  if(nrow(pop)>0) {
    if(!is.null(GEValues)) {
      pop[, z:=calcValue(env, gene, e, GEValues)]
    }
    if(pop[, length(unique(gene))>=2]) {
      pop[, summary(lm(z~as.factor(gene)))$adj.r.squared]
    } else {
      0
    }
  } else {
    NA
  }
}


#' Bootstrap ANOVA 
#' @param i integer vector identifiers of the individuals
#' @param idPair integer vector showing the membership of individuals in pairs 
#' or bigger classes
#' @param z numeric vector of phenotypes
#' @param bootstraps integer the number of bootstrap samples to perform 
#' (default 1000)
#' @param ... other columns that are appended to the resulting data.table if
#' as.dt = TRUE. This parameter is supplied for convinience only.
#' @return list iwth estimated rA and 95% confidence intervals
#' @details the function performs 1000 bootstraps
#' 
#' @import data.table 
#' @import boot
#' 
#' @export
rAboot <- function(i, idPair, z, bootstraps=1000, as.dt=FALSE, ...) {
  stats <- c('rA', 'sigmaG2', 'sigmae2', 'F')
  if(length(z) < 2) {
    if(as.dt) {
      return(data.table(stat=stats,
                        N=length(z), 
                        K=length(z)/2,
                        n=2,
                        est=rep(NA, length(stats)),
                        bCI=lapply(1:length(stats), function(i) c(NA, NA)),
                        tips=list(i),
                        bSample=lapply(1:length(stats), function(i) numeric()),
                        ...))
    } else {
      return(list(tips=list(i), bootstrap=NULL, 
                  bCI95lower=NA,
                  bCI95upper=NA,
                  rA=NA, 
                  CI95lower=NA, CI95upper=NA, 
                  sigma2G=NA, sigma2z=NA, n=NA, 
                  N=length(z), K=0))
    }
  } else {
    data <- data.table(gene=idPair, z)
    
    aovReport=rA(data=data, report=TRUE, includePop=FALSE)
    
    if(bootstraps>0)  {
      bootstrap <- boot::boot(data=data[, unique(gene)], statistic=function(idPP, ids) {
        unlist(rA(data=data[gene%in%idPP[ids]], report=TRUE, 
                  includePop=FALSE)[stats])
      }, R=bootstraps)  
      
      bCI <- try(boot::boot.ci(bootstrap, type='basic'), silent=TRUE) 
      
      bCI95lower <- ifelse(class(bCI)!='try-error'&is.list(bCI), bCI$basic[4], NA)
      bCI95upper <- ifelse(class(bCI)!='try-error'&is.list(bCI), bCI$basic[5], NA)
    } else {
      bootstrap <- NULL
      bCI <- NULL
      bCI95lower <- bCI95upper <- NA
    }
    
    
    if(as.dt) {
      data.table(stat=stats, 
                 N=aovReport$N, 
                 K=aovReport$K,
                 n=aovReport$n0,
                 est=unlist(aovReport[stats]),
                 CI=list(c(aovReport$CI95lower, aovReport$CI95upper),
                         c(NA, NA), c(NA, NA), c(NA, NA)),
                 bCI=lapply(1:length(stats), function(i) {
                   if(is.list(bootstrap)) {
                     bCI <- try(boot::boot.ci(bootstrap, type='basic', index=i), silent=TRUE) 
                     bCI95lower <- ifelse(class(bCI)!='try-error' & is.list(bCI), 
                                          bCI$basic[4], NA)
                     bCI95upper <- ifelse(class(bCI)!='try-error'&is.list(bCI), 
                                          bCI$basic[5], NA)
                     c(bCI95lower, bCI95upper)  
                   } else {
                     c(NA, NA)
                   }
                 }),
                 tips=list(i),
                 bSample=lapply(1:length(stats), function(i) {
                   if(is.list(bootstrap)) {
                     bootstrap$t[,i]  
                   } else {
                     NULL
                   }
                 }),
                 ...)
    } else {
      list(tips=list(i), bootstrap=list(bootstrap),
           bCI95lower=bCI95lower,
           bCI95upper=bCI95upper,
           rA=aovReport$H2aov, 
           CI95lower=aovReport$CI95lower, 
           CI95upper=aovReport$CI95upper, 
           sigma2G=aovReport$sigma2G, 
           sigma2z=aovReport$sigma2G+aovReport$sigma2E, 
           n=aovReport$n0, 
           N=aovReport$N, K=aovReport$K)    
    }
  }
}

#' Perform phylogenetic pair analysis
#' @param exclude a data.table specifying selection filters to be applied after 
#' the extraction of phylogenetic pairs. Every row specifies one such filter in
#' the form of R-expressions to be evaluated within the data.table of 
#' extracted phylogenetic pairs (see exctractPP for the description of columns 
#' available in one such data.table). For each row in the pp-table, for which
#' an expression evaluates to TRUE, the corresponding phylogenetic pair rows 
#' (the row itself and its partner row) get removed from the ANOVA analysis. 
#' The exclude data.table has three character vector columns as follows: 
#' name - meaningful name of the filter used as an index (key) in the data.table,
#' scopeAll - R-expression  evaluated before grouping by quantiles of tau.
#' scopeTau - R-expression evaluated within each tau-quantile group. 
#' Note that the filters are applied on the pp data.table or after a call to 
#' extractPP, so they cannot affect the formation of phylogenetic pairs. Note 
#' also that if a member of a phylogenetic pair gets excluded its 
#' pair-partner is also removed even if the filter expression does not evaluate 
#' to TRUE for it. For example the filter 
#' filters=data.table(name=c('all', 'no outliers in deltaz'), 
#'                    scopeAll=('FALSE'), 
#'                    scopeTau=c('FALSE', 
#'                               'deltaz>\{q=quantile(deltaz); q[4]+1.5*(q[4]-q[2])\}'),
#'                    key='name') 
#' would result in a PP analysis on all phylogenetic pairs, and a PP analysis on
#' the phylogenetic pairs which's phenotypic distance deltaz doesn't exceed 
#' .75%+1.5*IQR of deltazs within each tauQuantile group. By default, no 
#' filtering is done.
#' @param ... currently not used
#' 
#' @return a data.table with the estimated statistics for each tauQuantile (and filter).
#' 
#' @export
PP <- function(z, tree=NULL, dists=NULL, pp=NULL, seed=NA, 
               zName='z', treeName='tree', distsName='dists',
               tauQuantiles=c(V=.05, D=.1, O=.125, qu=.2, Q=.25, M=.5, A=1),
               bootstraps=0, 
               exclude=data.table(name=c('all'), 
                                  scopeAll=c('FALSE'), 
                                  scopeTau=c('FALSE'), 
                                  key='name'),
               verbose=FALSE, ...) {
  
  if(is.list(z)) {
    p <- z
    z <- p[[zName]]
    if(is.null(z)) {
      z <- p[['v']]
    }
    tree <- p[[treeName]]
    dists <- p[[distsName]]
    
    if(is.null(z)|(is.null(tree)&is.null(dists)&is.null(pp))) {
      stop('If a list is supplied as argument z, this list should contain a 
           vector of trait values named "z" or zName and either a phylo-object 
           named "tree" or treeName or a dists matrix named "dists" or 
           distsName.')
    }
  }
  
  if(is.null(z)) {
    if(is.null(pp)) {
      stop('If z is NULL, pp should be supplied and it should contain a numerical column called "z".')
    } else {
      z <- pp[, z]
    }
  } 
  
  if(!is.null(dists)) {
    if(ncol(dists)!=nrow(dists)) {
      stop('dists is not a square matrix.')
    }
    N <- nrow(dists)
  } else if(!is.null(tree)) {
    N <- length(tree$tip.label)
  }
  
  if(any(is.na(z))) {
    stop('NAs or NaNs found in z.')
  }
  
  if(is.null(exclude) | !is.data.table(exclude)) {
    exclude <- data.table(name=c('all'), 
                          scopeAll=c('TRUE'), 
                          scopeTau=c('TRUE'), 
                          key='name')
    warning('Incorrectly specified filter; Using default filtering.')
  } else if(is.data.table(exclude)) {
    if(!setequal(names(exclude), c('name', 'scopeAll', 'scopeTau'))) {
      stop('The column names in exclude should be "name", "scopeAll", "scopeTau".')
    }
  }
  
  if(is.null(pp)) {
    pp <- extractPP(tree=tree, tipDists=dists)
    setkey(pp, i)
  }
  
  # If a vector z is supplied as parameter, ensure that
  # the PP-analysis is done on this vector and not on 
  # previously set column z in pp.
  if('z' %in% names(pp)) {
    pp[, z:=NULL]
  }
  
  pp[, z:=z[i]]
  pp[, deltaz:=abs(z[1]-z[2]), by=idPair]
  
  if(!is.na(seed)) {
    set.seed(seed)
  }
  
  res <- list()
  res$tree <- tree
  res$z <- z
  
  if(!is.null(tree)) {
    res$N <- length(tree$tip.label) 
  } else if(!is.null(dists)) {
    res$N <- nrow(dists)
  }
  
  stats <- rbindlist(
    lapply(names(tauQuantiles), function(Qname) {
      if(verbose){
        cat(Qname, '...\n')
      }
      probs=seq(0, 1, by=tauQuantiles[Qname])
      pp[, paste0(Qname, 'tau'):={
        quants=quantile(tau, probs=probs)
        if(length(unique(quants))==1) {
          quants[1] <- .999*quants[1]
        }
        lower=names(quants)[match(unique(quants), quants)]
        upper=names(quants)[length(quants)-match(unique(quants), rev(quants))+1]
        labels=paste0(Qname,c('[', rep('(', length(lower)-2)),
                      lower[1:(length(lower)-1)],',',upper[2:length(upper)],']')
        
        cut(tau, breaks=unique(quants), labels=labels, include.lowest=TRUE)
      }]
      pp[, tauQuantile:=eval(parse(text=paste0(Qname, 'tau')))]
      
      dt <- rbindlist(lapply(exclude[, name], function(flt) {
        if(verbose) {
          cat('Processing filter', flt, '...\n')
        }
        pp[{
          exclude1 <- eval(parse(text=exclude[list(flt), scopeAll]))
          excludeIds1 <- unique(c(i[exclude1], j[exclude1]))
          !(i%in%excludeIds1)
        }, 
        .SD[{
          exclude2 <- eval(parse(text=exclude[list(flt), scopeTau]))
          excludeIds2 <- unique(c(i[exclude2], j[exclude2]))
          !(i%in%excludeIds2)
        }, rAboot(i.name, idPair, z, bootstraps=bootstraps, as.dt=TRUE, 
                  filter=flt, 
                  tauMean=mean(tau), tauMedian=median(tau), 
                  tMean = mean(t), tMedian=median(t),
                  deltazMean=mean(deltaz), deltazMedian=median(deltaz))], 
        keyby=tauQuantile]
      }))
      pp[, tauQuantile:=NULL]
      if(verbose) {
        print(dt)
      }
      dt
    }))
  
  res$pp <- pp
  res$tauQuantiles <- tauQuantiles
  res$seed <- seed
  #res$ruleOutXIQR <- ruleOutXIQR
  res$exclude <- exclude
  res$bootstraps <- bootstraps
  
  res$stats <- stats
  
  class(res) <- c('PP', class(res))
  res
}

#' Summarize the results of a PP-call in a data.table
#' @return a data.table containing the rA, sigma2G, sigma2e and F-statistics for
#' each quantile-region of the phylogenetic distance tau analyzed in the PP-call.
#' 
#' @import data.table
#' 
#' @export
summary.PP <- function(object, tauQuantileTypes = NULL, ...) {
  stats <- copy(object$stats)
  stats[, tauQuantileType:=factor(
    regmatches(as.character(tauQuantile), 
               regexpr('^[^\\[\\(]+', as.character(tauQuantile), perl=TRUE)),
    levels=c('V', 'D', 'O', 'qu', 'Q', 'M', 'A'))]
  
  if(!is.null(tauQuantileTypes)) {
    stats <- stats[tauQuantileType %in% tauQuantileTypes]
  }
  
  if(is.null(stats[['CI']])) {
    # In previous versions, the CI was only evaluated using bootstrap
    # here, we fix such objects by calling with 0 bootstraps
    if(is.null(object[['tauQuantiles']])) {
      rep <- PP(z=NULL, pp=object$pp, bootstraps=0)  
    } else {
      rep <- PP(z=NULL, pp=object$pp, bootstraps=0, 
                tauQuantiles=object$tauQuantiles)
    }
    stats <- merge(stats, rep$stats[, list(tauQuantile, stat, outliers, CI)], 
                   by=c('tauQuantile', 'stat', 'outliers'))
  }
  #if(!is.null(stats[['CI']])) {
  stats[, CI.lower:=sapply(CI, function(.) .[1])]
  stats[, CI.upper:=sapply(CI, function(.) .[2])]
  stats[, CI:=NULL]
  
  stats[, bCI.lower:=sapply(bCI, function(.) .[1])]
  stats[, bCI.upper:=sapply(bCI, function(.) .[2])]
  stats[, bCI:=NULL]
  
  stats[, filter:=factor(filter)]
  
  stats[, tips:=NULL]
  stats[, bSample:=NULL]
  stats[, n:=NULL]
  
  class(stats) <- c('summary.PP', class(stats))
  stats
}


#' A basic ggplot of a PP analysis
#' @param object an object of class 'PP' (see ?summary.PP)
#' @param tauQuantileType a character vector containing the types of quantiles
#' for which plots should be returned. If NA, all available types are plotted.
#' Accpetable types are c('V', 'D', 'O', 'qu', 'Q', 'M', 'A'), corresponding to 
#' the following fractions: c(0.05, 0.1, 0.125, 0.2, 0.25, 0.5, 1).
#' 
#' @importFrom ggplot2 ggplot facet_grid facet_wrap aes geom_boxplot geom_point geom_segment
#' @export
plot.PP <- function(object, abs=TRUE,
                    tauQuantileType=c('V', 'D', 'O', 'qu', 'Q', 'M', 'A'),
                    mapping='aes(tau, deltaz)',
                    filter='TRUE',
                    facets=tauQuantileType~.,
                    ...) {
  pp <- copy(object$pp)
  pp[, zj:=z[match(j, i)]]
  
  if(!abs) {
    pp[, deltaz:=z-z[match(j, i)]]
  }
  pplong <- rbindlist(
    lapply(tauQuantileType, function(tqt) {
      tqttau <- paste0(tqt, 'tau')
      if(length(match(tqttau, names(pp)))==1) {
        colIds <- match(c('i', 'j', 'idPair', 'tau', 'z', 'zj', 'deltaz', 't', tqttau),
                        names(pp))
        res <- pp[, colIds, with=FALSE]
        setnames(res, tqttau, 'tauQuantile')
        res[, tauQuantileType:=factor(tqt, levels=tauQuantileType)]
        res[eval(parse(text=filter))]
      } else {
        NULL
      }
    })
  )
  res <- ggplot(data=pplong, eval(parse(text=mapping))) 
  if(!is.null(facets)) {
    print(class(facets))
    res <- res + facet_grid(facets)
  }
  
  res
}

#' A box-plot of phenotypic distances against phylogenetic distances from a PP
#' analysis
#' @param object An object of class PP.
#' @param abs Logical indicating whether to plot absolute phenotypic distances.
#' @param filter A character evaluating as an R expression. Defaults to TRUE.
#' @param mapping  A character string. Defaults to 
#' @param facets A formula passed to facet_grid 'aes(tau, deltaz)'
#' 
#' @export
boxplot.PP <- function(object, abs=TRUE,
                       tauQuantileType=c('V', 'D', 'O', 'qu', 'Q', 'M', 'A'),
                       filter='TRUE',
                       mapping='aes(tau, deltaz)',
                       facets=tauQuantileType~.,
                       ...) {
  pp <- copy(object$pp)
  pp[, zj:=z[match(j, i)]]
  if(!abs) {
    pp[, deltaz:=z-z[match(j, i)]]
  }
  pplong <- rbindlist(
    lapply(tauQuantileType, function(tqt) {
      tqttau <- paste0(tqt, 'tau')
      if(length(match(tqttau, names(pp)))==1) {
        colIds <- match(c('i', 'j', 'idPair', 'tau', 'z', 'deltaz', 't', tqttau),
                        names(pp))
        res <- pp[, colIds, with=FALSE]
        setnames(res, tqttau, 'tauQuantile')
        res[, tauQuantileType:=factor(tqt, levels=tauQuantileType)]
        res[eval(parse(text=filter))]
      } else {
        NULL
      }
    })
  )
  res <- ggplot(pplong, eval(parse(text=mapping))) + 
    geom_boxplot(aes(group=tauQuantile))
  
  if(!is.null(facets)) {
    res <- res + facet_grid(facets)
  }
  res
}

#' Plot the summary of a PP analysis
#' 
#' @export
plot.summary.PP <- function(object, statistic='rA', 
                            tauQuantileType=NA, 
                            outliers=NA,
                            filter='TRUE',
                            CI='bootstrap',
                            facets=tauQuantileType~.) {
  summ <- object[stat==statistic[[1]]]
  tqt <- tauQuantileType
  outl <- outliers
  if(!is.na(tqt)) {
    summ <- summ[as.character(tauQuantileType)%in%as.character(tqt)]
  }
  if(!is.na(outl)) {
    summ <- summ[outliers%in%outl]
  }
  summ <- summ[eval(parse(text=filter))]
  
  res <- ggplot(data=summ) + 
    geom_point(aes(x=tauMean, y=est))
  
  if(tolower(CI)=='bootstrap' & summ[, all(!is.na(bCI.lower))]) {
    res <- res + geom_segment(aes(x=tauMean, xend=tauMean, 
                                  y=bCI.lower, yend=bCI.upper))
  } else {
    if(tolower(CI)=='bootstrap') {
      warning('No bootstrap CI available, using theoretical CI.')
    }
    res <- res + geom_segment(aes(x=tauMean, xend=tauMean, 
                                  y=CI.lower, yend=CI.upper))
  }
  
  if(!is.null(facets)) {
    res <- res + facet_grid(facets)
  }
  
  res
}

#' Scatter plot of phylogenetic pairs
#' @param z numeric vector of phenotypes
#' @param tree a phylo object with tip-labels corresponding to the entries in z
#' @param ppAnalysis a list resulting from a call to analyseCPPs on the same tree and z.
#' @param CPPthr numeric indicating the maximum phylogenetic distance separating a 
#' closest phylogenetic pair
#' @param ruleOutXIQR numeric a multiplier used to define outlier CPPs as defined in 
#' the referenced article.
#' @param zName,treeName characters used when the parameter z is a list; indicate the
#' corresponding names in the list
#' @param xlim,ylim,xlab,ylab graphical parameters passed to plot
#' 
#' @import data.table 
#' 
#' @export
scatterPlotPPs <- function(z, tree, ppAnalysis, CPPthr=10^-4, ruleOutXIQR=1.5, 
                           zName='z', treeName='tree', 
                           xlab=expression(lg(tau)~"[lg(subst. per site)]"), 
                           ylab=expression(group("|",Delta~lg(spVL),"|")),  
                           xlim=c(-6,0), ylim=c(0,5)) {
  if(is.list(z)) {
    p <- z
    z <- p[[zName]]
    if(is.null(z)) {
      z <- p[['v']]
    }
    tree <- p[[treeName]]
    
    if(is.null(z)|is.null(tree)) {
      stop('If a list is supplied as argument z, this list should contain a vector of trait values named "z" or zName and a phylo-object named "tree" or treeName')
    }
  }
  N <- length(tree$tip.label)
  rootTipDists <- POUMM:::nodeTimes(tree)[1:N]
  pp <- ppAnalysis$pp
  pp[, dRoot:=rootTipDists[i]]
  zDists <- ppAnalysis$zDists
  tipDists <- ppAnalysis$tipDists
  
  s <- sample(x=1:length(zDists), size=nrow(pp), replace=FALSE)
  
  plot(x=log10(tipDists[s]), y=zDists[s], type='p', pch=20, cex=0.25, col=adjustcolor('darkgrey', alpha.f=0.3), 
       xlab=xlab, ylab=ylab, 
       xlim=xlim, ylim=ylim)
  pp[, points(x=log10(d), y=deltaz, pch=20, cex=0.25, col=adjustcolor('darkgreen', alpha.f=0.3))]
  pp[tau<=CPPthr, points(x=log10(d), y=deltaz, pch=20, cex=0.25, col=adjustcolor('magenta', alpha.f=0.3))]
  pp[tau<=CPPthr & 
       deltaz>{
         q=quantile(unique(deltaz[tau<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}, 
     points(x=log10(d), y=deltaz, col='blue', pch=20, cex=0.4)]
}


#' Box plots of trait values along a tree
#' @param nGroups integer the number of groups
#' @param ... additional parameters passed to boxplot
#' 
#' @import data.table
#' 
#' @export
boxplotTraitAlongTree <- function(z, tree, nGroups=15, ...) {
  groups <- groupByRootDist(tree, nGroups=nGroups)
  midDistPoints <- groups$groupMeans
  rootTipDistGroups <- groups$rootTipDistGroups
  boxes <- lapply(midDistPoints, function(.) c())
  
  for(i in 1:length(rootTipDistGroups)) {
    boxes[[as.character(rootTipDistGroups[[i]])]] <-
      c(boxes[[as.character(rootTipDistGroups[[i]])]], z[i])
  }
  
  boxes <- boxes[sort(names(boxes))]
  names(boxes) <- as.character(midDistPoints)
  
  boxes <- lapply(boxes, function(.) as.numeric(.))
  
  params <- list(...)
  myParams <- list(varwidth=TRUE, 
                   names=names(boxes), cex=.4, border='grey40',
                   at=as.numeric(names(boxes)), 
                   main='', 
                   xlab='root-tip distance [subst. per site]',
                   ylab=expression(lg(spVL)))
  
  params <- c(params, myParams[setdiff(names(myParams), names(params))])
  
  do.call(boxplot, c(list(boxes), params))
}
