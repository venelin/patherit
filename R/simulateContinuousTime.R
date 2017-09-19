library(ape)
library(data.table)

#' @title Simulate within and between host dynamics in an epidemic
#' 
#' @param Ninit integer, initial population size
#' @param nu numeric, birth rate (constant unit time influx in the population)
#' @param mu numeric, natural death rate
#' @param pe numeric vector of frequences for each environment type
#' @param sde numeric, the standard deviation of the special environmental effects
#' @param pg.init numeric initial frequencies from which the pathogen genotype of the first case is drawn
#' @param GEValues numeric matrix with average value for individuals of specific immune system type infected by a 
#' specific pathogen. GEValues[i, j] is for individaul of immune type i and pathogen genotype j.
#' @param rateContact numeric, the rate at which an infected individual has contact with a susceptible.
#' @param rateInfect a function of two arguments, a numeric vector of trait values 
#' and a numeric rateContact parameter, returning a numeric vector of transmission rates corresponding 
#' z. 
#' @param rateDie a function of one argument, a numeric vector of trait values, 
#' returning a numeric vector of death rates corresponding to these values.
#' @param rateSample numeric, the rate at which infected individuals get diagnozed and recover
#' @param rateMutate a function of four arguments: the GEValues matrix (see above), a numeric vector es of special
#' environmental effects, an integer vector envs of environment type indices and an integer vector genes of
#' pathogen genotype indices, returning a numeric vector of the same size as the es with per-locus mutation rates.
#' @param rateTransTemplate a square matrix of size n x (n-1) containing the factor by which the per locus 
#' mutation rate must be multiplied in order to obtain a transition rate between any two pathogen genotypes. See
#' the function generateRateTransTemplate. rateTransTemplate[i,j] should be the factor for transition rate from state i to state j. By
#' definition of transition rate matrices the transition from a state to itself is always the negative sum of the 
#' transition rates to the other states. This parameter needs not to be specified explicitly if the 
#' parameter numsAllelesAtSites is specified (see below). 
#' @param numsAllelesAtSites integer vector specifying the number of alleles at each site, 
#' the length of the vector defining the number of QTLs. 
#' @param eUniqForEachG logical indicating should the special environmental effect in an individual be drawn for
#' each possible pathogen genotype or should there be a single special environmental effect for all genotypes.
#' @param selectWithinHost logical indicating if there is selection for higher trait values within a host. When 
#' TRUE a mutation results in a transition from one genotype to another only if it results in a higher trait value
#' for a host.
#' @param timeStep numeric, the time resolution at which infection, mutation, sampling and death events are 
#' sampled.
#' @param maxTime numeric, the maximum duration of a simulation before starting the graceful fade-away phase.
#' @param maxNTips numeric, the maximum number of sampled individuals before starting the graceful fade-away phase.
#' @param expandTimeAfterMaxNTips numeric greater or equal than 1 specifying how lone should the growth of the epidemic 
#' continue after reaching maxNTips or running the epidemic for maxTime/expandTimeAfterMaxNTips.
#' @param reportInterval positive integer, time interval at which some statistics such as total infected and 
#' recovered counts and frequencies of genotypes in the population are stored.
#' 
#' @import ape
#' @import data.table
#' 
#' @return A list of objects 
#' @export
simulateEpidemic <- function(
  Ninit=Inf, nu=0, mu=1/850, pe, sde, pg.init, GEValues, 
  rateContact, rateInfect, rateDie, rateSample, rateMutate, 
  rateTransTemplate=generateRateTransTemplate(numsAllelesAtSites),
  numsAllelesAtSites=NULL, eUniqForEachG=FALSE, selectWithinHost=FALSE, 
  timeStep=.1, maxTime=1200, maxNTips=8000, expandTimeAfterMaxNTips=1.5,
  reportInterval=12,  
  ...) {
  
  nGtps <- length(pg.init)
  nEnvs <- length(pe)
  
  largs <- list(pe, sde, pg.init, GEValues, 
                rateContact, rateInfect, rateDie, rateSample, 
                rateMutate, rateTransTemplate, numsAllelesAtSites, 
                eUniqForEachG, selectWithinHost, timeStep)
  names(largs) <- c(
    'pe', 'sde', 'pg.init', 'GEValues', 'rateContact', 'rateInfect', 'rateDie', 
    'rateSample', 'rateMutate', 'rateTransTemplate', 'numsAllelesAtSites',
    'eUniqForEachG', ' selectWithinHost', 'timeStep')
  print('Simulating epidemic on:') 
  for(i in names(largs)) {
    cat(i, ':\n')
    print(largs[[i]])
  }
  
  counts <- NULL
  countNames <- c(
    'time', 'Total', 'N', 'X', 'Y', 'Z', 'nTips', 'D', 'meanActive', 'sdActive',
    'qActive0%', 'qActive2.5%', 'qActive25%', 'qActive50%', 'qActive75%', 
    'qActive97.5%', 'qActive100%', 
    'R2adj.A', 'rA.A', 'b.AtInfectionA','b.A', 
    'R2adj.S', 'rA.S', 'b.AtInfectionS','b.S',
    'R2adj', 'rA', 'b.AtInfection', 'b', 
    'nCouplesA', 'nCouplesS', 'nCouples')
  for(k in 1:nEnvs) {
    countNames <- c(countNames, paste0('Y', k, 1:nGtps))
  }
  for(k in 1:nEnvs) {
    countNames <- c(countNames, paste0('Z', k, 1:nGtps))
  }
  
  iter <- 0
  l <- list(N=Ninit, nTips=0)
  countActive <- 0
  fadingEpidemic <- FALSE
  expandingTime <- FALSE
  done <- FALSE
  gen <- NULL
  timeEndGrowth <- maxTime
  
  while(!done) {  
    if(is.null(gen)|countActive==0) {
      l <- stepContTimenxNorm(
        l$N, nu, mu, NULL, pe=pe, sde=sde, pg.init=pg.init, GEValues=GEValues, 
        rateContact=rateContact, rateInfect=rateInfect, rateDie=rateDie, 
        rateSample=rateSample, 
        rateMutate=rateMutate, rateTransTemplate=rateTransTemplate, 
        eUniqForEachG=eUniqForEachG, selectWithinHost=selectWithinHost, 
        timeStep=timeStep, t=iter, fadingEpidemic=fadingEpidemic)
      if(l$hasNewInfections) {
        gen <- rbind(l$gen, l$genNew)
      } else {
        gen <- l$gen
      }
    } else {
      l <- stepContTimenxNorm(
        l$N, nu, mu, gen, pe=pe, sde=sde, pg.init=pg.init, GEValues=GEValues, 
        rateContact=rateContact, rateInfect=rateInfect, rateDie=rateDie, rateSample=rateSample, 
        rateMutate=rateMutate, rateTransTemplate=rateTransTemplate,
        eUniqForEachG=eUniqForEachG, selectWithinHost=selectWithinHost,
        edge=l$edge, edge.length=l$edge.length, nTips=l$nTips, timeStep=timeStep, t=iter, 
        fadingEpidemic=fadingEpidemic)
      
      if(l$hasNewInfections) {
        gen <- rbind(gen, l$genNew)
      }
    }
    
    countActive <- gen[, sum(active==1)]
    if(fadingEpidemic & countActive==0) {
      done=TRUE
    }
    if(!fadingEpidemic & !expandingTime & 
       (l$nTips>=maxNTips | iter*timeStep>maxTime/expandTimeAfterMaxNTips)) {
      timeEndGrowth <- iter*timeStep*expandTimeAfterMaxNTips
      expandingTime <- TRUE
      cat('Reached maxNTips, continuing growth unitl time', timeEndGrowth, '\n')
    }
    if(!fadingEpidemic & expandingTime & iter*timeStep>timeEndGrowth) {
      cat('Reached epidemic growth stop. Starting fading phase...\n')
      fadingEpidemic <- TRUE
      expandingTime <- FALSE
      iterActive <- iter
    }
    iter <- iter+1
    
    if((iter*timeStep)%%1==0) {
      epid <- list(gen=gen, timeStep=timeStep)
      cnt <- c(iter*timeStep, nrow(gen), 
               l$N, l$N- gen[, sum(alive==1)], gen[, sum(active==1)], 
               gen[, sum(alive==1&active==0)], l$nTips, 
               gen[, sum(alive==0)], 
               mean(l$zActive, na.rm=TRUE), sd(l$zActive, na.rm=TRUE), 
               quantile(l$zActive, probs=c(0, .025, .25, .5, .75, .975, 1), na.rm=TRUE), 
               
               R2adj(epidemic=epid, GEValues=GEValues, sampledOnly=FALSE, activeOnly=TRUE), 
               rA(epidemic=epid, data=NULL, GEValues=GEValues, sampledOnly=FALSE, activeOnly=TRUE), 
               b(epidemic=epid, GEValues=GEValues, sampledOnly=FALSE, activeOnly=TRUE, atInfection=TRUE), 
               b(epidemic=epid, GEValues=GEValues, sampledOnly=FALSE, activeOnly=TRUE, atInfection=FALSE),
               
               R2adj(epidemic=epid, GEValues=GEValues, sampledOnly=TRUE, activeOnly=FALSE), 
               rA(epidemic=epid, data=NULL, GEValues=GEValues, sampledOnly=TRUE, activeOnly=FALSE), 
               b(epidemic=epid, GEValues=GEValues, sampledOnly=TRUE, activeOnly=FALSE, atInfection=TRUE), 
               b(epidemic=epid, GEValues=GEValues, sampledOnly=TRUE, activeOnly=FALSE, atInfection=FALSE),
               
               R2adj(epidemic=epid, GEValues=GEValues, sampledOnly=FALSE, activeOnly=FALSE), 
               rA(epidemic=epid, data=NULL, GEValues=GEValues, sampledOnly=FALSE, activeOnly=FALSE), 
               b(epidemic=epid, GEValues=GEValues, sampledOnly=FALSE, activeOnly=FALSE, atInfection=TRUE), 
               b(epidemic=epid, GEValues=GEValues, sampledOnly=FALSE, activeOnly=FALSE, atInfection=FALSE),
               
               nrow(extractDRCouples(epidemic=epid, sampledOnly=FALSE, activeOnly=TRUE)), 
               nrow(extractDRCouples(epidemic=epid, sampledOnly=TRUE, activeOnly=FALSE)), 
               nrow(extractDRCouples(epidemic=epid, sampledOnly=FALSE, activeOnly=FALSE))
               )
      
      cnt <- c(cnt, gen[active == 1, {
        accum <- matrix(as.integer(NA), nGtps, nEnvs)
        for(k in 1:nEnvs) {
          for(j in 1:nGtps) {
            accum[j, k] <- sum(env == k & gene == j)
          }
        }
        accum
      }])
      
      cnt <- c(cnt, gen[sampled == 1, {
        accum <- matrix(as.integer(NA), nGtps, nEnvs)
        for(k in 1:nEnvs) {
          for(j in 1:nGtps) {
            accum[j, k] <- sum(env == k & gene == j)
          }
        }
        accum
      }])
      
      if(is.null(counts)) {
        counts <- cnt
      } else {
        counts <- rbind(counts, cnt)
      }
    
      if((iter*timeStep)%%reportInterval==0) {
        cat('### Time: ', iter*timeStep, '\n')
        names(cnt) <- countNames
        print(round(c(
          cnt[1:7], cnt[c('R2adj.A', 'R2adj.S', 'rA.A', 'rA.S')],
          gtpsDiv = sum(cnt[-(1:32)][1:((length(cnt)-32)/2)]!=0) /
            ((length(cnt)-32)/2)), digits=2))
      }
    }
  }
  
  colnames(counts) <- countNames
  
  c(l[c('N', 'edge', 'edge.length', 'nTips')], 
    list(gen=gen, time=timeStep*iterActive, timeFading=timeStep*iter, 
         timeStep=timeStep, rateContact=rateContact, 
         rateInfect=rateInfect, rateDie=rateDie, rateSample=rateSample, 
         rateMutate=rateMutate, eUniqForEachG=eUniqForEachG, 
         selectWithinHost=selectWithinHost, counts=counts,...))
}

#' norm-reaction plot
#' @export
normReactPlot <- function(genotypeXenv, GE, names=NULL, do.plot=T, ...) {
  normReactMat <- 
    rbind(as.data.table(cbind(genotypeXenv, GE))[ei==1, GE, keyby=gi][, GE],
          as.data.table(cbind(genotypeXenv, GE))[ei==2, GE, keyby=gi][, GE])
  
  if(do.plot) {
    matplot(x=1:nrow(normReactMat), y=normReactMat, type='p', 
            col=c(2:6,'brown'), pch=20, xaxt='n', 
            ...)
    if(!is.null(names)) {
      text(x=1, y=normReactMat[1,], labels=names, col=c(2:6, 'brown'), pos=4)
    }
    axis(1, at=1:nrow(normReactMat))
    matlines(x=(1:nrow(normReactMat))+c(.32, rep(0, nrow(normReactMat)-2), -0), 
             y=normReactMat, col=c(2:6, 'brown'), lty=1, lwd=1) 
  }
  normReactMat
}

#' @export
normReactMat <- function(genotypeXenv, GE) {
  normReactPlot(genotypeXenv, GE, do.plot=F)
}


#' Extract the transmission tree connecting sampled individuals in a simulated epidemic
#' @param epidemic list returned by simulateEpidemic
#' @param tips,idTips integer vectors (see details)
#' @param collapse.singles logical indicating whether non-bifurcating internal nodes should be collapsed using the ape function collapse.singles(). Default: TRUE
#' @details By default all individuals sampled during the active phase of the epidemic are 
#' included in the phylogeny. If specified, the argument tips has priority to the arugment 
#' idTips. The argument tips is sorted in increasing order and specifies the tip indices 
#' (numbers in the range 1:epidemic$nTips) to be retained in the phylogeny. If NULL (default), 
#' the argument idTips is used to identify which individuals to be included in the phylogeny. 
#' This is useful if, for example, one has searched all individuals who were sampled during 
#' a specified time interval of the epidemic.
#' 
#' @import ape
#' 
#' @export
extractTree <- function(epidemic, tips=NULL, idTips=NULL, collapse.singles=TRUE) {
  gen <- copy(epidemic$gen)
  edge.length <- epidemic$edge.length
  
  if(is.null(idTips)&is.null(tips)) {
    tips <- 1:epidemic$nTips
    edge <- epidemic$edge
    nTips <- epidemic$nTips 
  } else {
    if(is.null(tips)) {
      # get tips from idTips
      tips <- sort(gen[id%in%idTips, tip])
    }  
    nTips <- length(tips)
    edge <- epidemic$edge
 
    edge[, 3] <- 0
    etips <- match(-tips, edge[, 2])
    edge[etips, 3] <- 1
    nodeP <- edge[etips, 1]
    while(!(length(nodeP)==1&nodeP[1]==0) & any(edge[nodeP, 3]==0)) {
      edge[nodeP, 3] <- 1
      nodeP <- edge[nodeP, 1]
    } 
  }
  
  if(nTips <= 1)
    return(NULL)
  
  # the 3rd column in edge matrix is 1 if the edge participates on the route from a tip to the root.
  edge.length <- edge.length[edge[,3]==1]
  edge <- edge[edge[,3]==1, 1:2]
  
  cat('Generating tree: nTips=',nTips, ', number of uncollapsed edges=', nrow(edge), '\n')
  
  nTipsAll <- epidemic$nTips
  edge <- apply(edge, 1:2, function(x) if(x<0) -x else x+nTipsAll+1)
  nodesUnique <- sort(unique(as.vector(edge)))
  
  Nnode <- length(nodesUnique)-nTips
  
  edge[, 1] <- match(edge[, 1], nodesUnique)
  edge[, 2] <- match(edge[, 2], nodesUnique)
  
  tipl <- gen[tip%in%tips, list(id, tip)]
  setkey(tipl, tip)
  colnames(edge) <- NULL
  obj <- list(edge=edge, edge.length=edge.length, tip.label=tipl[,id], node.label=nodesUnique[-(1:nTips)]-nTipsAll-1, Nnode=Nnode)
  class(obj) <- 'phylo'
  if(collapse.singles) {
    obj <- collapse.singles(obj)
  }
  obj$edge.length <- obj$edge.length*epidemic$timeStep
  obj
}

#' Extract a sub-population of individuals involved in an epidemic
#' @param epidemic a list representing an epidemic (as returned by simulateEpidemic)
#' @param ids integer vector or NULL, specifying the id of individuals in epidemic$gen to be extracted
#' @param sampledOnly a logical indicating whether only recovered individuals should be extracted
#' @param activeOnly a logical indicating whether only currently infected and alive individuals should be extracted
#' @param tMin,tMax a numericals indicating the interval of observation times for extracted individuals. 
#' If sampledOnly is TRUE the observation time is the time of sampling (recovery), otherwise it is the time 
#' of infection. Default: tMin=0, tMax=Inf.
#' @param firstN integer indicating whether to return only the first firstN observed individuals (default is Inf)
#' @param lastN integer indicating whether to return only the most recently observed lastN individuals (default is Inf).
#' This filter is applied after the previous filters.
#' 
#' @return a data.table representing a subset of epidemic$gen after applying the filters
#' 
#' @import data.table
#' 
#' @export
extractPop <- function(epidemic, ids=NULL, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, firstN=Inf, lastN=Inf) {
  pop <- if(is.null(epidemic)) {
    warning('Parameter epidemic is NULL. Returning NULL.')
    NULL
  } else {
    epidemic$gen[{
      timeObs <- if(sampledOnly) {
        (tinf+tau)
      } else {
        tinf
      } 
      
      if(!is.null(ids)) {
        subset <- id%in%ids
      } else {
        subset <- rep(TRUE, length(id))
      }
      
      subset <- (subset & sampled >= ifelse(sampledOnly, 1, 0) & sampled<=1 & active>=ifelse(activeOnly, 1, 0) &
        timeObs*epidemic$timeStep >= tMin & timeObs*epidemic$timeStep <= tMax)
      
      N <- sum(subset)
      if(firstN < N) {
        timeFirstN <- sort(timeObs[subset])[firstN]
        subset <- subset & timeObs<=timeFirstN
      } 
      
      N <- sum(subset)
      if(lastN < N) {
        timeLastN <- sort(timeObs[subset])[N-lastN]
        subset <- subset & timeObs>=timeLastN
      }
      subset
    }]
  }
}

#' Extract donor-recipient couples from an epidemic
#' @param epidemic a list representing an epidemic (as returned by simulateEpidemic)
#' @param ids integer vector or NULL, specifying the id of individuals in epidemic$gen to be extracted
#' @param sampledOnly a logical indicating whether only recovered individuals should be extracted
#' @param activeOnly a logical indicating whether only currently infected and alive individuals should be extracted
#' @param tMin,tMax a numericals indicating the interval of observation times for extracted individuals. 
#' If sampledOnly is TRUE the observation time is the time of sampling (recovery), otherwise it is the time 
#' of infection. Default: tMin=0, tMax=Inf.
#' @param firstN integer indicating whether to return only the first firstN observed individuals (default is Inf)
#' @param lastN integer indicating whether to return only the most recently observed lastN individuals (default is Inf).
#' This filter is applied after the previous filters.
#'
#' @return a data.table with rows corresponding to the extracted couples with the following columns:
#'
#' idD: id of the donor
#' id: id of the recipient
#' envd: environment type of the donor
#' gd: transmitted strain from the donor to the recipient
#' ed: donor special environmental effect at the moment of transmission
#' env: environment type of the recipient
#' gene: strain in the recipient at the moment of sampling or at the current moment
#' e: special environmental effects of the recipient
#' tauR: age of infection in the recipient (in timeStep units)
#' tauDAtInf: age of infection in the donor at the moment of transmission (in timeStep units)
#' tauD: time (in timeStep units) in the donor from the moment of transmission until the moment 
#' of sampling or the current moment 
#' taum: =tauD+tauR
#' eD: special environmental effect in the donor at the moment of sampling
#' gD: strain in the donor at the moment of sampling
#' 
#' @import data.table
#' 
#' @export
extractDRCouples <- function(epidemic, ids=NULL, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, firstN=Inf, lastN=Inf) {
  if(is.null(epidemic)) {
    warning('Parameter epidemic is NULL. Returning NULL.')
    NULL
  } else {
    couples <- epidemic$gen[{
      timeObs <- if(sampledOnly) {
        (tinf+tau)
      } else {
        tinf
      } 
      
      if(!is.null(ids)) {
        subset <- id%in%ids
      } else {
        subset <- rep(TRUE, length(id))
      }
      
      subset <- (subset & sampled >= ifelse(sampledOnly, 1, 0) & sampled<=1 & 
                   active>=ifelse(activeOnly, 1, 0) &
                   timeObs*epidemic$timeStep >= tMin & 
                   timeObs*epidemic$timeStep <= tMax)
      
      N <- sum(subset)
      if(firstN < N) {
        timeFirstN <- sort(timeObs[subset])[firstN]
        subset <- subset & timeObs<=timeFirstN
      } 
      
      N <- sum(subset)
      if(lastN < N) {
        timeLastN <- sort(timeObs[subset])[N-lastN]
        subset <- subset & timeObs>=timeLastN
      }
      
      subset <- subset & (idD%in%id[subset])
      
      subset
    }, list(idD, id, envd, gd, ed, env, gene, e, tauR=tau, tauDAtInf=taud)]
    
    setkey(epidemic$gen, id)

    couples[, tauD:=epidemic$gen[J(couples[, idD]), tau]-tauDAtInf]
    couples[, taum:=tauD+tauR]
    
    couples[, gD0:=gd]
    couples[, gR0:=gd]

    if(nrow(couples)>0) {
      couples[, eR0:=eSpec(e, gd)]
    } else {
      couples[, eR0:=c()]
    }
    couples[, eD0:=ed]
    
    couples[, gD:=epidemic$gen[J(couples[, idD]), gene]]
    couples[, gR:=gene]
    if(nrow(couples)>0) {
      couples[, eD:=epidemic$gen[J(couples[, idD]), eSpec(e, gene)]]
      couples[, eR:=eSpec(e, gR)]
    } else {
      couples[, eD:=c()]
      couples[, eR0:=c()]
    }
      
    
    couples
  }
}

#' Decompose trait values according to the model z=G+I+E+epsilon
#' @param data a data.table returned by extractPop
#' @param GEValues a matrix of GEValues
#' @param copy logical indicating whether a copy of the original data.table should be returned (TRUE by default)
#' 
#' @import data.table
#' 
#' @export
decomposeTrait <- function(data, GEValues, copy=TRUE) {
  if(copy) {
    data <- copy(data)
  } else {
    data <- data
  }
    
  data[, z:=calcValue(env, gene, e, GEValues)]
  
  mu <- data[, mean(z)]
  
  # calculate genotypic values (grouping by genotype)
  data[, G:=mean(z), by=gene]
  
  # calculate environmental values (grouping by env-type)
  data[, E:=mean(z)-mu, by=env]
  
  # genotype by environment interaction
  data[, I:=mean(z-G-E), by=list(env, gene)]
  
  # special environmental effects
  data[, epsilon:=z-G-E-I]
  
  data
}


#' Variance decomposition of trait values 
#' 
#' @import data.table
#'@export
decomposeVar <- function(data) {
    data[, c(varz=var(z), varG=var(G), varE=var(E), varI=var(I), varEpsilon=var(epsilon), covGE=cov(G,E))]
}


#' @import data.table
#' @export
estimMean <- function(epidemic, GEValues, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, lastN=Inf) {
  if(is.null(epidemic)) {
    warning('Parameter epidemic is NULL. Returning NA.')
    NA
  } else {
    gen <- copy(epidemic$gen[sampled >= ifelse(sampledOnly, 1, 0) & sampled<=1 & active>=ifelse(activeOnly, 1, 0), ])
    
    if(sampledOnly) {
      gen[, timeObs:=(tinf+tau)*epidemic$timeStep]
    } else {
      gen[, timeObs:=tinf*epidemic$timeStep]
    }
    gen <- gen[timeObs>=tMin&timeObs<=tMax]
    if(lastN<nrow(gen)) {
      gen <- gen[(nrow(gen)-lastN+1):nrow(gen)]
    }
    if(nrow(gen)>0) {
      gen[, z:=calcValue(env, gene, e, GEValues)]
      gen[, mean(z, na.rm=TRUE)]
    }
  }
}

#' @import data.table
#' @export
estimSD <- function(epidemic, GEValues, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, lastN=Inf) {
  if(is.null(epidemic)) {
    warning('Parameter epidemic is NULL. Returning NA.')
    NA
  } else {
    gen <- copy(epidemic$gen[sampled >= ifelse(sampledOnly, 1, 0) & sampled<=1 & active>=ifelse(activeOnly, 1, 0), ])
    
    if(sampledOnly) {
      gen[, timeObs:=(tinf+tau)*epidemic$timeStep]
    } else {
      gen[, timeObs:=tinf*epidemic$timeStep]
    }
    gen <- gen[timeObs>=tMin&timeObs<=tMax]
    if(lastN<nrow(gen)) {
      gen <- gen[(nrow(gen)-lastN+1):nrow(gen)]
    }
    if(nrow(gen)>0) {
      gen[, z:=calcValue(env, gene, e, GEValues)]
      gen[, sd(z, na.rm=TRUE)]
    }
  }
}

#' @import data.table
#' @export
estimQuantile <- function(epidemic, GEValues, prob=.5, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, lastN=Inf) {
  if(is.null(epidemic)) {
    warning('Parameter epidemic is NULL. Returning NA.')
    NA
  } else {
    gen <- copy(epidemic$gen[sampled >= ifelse(sampledOnly, 1, 0) & sampled<=1 & active>=ifelse(activeOnly, 1, 0), ])
    
    if(sampledOnly) {
      gen[, timeObs:=(tinf+tau)*epidemic$timeStep]
    } else {
      gen[, timeObs:=tinf*epidemic$timeStep]
    }
    gen <- gen[timeObs>=tMin&timeObs<=tMax]
    if(lastN<nrow(gen)) {
      gen <- gen[(nrow(gen)-lastN+1):nrow(gen)]
    }
    if(nrow(gen)>0) {
      gen[, z:=calcValue(env, gene, e, GEValues)]
      gen[, quantile(z, na.rm=TRUE)]
    }
  }
}

# simulating evolution on continuous time scale 

# sde : numeric or vector of standard deviations of the environmental effects for each allele
# N : number of individuals
# nGtps : number of alleles 
# sde : a numeric or a numeric matrix, sde[env, g]: standard dev of the special env. effect for <env, g>
#' @import data.table
newgennxNorm <- function(pe, sde, pg.init, N, eUniqForEachG=FALSE) {
  nGtps <- length(pg.init)
  if(length(sde)==1) {
    sde <- matrix(sde, length(pe), nGtps)
  }
  env <- sample(1:length(pe), N, replace=TRUE, prob=pe)
  if(eUniqForEachG) {
    cbind(env, gene=sample(1:nGtps, N, replace=TRUE, prob=pg.init), 
          matrix(sapply(1:nGtps, function(g) rnorm(N, 0, sde[cbind(env, g)])), N))
  } else {
    cbind(env, gene=sample(1:nGtps, N, replace=TRUE, prob=pg.init),
          matrix(rep(rnorm(N, 0, sde[cbind(env, 1)]), nGtps), N))
  }
}

#' @import data.table
#' @useDynLib patherit
sampleEvent <- function(rates, timeStep) {
  # total rate at which an event occurs during this timeStep
  lambdas <- rowSums(rates)
  
  # probabilities that no event has happened during the timeStep
  probs0 <- exp(-lambdas*timeStep) 
 
  # set intervals of propbabilities for no-event, sampling, dieing, transmitting and mutating to each other strain
  probs <- cbind(probs0, rates/lambdas*(1-probs0))
  
  # for(i in 2:ncol(probs)) {
  #  probs[, i] <- probs[, i-1]+probs[, i]
  # }
  # a faster C++ version of the above:
  cumulativeRows(probs)
  
  random <- runif(nrow(rates))
  probs_random <- (probs-random) >= 0
  
  events <- rep(0, nrow(rates))
  # for(e in ncol(probs_random):1) {
  #   events[probs_random[, e]] <- e-1
  # }
  # a faster C++ version of the above:
  decideEvents(events, probs_random)
  
  events
}

#' @import data.table
infectnxNorm <- function(gen, pe, sde, pg.init, donors, eUniqForEachG=FALSE) {
  nGtps <- length(pg.init)
  nNew <- sum(donors)
  if(nNew) {
    gen.r <- newgennxNorm(pe, sde, pg.init, sum(donors), eUniqForEachG)
    gen.r[, 'gene'] <- gen[donors, gene]
    data.table(env=as.integer(gen.r[, 'env']), gene=as.integer(gen.r[, 'gene']), 
               e=lapply(seq_len(nrow(gen.r)), function(i) gen.r[i, 2+(1:nGtps)]), idD=gen[donors, id])
  } else {
    NULL
  }
}

#' @import data.table
#' @export
eSpec <- function(e, gene) {
  if(is.list(e)) {
    eMat <- do.call(rbind, e)
    es <- rep(NA, length(gene))
    for(g in seq_len(ncol(eMat))) {
      es[gene==g] <- eMat[gene==g, g]
    }
    es
  } else if(is.matrix(e)) {
    es <- rep(NA, length(gene))
    for(g in seq_len(ncol(e))) {
      es[gene==g] <- e[gene==g, g]
    }
    es
  } else {
    e
  }
}

#' calculate phenotypic values
#' @param env integer vector containing environmental types for each individual
#' @param gene integer vectory containing genotype indices for each individual
#' @param e numeric vector containing the special environmental effects for each individual
#' @param GEValues matrix containing the genotype by environment expected phenotypic values
#' @return  vector of phenotypic values
#' @import data.table
#' @export
calcValue <- function(env, gene, e, GEValues) {
  GEValues[cbind(env, gene)] + eSpec(e, gene)
}

#' Possible genotypes based on numbers of alleles at QTLs.
#' 
#' @param numsAllelesAtSites integer vector specifying the number of alleles at each site, 
#' the length of the vector defining the number of QTLs.
#' @import data.table
#' @export
generateGenotypes <- function(numsAllelesAtSites) {
  if(length(numsAllelesAtSites)==0) {
    stop("numsAllelesAtSites should be an integer vector of numbers of alleles bigger than 1.")
  } 
  if(any(numsAllelesAtSites<=1)) {
    stop("numsAllelesAtSites should be an integer vector of numbers of alleles bigger than 1.")
  }
  
  genotypes <- as.matrix(expand.grid(lapply(rev(numsAllelesAtSites), 
                                            function(all) 1:all))[, length(numsAllelesAtSites):1])
  
  
  # a matrix of genotypes: each row is a genotype containing the names of the alleles at each site
  genotypes <- apply(genotypes, 1, function(r) as.vector(r))
  
  if(is.matrix(genotypes)) {
    genotypes <- t(genotypes)
  } else {
    genotypes <- as.matrix(genotypes)
  }
  
}

#' Generate a state transition template matrix based on numbers of alleles at QTLs.
#' 
#' @param numsAllelesAtSites integer vector specifying the number of alleles at each site, 
#' the length of the vector defining the number of QTLs.
#' @param nameRows logical indicating whether the returned matrix should have named rows
#' (for explicative purpose only, default FALSE). The column names are different for each row, 
#' i.e. the row-names excluding the genotype corresponding to that specific row.
#' @details each element of the matrix is a multiplier for the corresponding genotype transition rate 
#' (see function simulateEpidemic in the patherit package). This multiplier is nonzero only for transitions
#' that need exactly one mutation (haming distance = 1) in order to convert one genotype into another.
#' The 
#' 
#' @import data.table
#' 
#' @export
generateRateTransTemplate <- function(numsAllelesAtSites, nameRows=FALSE) {
  genotypes <- generateGenotypes(numsAllelesAtSites)
  #create the transition template matrix
  transTemplate <- matrix(0, nrow=nrow(genotypes), ncol=nrow(genotypes)-1)

  for(i in 1:nrow(genotypes)) 
    for(j in (1:nrow(genotypes))[-i]) {
      k <- j-as.integer(j>=i)
      hammingDist <- sum(genotypes[i, ]!=genotypes[j, ])
      if(hammingDist==1) {
        transTemplate[i, k] <- 1/(numsAllelesAtSites[genotypes[i, ]!=genotypes[j, ]]-1)
      }
    }
  if(nameRows) {
    rownames(transTemplate) <- apply(genotypes, 1, function(g) do.call(paste0, as.list(g)))  
  }
  
  transTemplate
}

#' @import data.table
# state-transition rates
rateStateTransition <- function(
  z, rates, eMatrix, es, envs, genes, GEValues, rateTransTemplate, selectWithinHost) {
  if(selectWithinHost) {
    gtps <- seq_len(ncol(eMatrix))
    gtpsMatrix <- matrix(gtps, nrow=nrow(eMatrix), ncol=ncol(eMatrix), byrow=TRUE)
    
    gother <- matrix(NA, nrow=nrow(eMatrix), ncol=ncol(eMatrix)-1)
    envsother <- matrix(NA, nrow=nrow(eMatrix), ncol=ncol(eMatrix)-1)
    for(g in gtps) {
      geneg <- genes==g
      gother[geneg, ] <- gtpsMatrix[geneg, -g]
      envsother[geneg, ] <- eMatrix[geneg, -g]
    }
    
    # is the mutation to other beneficial?
    sgnother <- sapply(1:(ncol(gother)), function(i) {
      as.double(GEValues[cbind(envs, gother[, i])] + envsother[, i] > z)  
    })
    
    matrix(sgnother * rateTransTemplate[genes, ] * rates, nrow = length(z))
  } else {
    matrix(rateTransTemplate[genes, ] * rates, nrow = length(z))
  }
}

# X: current number of susceptible individuals
# nu: birth rate
# mu: natural per capita death rate
# gen: sub-population of infected and sampled (including those who already died)
# pe: environmental type probabilities
# sde: standard deviation of the special environmental effects
# pg.init: genotype probabilities for the first case
# nGtps: number of genotypes, should be the same 
# GE: genotypic x environmental values
# rateContact
# rateInfect
# rateDie
# rateSample
# rateMutate
# rateTransTemplate
# edge
# edge.length
# nTips
# timeStep
# t
# eUniqForEachG
# selectWithinHost
# fadingEpidemic
#' @import data.table
stepContTimenxNorm <- function(N=Inf, nu=0, mu,
                               gen=NULL, pe, sde, pg.init, GEValues, rateContact, rateInfect, rateDie, rateSample, 
                               rateMutate, rateTransTemplate, 
                               eUniqForEachG=FALSE, selectWithinHost=FALSE, 
                               edge=NULL, edge.length=c(), nTips=0, timeStep=1, t=0,
                               fadingEpidemic=FALSE) {
  nGtps <- length(pg.init)
  newInfections <- FALSE
  hasGen <- FALSE
  if(is.null(gen)) {
    # currently the first infected individual and its donor are always of immune system type 1
    gen <- newgennxNorm(c(1, 0), sde, pg.init, 1, eUniqForEachG)
    if(length(sde==1)) {
      envd <- as.integer(1)
      ed <- rnorm(1, 0, sde)
      zd=GEValues[envd, gen[,'gene']]+ed
    } else {
      envd <- as.integer(1)
      ed <- rnorm(1, 0, sde[1, gen[, 'gene']])
      zd=GEValues[envd, gen[,'gene']]+ed
    }
    
    gen <- do.call(data.table, 
                   c(list(env=as.integer(gen[, 'env']), gene=as.integer(gen[, 'gene']), 
                     e=lapply(seq_len(nrow(gen)), function(i) gen[i, 2+(1:nGtps)]), 
                     id=as.integer(1), idD=as.integer(0), 
                     tinf=t, tau=0, nodeP=as.integer(0), tauP=0, tip=0, 
                     alive=as.integer(1), active=as.integer(1), sampled=as.integer(0), 
                     envd=envd, gd=as.integer(gen[, 'gene']), ed=ed, 
                     zd=zd, taud=as.double(NA)), nrecips=0, 
                     eventTime=as.character(t), eventCode='0')) 
    
    nTips=0
    edge <- matrix(0, nrow=0, ncol=3)
    edge.length <- c()
    hasGen <- TRUE
  } 
  
  # number of susceptible individuals
  X <- N - gen[, sum(alive==1)]
  
  # let some of the alive and recovered (sampled) individuals die of natural death
  # number of recovered individuals still alive
  Z <- gen[, sum(alive==1&sampled>=1)]
  Zdie <- rep(FALSE, nrow(gen))
  Zdie[gen[, alive==1&sampled>=1]] <- as.logical(rbinom(Z, 1, 1-exp(-mu*timeStep)))
  gen[Zdie, eventTime:=paste0(eventTime, ',', t)]
  gen[Zdie, eventCode:=paste0(eventCode, ',', 2)]
  gen[Zdie, alive:=as.integer(0)]
  
  # The rate of risky contacts is the rate of contacts times the 
  # frequency of susceptibles in the active population
  if(is.infinite(N)) {
    rateRiskyContact <- rateContact
  } else {
    rateRiskyContact <- rateContact*X/N
  }
  
  # will use this a lot, so cache it.
  isActive <- gen[, active==1]
  
  #increment times for all active
  # tau: age of infection
  # tauP: time since previous transmission
  gen[isActive, tau:=tau+1]
  gen[isActive, tauP:=tauP+1]
  
  # sample events that happen at this timeStep for some of the actives
  envs <- gen[isActive, env] # host immune system types
  genes <- gen[isActive, gene] # viral genotypes
  
  # host-specific effects for each viral genotype
  eMatrix <- gen[isActive, do.call(rbind, e)] 
  # host-specific effects for current viral genotypes
  es <- gen[isActive, eSpec(eMatrix, gene)]
  
  # current phenotypic values
  zs <- GEValues[cbind(envs, genes)]+es
  
  ratesSample <- rep(rateSample, length(zs))
  ratesDie <-  rateDie(zs, mu)        
  ratesInfect <- rateInfect(zs, rateRiskyContact)   
  
  ratesMutate <- rateMutate(GEValues, es, envs, genes)
  
  ratesTrans <- rateStateTransition(zs, ratesMutate, eMatrix, es, envs, genes, 
                                    GEValues, rateTransTemplate, selectWithinHost)
  
  rates <- cbind(ratesSample, ratesDie, ratesInfect, ratesTrans)
  
  events <- sampleEvent(rates, timeStep)
  # event codes are:
  # 0: no event; 1: sampling; 2: death; 3: infect other; 4, 5, 6,... mutation
  
  #cat(sum(events>=4), '/', length(events), '\n')
  
  # update N
  if(!is.infinite(N)) {
    N <- N +                        # current total number
      rpois(1, nu*timeStep) +       # new individuals (birth or migration)
      (-rpois(1, mu*X*timeStep)) +  # those susceptibles who died naturally
      (-sum(events==2)) +           # those infected who died without getting recovered
      (-sum(Zdie))                  # those recovered who died naturally
  }
  
  if(any(events>0)) {
    activeEvents <- rep(FALSE, nrow(gen))
    activeEvents[isActive] <- events>0
    gen[activeEvents, eventTime:=paste0(eventTime, ',', t)]
    gen[activeEvents, eventCode:=paste0(eventCode, ',', events[events>0])]
    
    sampledNew <- rep(FALSE, nrow(gen))
    sampledNew[isActive] <- (events==1)
    
    deadNew <- rep(FALSE, nrow(gen))
    deadNew[isActive] <- (events==2)
    
    donorsNew <- rep(FALSE, nrow(gen))
    donorsNew[isActive] <- (events==3)
    
    mutatedNew <- rep(0, nrow(gen))
    mutatedNew[isActive][events>3] <- events[events>3]-3
    
    if(sum(sampledNew) > 0) {
      if(!fadingEpidemic) {
        gen[sampledNew, sampled:=as.integer(1)]
        gen[sampledNew, tip:=(nTips+(1:sum(sampledNew)))]
        
        edgeNew <- gen[sampledNew, list(nodeP, -tip)]
        edgeLengthNew <- gen[sampledNew, tauP]
        
        nTips <- nTips+nrow(edgeNew)
        edge <- rbind(edge, cbind(as.matrix(edgeNew), rep(1, nrow(edgeNew))))
        edge.length <- c(edge.length, edgeLengthNew)
        
        nodeP <- edgeNew[, nodeP]
        while(!(length(nodeP)==1&nodeP[1]==0) & any(edge[nodeP, 3]==0)) {
          edge[nodeP, 3] <- 1
          # this works because for non-tips edge[, 2] is the same as the line-number in edge
          nodeP <- edge[nodeP, 1]
        } 
        gen[sampledNew, active:=as.integer(0)]  
      } else {
        gen[sampledNew, sampled:=as.integer(2)]
        gen[sampledNew, active:=as.integer(0)] 
      }
    }

    if(any(deadNew)) {
      gen[deadNew, active:=as.integer(0)]
      gen[deadNew, alive:=as.integer(0)]
    }
    
    # hosts in which the virus mutated at this time-step and new viral genotypes for these mutations
    if(any(mutatedNew>0)) {
      oldGene <- gen[mutatedNew>0, gene]
      newGene <- apply(cbind(oldGene, mutatedNew[mutatedNew>0]), 1, 
                       function(mut) 
                         as.integer(ifelse(mut[1]>mut[2], mut[2], mut[2]+1)))
      gen[mutatedNew>0, gene:=newGene]
    }
    
    # infect others; every infection is a node
    if(any(donorsNew)){
      gen[donorsNew, nrecips:=nrecips+1]
      
      if(!fadingEpidemic) {          
        infectNew <- infectnxNorm(gen, pe, sde, pg.init, donorsNew, eUniqForEachG) 
        
        if(!is.null(infectNew)) {
          newInfections <- TRUE
          infectNew[, id:=as.integer(nrow(gen)+(1:nrow(infectNew)))]
          infectNew[, tinf:=t]
          infectNew[, tau:=0]
          infectNew[, tauP:=0]
          infectNew[, nodeP:=as.integer(nrow(edge)+(1:nrow(infectNew)))]
          infectNew[, tip:=as.integer(0)]
          infectNew[, alive:=as.integer(1)]
          infectNew[, active:=as.integer(1)]
          infectNew[, sampled:=as.integer(0)]
          
          envD <- gen[infectNew[, idD], env]
          eD <- gen[infectNew[, idD], eSpec(e, gene)]
          gD <- gen[infectNew[, idD], gene]
          zD <- GEValues[cbind(envD, gD)] + eD
          tauD <- gen[infectNew[, idD], tau]
          
          infectNew[, envd:=envD]
          infectNew[, ed:=eD]
          infectNew[, gd:=gD]
          infectNew[, zd:=zD]
          infectNew[, taud:=tauD]
          infectNew[, nrecips:=0]
          infectNew[, eventTime:=as.character(t)]
          infectNew[, eventCode:='0']
          
          #create new edge
          edgeNew <- cbind(gen[id%in%infectNew[, idD], nodeP], nrow(edge)+(1:nrow(infectNew)), rep(0, nrow(infectNew)))
          edgeLengthNew <- gen[id%in%infectNew[, idD], tauP]
          
          gen[id%in%infectNew[, idD], nodeP:=as.integer(nrow(edge)+(1:nrow(infectNew)))]
          gen[id%in%infectNew[, idD], tauP:=0]
          
          edge <- rbind(edge, edgeNew)
          edge.length=c(edge.length, edgeLengthNew)
          
        }
      }
    }
  }
    
  res <- list(N=N, hasGen=hasGen, hasNewInfections=newInfections, edge=edge, edge.length=edge.length, nTips=nTips, zActive=zs)
  if(hasGen) {
    res$gen <- gen
  }
  if(newInfections) {
    res$genNew <- infectNew
  }
  res
}