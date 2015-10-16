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
#' the function rateTransTemplate_32 for a 2-loci scenario with 3 alleles for the first locus and two alleles for 
#' the second locus. rateTransTemplate[i,j] should be the factor for transition rate from state i to state j. By
#' definition of transition rate matrices the transition from a state to itself is always the negative sum of the 
#' transition rates to the other states. 
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
#' @return A list of objects 
#' @export
simulateEpidemic <- function(Ninit=Inf, nu=0, mu=1/850, 
                             pe, sde, pg.init, GEValues, rateContact, rateInfect, rateDie, rateSample, 
                             rateMutate, rateTransTemplate,
                             eUniqForEachG=FALSE, selectWithinHost=FALSE, timeStep=.1, 
                             maxTime=1200, maxNTips=8000, expandTimeAfterMaxNTips=1.5,
                             reportInterval=12,  ...) {
  nGtps <- length(pg.init)
  nEnvs <- length(pe)
  
  largs <- list(pe, sde, pg.init, GEValues, rateContact, rateInfect, rateDie, rateSample, 
                rateMutate, rateTransTemplate, 
                eUniqForEachG, selectWithinHost, timeStep)
  names(largs) <- c('pe', 'sde', 'pg.init', 'GEValues', 'rateContact', 'rateInfect', 'rateDie', 'rateSample', 
                    'rateMutate', 'rateTransTemplate',
                    'eUniqForEachG', ' selectWithinHost', 'timeStep')
  print('Simulating epidemic on:') 
  for(i in names(largs)) {
    cat(i, ':\n')
    print(largs[[i]])
  }
  
  counts <- NULL
  countNames <- c('time', 'Total',
                  'N', 'X', 'Y', 'Z', 'nTips', 'D', 'meanActive', 'sdActive',
                  'qActive0%', 'qActive2.5%', 'qActive25%', 'qActive50%', 'qActive75%', 'qActive97.5%', 'qActive100%', 
                  'H2A', 'H2bAtInfectionA', 'betaAtInfectionA',
                  'H2S', 'H2bAtInfectionS', 'betaAtInfectionS', 
                  'H2', 'H2bAtInfection', 'betaAtInfection', 
                  'H2bA', 'betaA',
                  'H2bS', 'betaS',
                  'H2b', 'beta', 
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
      l <- stepContTimenxNorm(l$N, nu, mu, NULL, pe=pe, sde=sde, pg.init=pg.init, GEValues=GEValues, 
                              rateContact=rateContact, rateInfect=rateInfect, rateDie=rateDie, rateSample=rateSample, 
                              rateMutate=rateMutate, rateTransTemplate, 
                              eUniqForEachG=eUniqForEachG, selectWithinHost=selectWithinHost, 
                              timeStep=timeStep, t=iter, fadingEpidemic=fadingEpidemic)
      if(l$hasNewInfections) {
        gen <- rbind(l$gen, l$genNew)
      } else {
        gen <- l$gen
      }
    } else {
      l <- stepContTimenxNorm(l$N, nu, mu, gen, pe=pe, sde=sde, pg.init=pg.init, GEValues=GEValues, 
                              rateContact=rateContact, rateInfect=rateInfect, rateDie=rateDie, rateSample=rateSample, 
                              rateMutate=rateMutate, rateTransTemplate,
                              eUniqForEachG=eUniqForEachG, selectWithinHost=selectWithinHost,
                              edge=l$edge, edge.length=l$edge.length, nTips=l$nTips, timeStep=timeStep, t=iter, 
                              fadingEpidemic=fadingEpidemic)
      if(l$hasNewInfections) {
        gen <- rbind(gen, l$genNew)
      }
    }
    
    countActive <- nrow(gen[active==1])
    if(fadingEpidemic & countActive==0) {
      done=TRUE
    }
    if(!fadingEpidemic & !expandingTime & (l$nTips>=maxNTips | iter*timeStep>maxTime/expandTimeAfterMaxNTips)) {
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
               l$N, l$N-nrow(gen[alive==1]), nrow(gen[active==1]), nrow(gen[alive==1&active==0]), l$nTips, nrow(gen[alive==0]), 
               mean(l$zActive), sd(l$zActive), quantile(l$zActive, probs=c(0, .025, .25, .5, .75, .975, 1)), 
               estimH2(epid, GEValues, FALSE, TRUE), estimH2b(epid, GEValues, TRUE, FALSE, TRUE), estimBeta(epid, GEValues, TRUE, FALSE, TRUE), 
               estimH2(epid, GEValues, TRUE, FALSE), estimH2b(epid, GEValues, TRUE, TRUE, FALSE), estimBeta(epid, GEValues, TRUE, TRUE, FALSE), 
               estimH2(epid, GEValues, FALSE, FALSE), estimH2b(epid, GEValues, TRUE, FALSE, FALSE), estimBeta(epid, GEValues, TRUE, FALSE, FALSE), 
               estimH2b(epid, GEValues, FALSE, FALSE, TRUE), estimBeta(epid, GEValues, FALSE, FALSE, TRUE),
               estimH2b(epid, GEValues, FALSE, TRUE, FALSE), estimBeta(epid, GEValues, FALSE, TRUE, FALSE),
               estimH2b(epid, GEValues, FALSE, FALSE, FALSE), estimBeta(epid, GEValues, FALSE, FALSE, FALSE),
               nrow(extractDRCouples(epid, FALSE, TRUE)), nrow(extractDRCouples(epid, TRUE, FALSE)), nrow(extractDRCouples(epid, FALSE, FALSE))
               )
      
      for(k in 1:nEnvs) {
        for(j in 1:nGtps) {
          cnt <- c(cnt, nrow(gen[active==1&env==k&gene==j]))  
        }
      }
      for(k in 1:nEnvs) {
        for(j in 1:nGtps) {
          cnt <- c(cnt, nrow(gen[sampled==1&env==k&gene==j]))  
        }
      }
      
      if(is.null(counts)) {
        counts <- cnt
      } else {
        counts <- rbind(counts, cnt)
      }
    
      if((iter*timeStep)%%reportInterval==0) {
        cat('### Time: ', iter*timeStep, '\n')
        names(cnt) <- countNames
        print(c(cnt[1:7], cnt[c('H2A', 'betaA', 'nCouplesA', 'H2S', 'H2bAtInfectionA', 'betaAtInfectionA', 'H2bA')], 
                cnt[-(1:35)]/as.vector(matrix(cnt[c('Y', 'Z')], nrow=(length(cnt)-35)/2, ncol=2, byrow=TRUE))))
      }
    }
  }
  
  colnames(counts) <- countNames
  
  c(l[c('N', 'edge', 'edge.length', 'nTips')], 
    list(gen=gen, time=timeStep*iterActive, timeFading=timeStep*iter, timeStep=timeStep, rateContact=rateContact, 
         rateInfect=rateInfect, rateDie=rateDie, rateSample=rateSample, rateMutate=rateMutate, eUniqForEachG=eUniqForEachG, 
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

#' Transmission chain between sampled individuals in an epidemic
#' @param epidemic list returned by simulateEpidemic
#' @param tMin,tMax numeric, time interval for which to measure the heritability. if sampledOnly is TRUE this would be 
#'    the times of sampling for the individuals; Currently these two paraemters are not implemented.
#' @export
makeTree <- function(epidemic, tMin=0, tMax=Inf) {
  gen <- copy(epidemic$gen)
  gen[, timeObs:=(tinf+tau)*epidemic$timeStep]
  
  edge <- epidemic$edge
  edge.length <- epidemic$edge.length
  nTips <- epidemic$nTips
  timeStep <- epidemic$timeStep
  if(nTips <= 1)
    return(NULL)
  cat('Generating tree: nTips=',nTips, ', number of uncollapsed edges=', nrow(edge), '\n')
  if(ncol(edge)==3) {
    edge.length <- edge.length[edge[,3]==1]
    edge <- edge[edge[,3]==1, 1:2]
  }
  edge <- apply(edge, 1:2, function(x) if(x<0) -x else x+nTips+1)
  nodesUnique <- sort(unique(as.vector(edge)))
  
  Nnode <- length(nodesUnique)-nTips
  
  edge[, 1] <- match(edge[, 1], nodesUnique)
  edge[, 2] <- match(edge[, 2], nodesUnique)
  
  tipl <- gen[tip!=0, list(id, tip, timeObs)]
  setkey(tipl, tip)
  colnames(edge) <- NULL
  obj <- list(edge=edge, edge.length=edge.length, tip.label=tipl[,id], node.label=nodesUnique[-(1:nTips)]-nTips-1, Nnode=Nnode)
  class(obj) <- 'phylo'
  if(tMin!=0|!is.infinite(tMax)) {
    cat('Dropping tips outside of time interval :', tMin, ',', tMax, '\n')
    
    #obj <- drop.tip(obj, tip=c(tipl[timeObs<tMin, tip], tipl[timeObs>tMax, tip]), trim.internal=TRUE)
  }
  obj <- collapse.singles(obj)
  obj$edge.length <- obj$edge.length*timeStep
  obj
}

#' @export
extractDRCouples <- function(epidemic, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf) {
  gen <- copy(epidemic$gen[sampled >= ifelse(sampledOnly, 1, 0) & sampled<=1 & active>=ifelse(activeOnly, 1, 0), ])
  
  if(sampledOnly) {
    gen[, timeObs:=(tinf+tau)*epidemic$timeStep]
  } else {
    gen[, timeObs:=tinf*epidemic$timeStep]
  }
  gen <- gen[timeObs>=tMin&timeObs<=tMax]
  gen[idD%in%id, list(idD, id, envd, gd, ed, env, gene, e)]
}

#' Regression slope of recipient on donor values in an epidemic
#' @param sampledOnly logical, indicating if only sampled individuals should be included in the heritability calculation
#' @param tMin,tMax numeric, time interval for which to measure the heritability. if sampledOnly is TRUE this would be 
#'    the times of sampling for the individuals; otherwise, this would be the times of infection. 
#' @export
estimBeta <- function(epidemic, GEValues, atInfection=FALSE, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, 
                      lastN=Inf, report=FALSE) {
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
      setkey(gen, id)
      
      couples <- gen[idD%in%id, list(idD, id, envd, gd, ed, env, gene, e, tauR=tau, tauDAtInf=taud)]
      couples[, tauD:=gen[J(couples[, idD]), tau]-tauDAtInf]
      couples[, taum:=tauD+tauR]
      
      if(nrow(couples) > 0) {
        if(atInfection) {
          couples[, zD0:=calcValue(envd, gd, ed, GEValues)] 
          couples[, zR0:=calcValue(env, gd, e, GEValues)]
          beta=couples[, cov(zD0, zR0)/var(zD0)]
          if(report) {
            couples[, GED0:=GEValues[cbind(envd, gd)]]
            couples[, GER0:=GEValues[cbind(env, gd)]]
            couples[, eD0:=ed]
            couples[, eR0:=eSpec(e, gd)]
            list(beta=beta, couples=couples)
          } else {
            beta
          }
        } else {
          couples[, zD:=gen[J(couples[, idD]), calcValue(env, gene, e, GEValues)]]
          couples[, zR:=calcValue(env, gene, e, GEValues)]
          beta=couples[, cov(zD, zR)/var(zD)]
          if(report) {
            couples[, GED:=gen[J(couples[, idD]), GEValues[cbind(env, gene)]]]
            couples[, GER:=GEValues[cbind(env, gene)]]
            couples[, eD:=gen[J(couples[, idD]), eSpec(e, gene)]]
            couples[, eR:=eSpec(e, gene)]
            list(beta=beta, couples=couples)
          } else {
            beta
          }
        }  
      } else {
        NA
      }           
    } else {
      NA
    }
  } 
}

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

#' Breeding heritability of a pathogen trait
#' @param epidemic a list of objects returned from simulateEpidemic
#' @param GEValues genotype-environment trait values
#' @param sampledOnly logical, indicating if only sampled individuals should be included in the heritability calculation
#' @param tMin,tMax numeric, time interval for which to measure the heritability. if sampledOnly is TRUE this would be 
#' the times of sampling for the individuals; otherwise, this would be the times of infection.
#' @export
estimH2b <- function(epidemic, GEValues, atInfection=FALSE, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, lastN=Inf, report=FALSE) {
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
      setkey(gen, id)
      
      couples <- gen[idD%in%id, list(idD, id, envd, gd, ed, env, gene, e, tauR=tau, tauDAtInf=taud)]
      couples[, tauD:=gen[J(couples[, idD]), tau]-tauDAtInf]
      couples[, taum:=tauD+tauR]
      
      if(nrow(couples) > 0) {
        if(atInfection) {
          couples[, zD0:=calcValue(envd, gd, ed, GEValues)]
          couples[, zR0:=calcValue(env, gd, e, GEValues)]
          couples[, bD0:=mean(zR0), by=gd]
          H2b0 <- couples[, var(bD0)/var(zD0)]
          if(report) {
            list(H2b0=H2b0, couples=couples)
          } else {
            H2b0
          }
        } else {
          couples[, zD:=gen[J(couples[, idD]), calcValue(env, gene, e, GEValues)]]
          couples[, zR:=calcValue(env, gene, e, GEValues)]
          couples[, bD:=mean(zR), by=gd]
          H2b <- couples[, var(bD)/var(zD)]
          if(report) {
            list(H2b=H2b, couples=couples)
          } else {
            H2b
          }
        }  
      } else {
        NA
      }           
    } else {
      NA
    }
  }
}

#' Broad-sense heritability of a pathogen trait
#' @param epidemic a list of objects returned from simulateEpidemic
#' @param GEValues genotype-environment trait values
#' @param sampledOnly logical, indicating if only sampled individuals should be included in the heritability calculation
#' @param tMin,tMax numeric, time interval for which to measure the heritability. if sampledOnly is TRUE this would be 
#' the times of sampling for the individuals; otherwise, this would be the times of infection.
#' @return numeric indicating the estimated heritability for infected individuals in the population. 
#'  
#' @export
estimH2 <- function(epidemic, GEValues, sampledOnly=TRUE, activeOnly=FALSE, tMin=0, tMax=Inf, lastN=Inf) {
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
    if(nrow(gen) >0) {
      gen[, z:=calcValue(env, gene, e, GEValues)]
      gen[, G:=mean(z), by=gene]
      gen[, var(G)/var(z)]  
    } else {
      NA
    }
  }
}

# simulating evolution on continuous time scale 

# sde : numeric or vector of standard deviations of the environmental effects for each allele
# N : number of individuals
# nGtps : number of alleles 
# sde : a numeric or a numeric matrix, sde[env, g]: standard dev of the special env. effect for <env, g>
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

sampleEvent <- function(rates, timeStep) {
  lambdas <- rowSums(rates)
  probs0 <- exp(-lambdas*timeStep)
 
  probs <- cbind(probs0, rates/lambdas*(1-probs0))
  for(i in 2:ncol(probs)) {
    probs[, i] <- probs[, i-1]+probs[, i]
  }
  
  random <- runif(nrow(rates))
  probs_random <- (probs-random)>=0
  
  events <- rep(0, nrow(rates))
  for(e in ncol(probs_random):1) {
    events[probs_random[, e]] <- e-1
  }
  events
}

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
#' @export
calcValue <- function(env, gene, e, GEValues) {
  GEValues[cbind(env, gene)] + eSpec(e, gene)
}

#'
#' state transition template matrix (main diagonal omitted)
#' element [i, j] is a multiplier for the rate of transition 
#' from state i to state j+k, where k=1 if j>=i and k=0 otherwise
#' @export
rateTransTemplate_32 <- function() {
  rbind(
    c(1, 1/2, 0, 1/2, 0),
    c(1, 0, 1/2, 0, 1/2),
    c(1/2, 0, 1, 1/2, 0),
    c(0, 1/2, 1, 0, 1/2),
    c(1/2, 0, 1/2, 0, 1),
    c(0, 1/2, 0, 1/2, 1))
}

# state-transition rates
rateStateTransition <- function(z, rates, eMatrix, es, envs, genes, GEValues, rateTransTemplate, selectWithinHost) {
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
    
    matrix(sgnother*rateTransTemplate[genes, ]*rates, nrow=length(z))
  } else {
    matrix(rateTransTemplate[genes, ]*rates, nrow=length(z))
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
  X <- N-nrow(gen[alive==1])
  
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
  envs <- gen[isActive, env]
  genes <- gen[isActive, gene]
  
  eMatrix <- gen[isActive, do.call(rbind, e)]
  es <- gen[isActive, eSpec(eMatrix, gene)]
  
  zs <- GEValues[cbind(envs, genes)]+es
  
  ratesSample <- rep(rateSample, length(zs))
  ratesDie <-  rateDie(zs, mu)        
  ratesInfect <- rateInfect(zs, rateRiskyContact)   
  
  ratesMutate <- rateMutate(GEValues, es, envs, genes)
  
  ratesTrans <- rateStateTransition(zs, ratesMutate, eMatrix, es, envs, genes, 
                                    GEValues, rateTransTemplate, selectWithinHost)
  
  rates <- cbind(ratesSample, ratesDie, ratesInfect, ratesTrans)
  
  events <- sampleEvent(rates, timeStep)
  
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
    
    if(any(sampledNew)) {
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
            nodeP <- edge[nodeP, 1]
          } 
          gen[sampledNew, active:=as.integer(0)]  
        } else {
          gen[sampledNew, sampled:=as.integer(2)]
          gen[sampledNew, active:=as.integer(0)] 
        }
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