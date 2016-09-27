# bijection NxN->N
cantorPairingNumber <- function(a, b, ordered=TRUE) {
  a <- as.integer(a)
  b <- as.integer(b)
  if(!ordered) {
    c <- a
    c[a>b] <- b[a>b]
    b[a>b] <- a[a>b]
    a <- c
  }
  ((a+b)*(a+b+1)/2+b)
}

#' Extract tip pairs from a phylogeny.
#' @param tree a phylo object
#' @param threshold numeric indicating the maximum pair distance that should be allowed in the returned pairs
#' @param firstN integer indicating whether only the nearest firstN pairs should be returned
#' @return a data.table with four columns:
#' i, j : integers - the members of each tip-pair. For each entry (i,j) a symmetric entry (j,i) is present;
#'    to obtain the corresponding tip labels in the tree use tree$tip.label[i] and tree$tip.label[j] respectively.
#' tau: the phylogenetic distance between i and j
#' idPair: the unique identifier of each pair.
#' @export
extractTipPairs <- function(tree, threshold=Inf, firstN=Inf) {
  N <- length(tree$tip.label)
  tipDists <- dist.nodes(tree)[1:N, 1:N]
  diag(tipDists) <- Inf
  pairs <- as.data.table(do.call(rbind, lapply(1:N, function(i) {
    cbind(i, (1:N)[-i], tipDists[i, -i])
  })))
  setnames(pairs, c('i', 'j', 'tau'))
  pairs[, idPair:=cantorPairingNumber(i, j, ordered=FALSE)]
  if(!is.infinite(threshold)) {
    pairs[tau<=threshold][order(idPair)]
  } else if(!is.infinite(firstN)) {
    pairs[tau<=order(tau)[min(length(tau), firstN)]][order(idPair)]
  } else {
    pairs[order(idPair)]
  }
}

#' Extract phylogenetic pairs from a phylogeny or mutually closest pairs from a distance matrix.
#' @param tree a phylo object; Can be left NULL, in which case the tipDists matrix should be specified;
#' @param tipDists a N by N numeric matrix containing the tip-distances in the tree; can be NULL, in which
#' case it is going to be calculated; if specified, the tree-parameter will only be used for calculation
#' of vcvMat and no validation of tip-distances in the tree with the tipDists matrix will be done.
#' @param threshold numeric indicating the maximum pair distance that should be allowed in the returned pairs
#' @param firstN integer indicating whether only the nearest firstN pairs should be returned
#' @param vcvMat a N by N matrix where N is the number of tips in the tree. The diagonal elements of this
#' matrix represent the root-tip distances; the off-diagonal elements represent the distance from the
#' root to the most recent common ancestor of each couple of tips. This matrix is calculated by the 
#' vcv function from the ape package. If tree is specified, the matrix can be specified only for the 
#' sake of saving calculation time when a call to vcv has already been executed on the tree; 
#' If the tree is not specified (i.e. a tipDists-matrix is given), then it is assumed that the tree is 
#' ultrametric with length equal to t=max(tipDists)/2 and vcvMat[i,j] is calculated as 
#' vcvMat[i,j]=t-tipDists[i,j]/2.
#' @return a data.table with five columns:
#' i, j : integers - the members of each phylogenetic pair. For each entry (i,j) a symmetric entry (j,i) is present; 
#'    to obtain the corresponding tip labels in the tree use tree$tip.label[i] and tree$tip.label[j] respectively.
#' tau: the phylogenetic distance between i and j
#' t: the root-tip distance of the mrca of i and j
#' idPair: the unique identifier of each pair
#' 
#' @details Phylogenetic pairs represent pairs of tips each of which is the other's nearest neighbor-tip in
#' the phylogenty according to patristic (phylogenetic distance)
#' @export
extractPP <- function(tree=NULL, tipDists=NULL, threshold=Inf, firstN=Inf, vcvMat=NULL) {
  if(!is.null(tipDists)) {
    if(ncol(tipDists)!=nrow(tipDists)) {
      stop('tipDists is not a square matrix.')
    }
    N <- nrow(tipDists)
    if(is.null(tree)) {
      if(is.null(vcvMat)) {
        # assuming ultrametric tree of length max(tipDists)/2
        vcvMat <- max(tipDists)/2-tipDists/2
      }
    } else {
      if(is.null(vcvMat)) {
        vcvMat <- vcv(tree)
      }
    }
  } else {
    if(is.null(tree)) {
      stop('Parameter tree has to be specified in case of NULL tipDists.')
    }
    N <- length(tree$tip.label) 
    if(is.null(vcvMat)) {
      vcvMat=vcv(tree)
    }
  }
  
  ttips <- diag(vcvMat)
  pairs <- t(sapply(1:N, function(i) {
    vcvMati <- vcvMat[, i]
    if(is.null(tipDists)) {
      tipDistsi <- ttips[i]+ttips-2*vcvMati
    } else {
      tipDistsi <- tipDists[i, ]
    }
    
    tipDistsi[i] <- Inf  
    j <- which.min(tipDistsi)[1]; 
    c(i, j, tipDistsi[j], vcvMat[i,j])
  }))
  colnames(pairs) <- c('i', 'j', 'tau', 't')
  pp <- data.table(i=as.integer(pairs[, 'i']), j=as.integer(pairs[, 'j']), 
                   tau=pairs[, 'tau'], t=pairs[, 't'], 
                   idPair=apply(pairs, 1, function(p) {
                     if(p['i'] == pairs[p['j'], 'j']) {
                       cantorPairingNumber(p['i'], p['j'], ordered=FALSE)
                     } else {
                       NA
                     }
                   }))[!is.na(idPair)]
  
  if(!is.infinite(threshold)) {
    pp[tau<=threshold][order(idPair)]
  } else if(!is.infinite(firstN)) {
    pp[tau<=order(tau)[min(length(tau), firstN)]][order(idPair)]
  } else {
    pp[order(idPair)]
  }
}


#' Bootstrap ANOVA 
#' @param i integer vector identifiers of the individuals
#' @param idPair integer vector showing the membership of individuals in pairs or bigger classes
#' @param z numeric vector of phenotypes
#' @param bootstraps integer the number of bootstrap samples to perform (default 1000)
#' @return list iwth estimated rA and 95% confidence intervals
#' @details the function performs 1000 bootstraps
#' 
#' @export
rAboot <- function(i, idPair, z, bootstraps=1000) {
  if(length(z) < 2) {
    return(list(tips=list(i), bootstrap=NULL, 
                bCI95lower=NA,
                bCI95upper=NA,
                rA=NA, 
                CI95lower=NA, CI95upper=NA, 
                sigma2G=NA, sigma2z=NA, n=NA, 
                N=length(z), K=0))
  } else {
    data <- data.table(gene=idPair, z)
    bootstrap <- boot(data=data[, unique(idPair)], statistic=function(idPP, ids) {
      rA(data=data[idPair%in%idPP[ids]])
    }, R=bootstraps)
    aovReport=rA(data=data, report=TRUE)
    
    bCI <- try(boot.ci(bootstrap, type='basic'), silent=TRUE) 
    
    bCI95lower <- ifelse(class(bCI)!='try-error', bCI$basic[4], NA)
    bCI95upper <- ifelse(class(bCI)!='try-error', bCI$basic[5], NA)
    
    list(tips=list(i), bootstrap=list(bootstrap),
         bCI95lower=bCI95lower,
         bCI95upper=bCI95upper,
         rA=aovReport$H2aov, 
         CI95lower=aovReport$CI95lower, CI95upper=aovReport$CI95upper, 
         sigma2G=aovReport$sigma2G, sigma2z=aovReport$sigma2G+aovReport$sigma2E, n=aovReport$n0, 
         N=aovReport$N, K=aovReport$K)  
  }
  
}

#' Intra-class correlation (ICC) analysis of (closest) pairs in a phylogeny (CPP) or in a distance
#' matrix
#' @param z numeric vector of phenotypes
#' @param tree a phylo object with tip-labels corresponding to the entries in z; can be NULL, in which
#' case the dists matrix is considered; If both tree and dists are specified, dists is used.
#' @param dists a N by N numerical matrix or NULL; If both tree and dists are specified, dists is used.
#' @param CPPthr numeric indicating the maximum phylogenetic distance separating a 
#' closest phylogenetic pair
#' @param seed integer indicating if set.seed should be called with the specific seed. Specify NA to 
#' avoid calling set.seed. Defaults to 10.
#' @param ruleOutXIQR numeric a multiplier used to define outlier CPPs as defined in the referenced
#' article.
#' @param zName,treeName,distsName characters used as member-names when the parameter z is a list; default
#' values are 'z', 'tree' and 'dists'. 
#' 
#' @details The names of z are not considered. It is assumed that the elements of z are in the order 1:N
#' in tree or in the order of rows (columns) in the dists matrix
#' @export
analyseCPPs <- function(z, tree=NULL, dists=NULL, CPPthr=10^-4, seed=10, ruleOutXIQR=1.5, zName='z', treeName='tree', distsName='dists') {
  if(is.list(z)) {
    p <- z
    z <- p[[zName]]
    if(is.null(z)) {
      z <- p[['v']]
    }
    tree <- p[[treeName]]
    dists <- p[[distsName]]
    
    if(is.null(z)|(is.null(tree)&is.null(dists))) {
      stop('If a list is supplied as argument z, this list should contain a vector of trait values named "z" or zName and either a phylo-object named "tree" or treeName or a dists matrix named "dists" or distsName.')
    }
  }

  
  if(!is.null(dists)) {
    if(ncol(dists)!=nrow(dists)) {
      stop('dists is not a square matrix.')
    }
    N <- nrow(dists)
    names(z) <- rownames(dists)
  } else {
    N <- length(tree$tip.label)
    names(z) <- tree$tip.label
  }
  
  pp <- extractPP(tree=tree, tipDists=dists)
  setkey(pp, i)
  
  # an environment is needed to access z if there is a column z already in pp.
  env <- environment()
  pp[, z:=env$z[i]]
  
  pp[!is.na(idPair), deltaz:=abs(z[1]-z[2]), by=idPair]

  # tree.noCPP <- drop.tip(tree, tip=pp[tau<=CPPthr, i])
  # v.noCPP <- z[tree.noCPP$tip.label]
  # 
  # tree.noOutl <- drop.tip(tree, 
  #                                tip=pp[tau<=CPPthr & deltaz>{
  #                                    q=quantile(unique(deltaz[tau<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])
  #                                  }, i])
  # 
  # v.noOutl <- z[tree.noOutl$tip.label]
  # 
  if(!is.na(seed)) {
    set.seed(seed)
  }
  
  analysis.CPP <- c(list(pp=pp[tau<=CPPthr]), with(pp[tau<=CPPthr], rAboot(i, idPair, z)))
  analysis.PP <- c(list(pp=pp), with(pp, rAboot(i, idPair, z)))
  
  analysis.CPP.noOutl <- c(list(pp=pp[tau<=CPPthr & 
                               deltaz<={
                                 q=quantile(unique(deltaz[tau<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}]),
                  with(pp[tau<=CPPthr & 
                            deltaz<={
                              q=quantile(unique(deltaz[tau<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}], 
                       rAboot(i, idPair, z)))
  
  analysis.PP.noOutl <- c(list(pp=pp[deltaz<={
    q=quantile(unique(deltaz[tau<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}]),
    with(pp[deltaz<={q=quantile(unique(deltaz[tau<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}],
         rAboot(i, idPair, z)))
  
  list(N=N, pp=pp, 
       #tree.noCPP=tree.noCPP, z.noCPP=v.noCPP,
       #tree.noOutl=tree.noOutl, z.noOutl=v.noOutl, 
       analysis.CPP=analysis.CPP, analysis.PP=analysis.PP, 
       analysis.CPP.noOutl=analysis.CPP.noOutl, analysis.PP.noOutl=analysis.PP.noOutl)
}

#' Scatter plot of phylogenetic pa
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
#' @export
scatterPlotPPs <- function(z, tree, ppAnalysis, CPPthr=10^-4, ruleOutXIQR=1.5, zName='z', treeName='tree', xlab=expression(lg(tau)~"[lg(subst. per site)]"), ylab=expression(group("|",Delta~lg(spVL),"|")),  xlim=c(-6,0), ylim=c(0,5)) {
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
  rootTipDists <- nodeTimes(tree)[1:N]
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
