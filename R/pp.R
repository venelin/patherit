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
#' d: the phylogenetic distance between i and j
#' idPair: the unique identifier of each pair.
#' @export
extractTipPairs <- function(tree, threshold=Inf, firstN=Inf) {
  N <- length(tree$tip.label)
  tipDists <- dist.nodes(tree)[1:N, 1:N]
  diag(tipDists) <- Inf
  pairs <- as.data.table(do.call(rbind, lapply(1:N, function(i) {
    cbind(i, (1:N)[-i], tipDists[i, -i])
  })))
  setnames(pairs, c('i', 'j', 'd'))
  pairs[, idPair:=cantorPairingNumber(i, j, ordered=FALSE)]
  if(!is.infinite(threshold)) {
    pairs[d<=threshold][order(idPair)]
  } else if(!is.infinite(firstN)) {
    pairs[d<=order(d)[min(length(d), firstN)]][order(idPair)]
  } else {
    pairs[order(idPair)]
  }
}

#' Extract phylogenetic pairs from a phylogeny.
#' @param tree a phylo object
#' @param threshold numeric indicating the maximum pair distance that should be allowed in the returned pairs
#' @param firstN integer indicating whether only the nearest firstN pairs should be returned
#' @param vcvMat a N by N matrix where N is the number of tips in the tree. The diagonal elements of this
#' matrix represent the root-tip distances; the off-diagonal elements represent the distance from the
#' root to the most recent common ancestor of each couple of tips. This matrix is calculated by the 
#' vcv function from the ape package and can be passed as parameter only for the sake of saving 
#' calculation time if a call to vcv has already been executed on the tree.
#' @return a data.table with four columns:
#' i, j : integers - the members of each phylogenetic pair. For each entry (i,j) a symmetric entry (j,i) is present; 
#'    to obtain the corresponding tip labels in the tree use tree$tip.label[i] and tree$tip.label[j] respectively.
#' tau: the phylogenetic distance between i and j
#' t: the root-tip distance of the mrca of i and j
#' idPair: the unique identifier of each pair
#' 
#' @details Phylogenetic pairs represent pairs of tips each of which is the other's nearest neighbor-tip in
#' the phylogenty according to patristic (phylogenetic distance)
#' @export
extractPP <- function(tree, threshold=Inf, firstN=Inf, vcvMat=vcv(tree)) {
  N <- length(tree$tip.label)
  
  ttips <- diag(vcvMat)
  pairs <- t(sapply(1:N, function(i) {
    vcvMati <- vcvMat[, i]
    tipDistsi <- ttips[i]+ttips-2*vcvMati
    tipDistsi[i] <- Inf
    j <- which.min(tipDistsi)[1]; 
    c(i, j, tipDistsi[j], vcvMat[i,j])
  }))
  colnames(pairs) <- c('i', 'j', 'tau', 't')
  pp <- data.table(i=as.integer(pairs[, 'i']), j=as.integer(pairs[, 'j']), 
                   tau=pairs[, 'tau'], t=pairs[, 't'], 
                   idPair=apply(pairs, 1, function(p) if(p['i'] == pairs[p['j'], 'j']) cantorPairingNumber(p['i'], p['j'], ordered=FALSE) else NA))[!is.na(idPair)]
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
  data <- data.table(gene=idPair, z)
  bootstrap <- boot(data=data[, unique(idPair)], statistic=function(idPP, ids) {
    rA(data=data[idPair%in%idPP[ids]])
  }, R=bootstraps)
  aovReport=rA(data=data, report=TRUE)
  list(tips=list(i), bootstrap=list(bootstrap),
       bCI95lower=boot.ci(bootstrap, type='basic')$basic[4], 
       bCI95upper=boot.ci(bootstrap, type='basic')$basic[5],
       rA=aovReport$H2aov, 
       CI95lower=aovReport$CI95lower, CI95upper=aovReport$CI95upper, 
       sigma2G=aovReport$sigma2G, sigma2z=aovReport$sigma2G+aovReport$sigma2E, n=aovReport$n0, 
       N=aovReport$N, K=aovReport$K)
}

#' (Closest) phylogenetic pair analysis
#' @param z numeric vector of phenotypes
#' @param tree a phylo object with tip-labels corresponding to the entries in z
#' @param CPPthr numeric indicating the maximum phylogenetic distance separating a 
#' closest phylogenetic pair
#' @param seed integer indicating if set.seed should be called with the specific seed. Specify NA to 
#' avoid calling set.seed. Defaults to 10.
#' @param ruleOutXIQR numeric a multiplier used to define outlier CPPs as defined in the referenced
#' article.
#' @param zName,treeName characters used when the parameter z is a list; indicate the 
#' @export
analyseCPPs <- function(z, tree, CPPthr=10^-4, seed=10, ruleOutXIQR=1.5, zName='z', treeName='tree') {
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
  names(z) <- tree$tip.label
  
  N <- length(tree$tip.label)

  tipDists <- dist.nodes(tree)[1:N,1:N]
  tipDists <- sapply(1:N, function(i) sapply(1:N, function(j) if(i>j) tipDists[i,j] else NA))
  tipDists <- tipDists[!is.na(tipDists)]
  zDists <- sapply(1:N, function(i) sapply(1:N, function(j) if(i>j) abs(z[i]-z[j]) else NA))
  zDists <- zDists[!is.na(zDists)]
  pp <- extractPP(tree)
  setkey(pp, i)
  
  # a copy is needed only for the next operation
  v <- z
  pp[, z:=v[i]]
  
  pp[!is.na(idPair), deltaz:=abs(z[1]-z[2]), by=idPair]

  tree.noCPP <- drop.tip(tree, tip=pp[d<=CPPthr, i])
  v.noCPP <- v[tree.noCPP$tip.label]
  
  tree.noOutl <- drop.tip(tree, 
                                 tip=pp[d<=CPPthr & deltaz>{
                                     q=quantile(unique(deltaz[d<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])
                                   }, i])
  
  v.noOutl <- v[tree.noOutl$tip.label]
  
  if(!is.na(seed)) {
    set.seed(seed)
  }
  
  analysis.CPP <- c(list(pp=pp[d<=CPPthr]), with(pp[d<=CPPthr], rAboot(i, idPair, z)))
  analysis.PP <- c(list(pp=pp), with(pp, rAboot(i, idPair, z)))
  
  analysis.CPP.noOutl <- c(list(pp=pp[d<=CPPthr & 
                               deltaz<={
                                 q=quantile(unique(deltaz[d<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}]),
                  with(pp[d<=CPPthr & 
                            deltaz<={
                              q=quantile(unique(deltaz[d<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}], 
                       rAboot(i, idPair, z)))
  
  analysis.PP.noOutl <- c(list(pp=pp[deltaz<={
    q=quantile(unique(deltaz[d<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}]),
    with(pp[deltaz<={q=quantile(unique(deltaz[d<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}],
         rAboot(i, idPair, z)))
  
  list(tipDists=tipDists, zDists=zDists, pp=pp, tree.noCPP=tree.noCPP, z.noCPP=v.noCPP,
       tree.noOutl=tree.noOutl, z.noOutl=v.noOutl, 
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
  pp[d<=CPPthr, points(x=log10(d), y=deltaz, pch=20, cex=0.25, col=adjustcolor('magenta', alpha.f=0.3))]
  pp[d<=CPPthr & 
       deltaz>{
         q=quantile(unique(deltaz[d<=CPPthr])); q[4]+ruleOutXIQR*(q[4]-q[2])}, 
     points(x=log10(d), y=deltaz, col='blue', pch=20, cex=0.4)]
}

#' Box plots of trait values along a tree
#' @param nGroups integer the number of groups
#' @param ... additional parameters passed to boxplot
#' @export
boxplotTraitAlongTree <- function(z, tree, nGroups=15, ...) {
  groups <- groupByRootDist(tree, nGroups=nGroups)
  rootTipDistGroups <- groups$rootTipDistGroups
  boxes <- list()

  for(i in 1:length(rootTipDistGroups)) {
    boxes[[as.character(rootTipDistGroups[[i]])]] <-
      c(boxes[[as.character(rootTipDistGroups[[i]])]], z[i])
  }

  boxes <- boxes[sort(names(boxes))]
  midDistPoints <- groups$groupMeans
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
