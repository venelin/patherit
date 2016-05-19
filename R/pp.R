#' Bootstrap ANOVA 
#' @param i integer vector identifiers of the individuals
#' @param idPair integer vector showing the membership of individuals in pairs or bigger classes
#' @param z numeric vector of phenotypes
#' @return list iwth estimated rA and 95% confidence intervals
#' @details the function performs 1000 bootstraps
#' 
#' @export
rAboot <- function(i, idPair, z) {
  data <- data.table(gene=idPair, z)
  bootstrap <- boot(data=data[, unique(idPair)], statistic=function(idPP, ids) {
    rA(data=data[idPair%in%idPP[ids]])
  }, R=1000)
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
#' @param analysis logical indicating whehter the sampled mcmc should be analysed.
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
  
  N <- length(tree$tip.label)

  tipDists <- dist.nodes(tree)[1:N,1:N]
  tipDists <- sapply(1:N, function(i) sapply(1:N, function(j) if(i>j) tipDists[i,j] else NA))
  tipDists <- tipDists[!is.na(tipDists)]
  zDists <- sapply(1:N, function(i) sapply(1:N, function(j) if(i>j) abs(v[i]-v[j]) else NA))
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
