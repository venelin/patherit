# Implementation of the POUMM likelihood and heritability estimators

NULL

#' Expected value of the Ornstein-Uhlenbeck process conditioned on initial value
#' @param g0 numeric, initial value
#' @param alpha the strength of the selection
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param t time at which the expected value should be estimated
#' @export
mu.poumm <- function(g0, alpha, theta, sigma, t=Inf) {
  ifelse(t==Inf, 
         theta,
         theta + (g0 - theta)*exp(-alpha*t))
}

#' POUMM variance covariance matrix of tips
#' 
#' @param tree a phylo object
#' @param alpha a non-negative numeric, default is 0
#' @param sigma a non-negative numeric, default is 1
#' @param sigmae a non-negative numeric, default is 0
#' @param corr logical indicating if a correlation matrix shall be returned
#' @param tanc numerical matrix with the time-distance from the root of the 
#' tree to the mrca of each tip-couple. If NULL it will be calculated.
#' 
#' @return a variance covariance or a correlation matrix of the tips in tree.
#' @references (Hansen 1997) Stabilizing selection and the comparative analysis
#'   of adaptation.
#' @export
cov.poumm <- function(tree, alpha=0, sigma=1, sigmae=0, corr=F, tanc=NULL) {
  N <- length(tree$tip.label)
  
  # distances from the root to the mrca's of each tip-couple.
  if(is.null(tanc)) {
    tanc <- vcv(tree)
  }
  if(alpha > 0) {
    varanc <- 0.5*(1-exp(-2*alpha*tanc))/alpha
  } else {
    # limit case alpha==0
    varanc <- tanc
  }
  
  # time from the root to each tip
  ttips <- diag(tanc) # nodeTimes(tree)[1:N]
  # distance-matrix between tips
  tij <- sapply(ttips, function(.) .+ttips)-2*tanc
  
  covij <- exp(-alpha*tij)*varanc*sigma^2 + diag(rep(sigmae^2, N))
  if(corr) {
    sqvar <- sqrt(diag(covij))
    sccor <- sapply(sqvar, function(.) .*sqvar)
    covij/sccor
  } else {
    covij
  }  
}

#' Distribution of the genotypic values under a POUMM fit
#' 
#' @param tree an object of class phylo
#' @param z A numeric vector of size length(tree$tip.label) representing the trait
#'     values at the tip-nodes.
#' @export
g.poumm <- function(z, tree, fit.poumm, divideEdgesBy=fit.poumm$divideEdgesBy) {
  N <- length(tree$tip.label)

  alpha <- fit.poumm$par[1]
  theta <- fit.poumm$par[2]
  sigma <- fit.poumm$par[3]
  sigmae <- fit.poumm$par[4]
  g0 <- attr(fit.poumm$value, 'grmax')

  tree$edge.length <- tree$edge.length/divideEdgesBy

  V.g <- cov.poumm(tree, alpha, sigma, sigmae)
  V.g_1 <- chol2inv(chol(V.g))
  mu.g <- mu.poumm(g0, alpha=alpha, theta=theta, sigma=sigma, t=nodeTimes(tree)[1:N])

  V.e <- diag(sigmae^2, nrow=N, ncol=N)
  V.e_1 <- V.e
  diag(V.e_1) <- 1/diag(V.e)
  mu.e <- z

  V.g.poumm <- try(chol2inv(chol(V.g_1+V.e_1)), silent = TRUE)
  if(class(V.g.poumm)=='try-error') {
    warning(V.g.poumm)
    list(V.g=V.g, V.g_1=V.g_1, mu.g=mu.g, V.e=V.e, V.e_1=V.e_1, mu.e=mu.e)
  } else {
    mu.g.poumm <- V.g.poumm%*%(V.g_1%*%mu.g+V.e_1%*%mu.e)

    list(V.g=V.g, V.g_1=V.g_1, mu.g=mu.g, V.e=V.e, V.e_1=V.e_1, mu.e=mu.e, V.g.poumm=V.g.poumm, mu.g.poumm=mu.g.poumm)
  }
}

#' Likelihood of observed node values given a tree
#' 
#' @description
#' Calculates the conditional probability density of the tree 
#' trait values at the tips and internal nodes of the tree 
#' given the tree and the trait-value at the root according to a model, such 
#' that the probability density function of a value 
#' at a descendent node gdesc given its ancestor-value ganc and edge-length l 
#' is given by pdf(gdesc, ganc, l, ...). 
#
#' @param g A numeric vector of size length(tree$tip.label)+tree$Nnode 
#'    representing the genetic contibution to the trait data at tips, root and 
#'    internal nodes. 
#' @param tree an object of class phylo
#' @param pdf A probability density function
#' @param log Logical indicating whether a log-likelihood should be returned 
#'    instead of a likelihood. Default is TRUE.
#' @param ... additional arguments passed to the pdf function
#' @return
#'    If !nodeprobs (default) log P(tree, g | pdf), otherwise the vector of the
#'    node-value's (log-)probabilities conditioned on their ancestor's values 
#'    otherwise (root not included)
#' @export
likGTree <- function(g, tree, pdf, log=T, ...) {
  
  if (is.null(tree$edge.length)) 
    stop("tree has no branch length")
  if (any(tree$edge.length < 0)) 
    stop("at least one branch length negative")
  
  N <- dim(tree$edge)[1]            # number of edges in the tree
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  el <- tree$edge.length
  probs <- sapply(1:N, function(i) pdf(g[des[i]], g[anc[i]], el[i], ..., log=T))
  
  ifelse(log, sum(probs), exp(sum(probs)))
}

#' Likelihood of observed genetic contributions at all nodes and environmental 
#' deviations at the tips given a tree.
#' 
#' @description
#' Calculates the probability density function of the genetic trait 
#' contributions at all nodes, g, and environmental 
#' trait deviations at the tip-nodes given the tree topology and assuming that 
#' the genetic contributions evolved on the tree according to an OU process 
#' with parameters alpha, theta and sigma and the environmental deviations 
#' represent iid samples from a normal distribution with mean mue and standard 
#' deviations sigmae.
#'
#' @param g A numeric vector of size length(tree$tip.label)+tree$Nnode 
#'    representing the genetic contibution to the trait data at tips, root and 
#'    internal nodes. Note that g[length(tree$tip.label)+1] is assumed to be 
#'    the root.
#' @param tree an object of class phylo
#' @param alpha the strength of the selective constraint
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param e NULL or a numeric vector of size length(tree$tip.label) representing the 
#'    environmental deviations at the tip-nodes. If NULL the arguments mue and
#'    sigmae are ignored.
#' @param mue Mean of the environmental deviation distribution.
#' @param sigmae standard deviation of the environmental deviation distribution.
#' @param log Logical indicating whether a log-likelihood should be returned 
#'    instead of a likelihood. Default is TRUE.
#' @param ...  additional arguments passed to the pdf function
#' @return
#'    (log-)probability of g and e given tree, alpha, theta, sigma, mue, sigmae.
#' @export
likGETreeOU <- function(g, tree, alpha, theta, sigma, 
                        e=NULL, mue=0, sigmae=1, log=T) {
  lTG <- likGTree(g, tree, dCondOU, alpha=alpha, theta=theta, sigma=sigma, log=T)
  lRoot <- dStOU(g[length(tree$tip.label)+1], alpha, theta, sigma, log=T)
  lE <- 0
  if(!is.null(e)) 
    lE <- dnorm(e, mue, sigmae, log=T)
  
  ifelse(log, sum(lTG, lRoot, lE), exp(sum(lTG, lRoot, lE)))
}

#' Likelihood of observed tip-values given a tree, assuming Ornstein-Uhlenbeck
#' process for the genetic contributions along the tree and normally distributed
#' environmental deviations.
#' 
#' @description
#' Calculates the (log-)probability density of trait values at the tip-nodes 
#' given the tree and assuming that the trait value at the tips is the sum 
#' of a genetic contribution, g, that evolved on the tree according to an OU 
#' process with parameters alpha, theta, sigma and an environmental deviation,
#' e, that is distributed normally and independently between the tips of the 
#' tree.
#' Note: Without additional assumptions for the distribution of the value at 
#' the root of the tree, the likelihood is not defined at alpha=0, although 
#' this corresponds to the limiting Brownian motion process with mean value
#' theta and unit time variance sigma^2.
#'
#' @param tree an object of class phylo
#' @param z A numeric vector of size length(tree$tip.label) representing the trait
#'     values at the tip-nodes.
#' @param alpha the strength of the selection
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param sigmae the standard deviation of the environmental deviation added to the
#'     genetic contribution at each tip, by default 0, meaning no environmental
#'     deviation.
#' @param distgr a character or a numeric describing the distribution or fixed value 
#'     to be assumed for the genetic contribution at the root. If a character the 
#'     acceptable values are : 
#'       'normal'- normal distribution with mean mugr and standard deviation sigmagr; 
#'       'maxlik' - fixed value that would maximise the conditional likelihood on gr.
#' @param mugr see distgr
#' @param sigmagr see distgr
#' @param log Logical indicating whether log-likelihood should be returned instead
#'    of likelihood, default is TRUE.
#' @param pruneInfo list returned by pruneTree(tree) to be passed in explicit calls to lik.poumm. 
#' @param usempfr,maxmpfr integer indicating if and how mpfr should be used for small parameter values 
#'    (any(c(alpha, sigma, sigmae) < 0.01)). Using the mpfr package can be forced by specifying an 
#'    integer greater or equal to 2. Setting usempfr=0 disables high precision likelihood calculation.
#'    Requires the Rmpfr package. Note that when using mpfr, the time for one likelihood calculation can
#'    increase more than 100-fold. 
#' @param precbits integer specifying precision bits for mpfr. Default is 512.
#' @param debug logical, if set to TRUE some debugging information is stored in a global 
#'    list called .likVTreeOUDebug
#' @export
lik.poumm <- function(z, tree, alpha, theta, sigma, sigmae=0, 
                       distgr=c('normal', 'maxlik'), 
                       mugr=theta, 
                       sigmagr=ifelse(alpha==0&sigma==0, 0, sigma/sqrt(2*alpha)),
                       log=T, pruneInfo=pruneTree(tree),
                       usempfr=1, maxmpfr=2, precbits=512, debug=F) {
  alphaorig <- alpha
  thetaorig <- theta
  sigmaorig <- sigma
  sigmaeorig <- sigmae
  mugrorig <- mugr
  sigmagrorig <- sigmagr
  vorig <- z
  
  if(any(tree$edge.length<=0)) {
    stop('All branch lengths in tree should be positive!')
  }
  if(any(is.na(c(alpha, theta, sigma, sigmae))) | 
    any(is.nan(c(alpha, theta, sigma, sigmae)))) {
    warning('Some parameters are NA or NaN')
    ifelse(log, -Inf, 0)
  } 
  if(any(is.infinite(c(alpha, sigma, sigmae))) | # case 9
       any(c(alpha, sigma, sigmae) < 0) |        # case 10
       (alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
       all(c(alpha, sigma, sigmae) == 0)        # case 8
       ) {
    ifelse(log, -Inf, 0)
  } else {
    logpi <- log(pi)
    loge2 <- log(2)
    
    availRmpfr <- require(Rmpfr) 
    
    if(as.integer(usempfr)>=2 |
         (usempfr & ((as.double(alpha) != 0 & as.double(alpha) < 0.01) |
                       (as.double(sigma) != 0 & as.double(sigma) < 0.01) |
                       (as.double(sigmae) != 0 & as.double(sigmae) < 0.01)))) {
      usempfr <- 2 # indicate for later code that we shall convert to mpfr!
    }
     
    done <- F
    # this loop governs the use of mpfr
    while(!done & usempfr <= maxmpfr) { 
      if(availRmpfr & usempfr >= 2) {
        if(is.double(alphaorig)) alpha <- mpfr(alphaorig, precbits*2^(usempfr-2))
        if(is.double(thetaorig)) theta <- mpfr(thetaorig, precbits*2^(usempfr-2))
        if(is.double(sigmaorig)) sigma <- mpfr(sigmaorig, precbits*2^(usempfr-2))
        if(is.double(sigmaeorig)) sigmae <- mpfr(sigmaeorig, precbits*2^(usempfr-2))
        if(is.double(mugr)) mugr <- mpfr(mugr, precbits*2^(usempfr-2))
        if(is.double(sigmagr)) sigmagr <- mpfr(sigmagr, precbits*2^(usempfr-2))
        if(is.double(z)) z <- mpfr(z, precbits*2^(usempfr-2))
      }
      

        N <- length(tree$tip.label)                # number of tips
        
        if(is.null(pruneInfo)) {
          pruneInfo <- pruneTree(tree)
        }
        
        edge <- tree$edge
        
        M <- pruneInfo$M
        endingAt <- pruneInfo$endingAt
        nodesVector <- pruneInfo$nodesVector
        nodesIndex <- pruneInfo$nodesIndex
        nLevels <- pruneInfo$nLevels
        unVector <- pruneInfo$unVector
        unIndex <- pruneInfo$unIndex
        unJ <- 1
        
        t <- tree$edge.length
        
        alphasigma2 <- alpha/sigma/sigma
        theta2 <- theta^2
        sigma2 <- sigma^2
        sigmae2 <- sigmae^2
        logsigma <- log(sigma)
        logsigmae <- log(sigmae)
        
        talpha <- t*alpha
        etalpha <- exp(talpha)
        if(alpha != 0) {
          fetalpha <- alpha/(1 - etalpha)
        } else {
          fetalpha <- -1/t
        }
        
        e2talpha <- etalpha*etalpha
        limitfe2talpha <- -0.5/t
        if(alpha != 0) {
          fe2talpha <- alpha/(1 - e2talpha)
        } else {
          fe2talpha <- limitfe2talpha
        }
        
        if(availRmpfr & usempfr & usempfr < maxmpfr) {
          if(any(fe2talpha >= 0) | any(fe2talpha < limitfe2talpha)) {
            usempfr <- usempfr + 1
            next 
          }
        }
        
        log_fe2talpha <- log(-fe2talpha)
        
        fe2talphasigma2 <- fe2talpha/sigma2
        
        # for sigmae=0
        r0 <- talpha + 0.5*log_fe2talpha - 0.5*logpi - logsigma
        
        # matrix for paramsIntegForks one pif for every node
        # the following code creates a matrix for the class of alpha,
        # i.e. could be mpfr as well as double
        pif <- rep(alpha*0, M*3)
        dim(pif) <- c(M, 3)
        
        for(i in 1:nLevels) {
          nodes <- nodesVector[(nodesIndex[i]+1):nodesIndex[i+1]]
          es <- endingAt[nodes]
          nodes <- c()
          edgeEnds <- edge[es, 2]
          
          if(edge[es[1], 2] <= N) {
            # all es pointing to tips
            if(sigmae==0) {
              # no environmental deviation
              g1 <- z[edgeEnds]  # tip values
              etalphag1thetatheta <- etalpha[es]*(g1-theta)+theta
              
              pif[edgeEnds, 1] <- fe2talphasigma2[es]
              pif[edgeEnds, 2] <- -2 * etalphag1thetatheta * fe2talphasigma2[es]
              pif[edgeEnds, 3] <- r0[es] + etalphag1thetatheta^2 * fe2talphasigma2[es]
            } else {
              z1 <- z[edgeEnds]  # tip values
              
              u <- -0.5/sigmae2
              v <- z1/sigmae2
              w <- -0.5*loge2 - 0.5*logpi  - z1^2/(2*sigmae2) - logsigmae 
              
              gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
              
              pif[edgeEnds, 1] <- u / gutalphasigma2
              pif[edgeEnds, 2] <- (etalpha[es]*(2*theta*u+v)-2*theta*u) / gutalphasigma2
              pif[edgeEnds, 3] <- -0.5*log(gutalphasigma2) -
                0.25*sigma2*v^2/(fe2talpha[es]-alpha + sigma2*u) +
                alpha*theta*(theta*u-etalpha[es]*(theta*u+v)) /
                ((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + talpha[es]+w
            }
          } else {
            # edges pointing to internal nodes, for which all children nodes have been visited
            u <- pif[edgeEnds, 1]
            v <- pif[edgeEnds, 2]
            w <- pif[edgeEnds, 3]
            
            gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]

            pif[edgeEnds, 1] <- u / gutalphasigma2
            pif[edgeEnds, 2] <- (etalpha[es]*(2*theta*u+v)-2*theta*u) / gutalphasigma2
            pif[edgeEnds, 3] <- -0.5*log(gutalphasigma2) -
              0.25*sigma2*v^2/(fe2talpha[es]-alpha + sigma2*u) +
              alpha*theta*(theta*u-etalpha[es]*(theta*u+v)) /
              ((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + talpha[es]+w
          } 
          
          #update parent pifs
          while(length(es)>0) {
            un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
            unJ <- unJ+1
            pif[edge[es[un], 1], ] <- pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
            es <- es[-un]
          }
        }
        
        # at the root: P[Vtips|tree, g_root, alpha, theta, sigma, sigmae] =
        #                  abc[1]*g_root^2 + abc[2]*g_root + abc[3]
        abc <- pif[N+1, ]
        
        if(availRmpfr & usempfr & (any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3]))) {
          usempfr <- usempfr + 1
          next 
        }
        done <- T
    }
    
    if(!is.character(distgr <- distgr[1])) {
      # fixed value
      loglik <- abc[1]*distgr[1]^2 + abc[2]*distgr[1] + abc[3]
    } else {        
      if(distgr[1] == 'normal') {
        mugr2 <- mugr^2
        if(sigmagr > 0) {
          sigmagr2 <- sigmagr^2
          def <- c(-0.5/sigmagr2, 
                   mugr/sigmagr2, 
                   -0.5*mugr2/sigmagr2-log(sigmagr)-0.5*(loge2+logpi))
          abc <- abc + def
          loglik <- abc[3] + 0.5*(log(pi)-log(-abc[1])) - abc[2]^2/(4*abc[1])
        } else {
          #sigmagr = 0, i.e. gr=mugr with probability 1
          loglik <- abc[1]*mugr[1]^2 + abc[2]*mugr[1] + abc[3]
        }
      } else if(distgr[1] == 'maxlik') {
        # maxlik: derivative : 2*abc[1]*gr+abc[2] == 0
        grmax <- -0.5*abc[2]/abc[1]
        loglik <- abc[1]*grmax^2 + abc[2]*grmax + abc[3]
      } 
    }  
    if(debug) {
      debugdata <- data.table(alpha=list(alpha), theta=list(theta), sigma=list(sigma), sigmae=list(sigmae),
                              alphasigma2=list(alphasigma2), theta2=list(theta2), sigma2=list(sigma2),
                              sigmae2=list(sigmae2), logsigma=list(logsigma), logsigmae=list(logsigmae), 
                              t=list(t),
                              talpha=list(talpha), etalpha=list(etalpha), e2talpha=list(e2talpha), 
                              fe2talpha=list(fe2talpha), log_fe2talpha=list(log_fe2talpha), 
                              r0=list(r0), pif=list(pif), abc=list(abc), 
                              distgr=list(distgr), mugr=list(mugr), sigmagr=list(sigmagr),
                              def=list(c(-0.5/sigmagr^2, 
                                         mugr/sigmagr^2, -0.5*mugr^2/sigmagr^2-log(sigmagr)-0.5*(loge2+logpi))),
                              loglik=list(loglik), grmax=list(-0.5*abc[2]/abc[1]), 
                              availRmpfr=list(availRmpfr), usempfr=list(usempfr), precbits=list(precbits))
      if(!exists('.likVTreeOUDebug')) {
        .likVTreeOUDebug <<- debugdata
      } else {
        .likVTreeOUDebug <<- rbindlist(list(.likVTreeOUDebug, debugdata))
      }
    }

    if(is.double(alphaorig) & is.double(thetaorig) & is.double(sigmaorig) & 
         is.double(sigmaeorig) & is.double(vorig)) {
      loglik <- as.double(loglik)
    }
    
    value <- if(log) loglik else exp(loglik)
    if(exists('grmax'))
      attr(value, 'grmax') <- as.double(grmax)
    value
  }
}


#' Calculate alpha given heritability H2, sigma and sigmae
#' @export
alpha.poumm <- function(H2, sigma, sigmae) {
  sigma^2*(1-H2)/(2*H2*sigmae^2)
}

#' Calculate sigmaOU given heritability at time, alpha and sigmae
#' @param H2 numeric, the heritability at time t
#' @param alpha numeric, selection strength of the OU process
#' @param sigmae numeric, environmental phenotypic deviation at the tips
#' @param t numeric, time since the beginning of the OU process, default Inf
#' 
#' 
#' @export
sigmaOU.poumm <- function(H2, alpha, sigmae, t=Inf) {
  if(H2>1 | H2<0 | alpha < 0 | sigmae < 0 | t < 0) {
    stop("H2(t) should be in [0, 1], alpha, sigmae and t should be non-negative.")
  }
  if(is.infinite(t)) {
    if(alpha==0) {
      # BM
      sqrt(sigmae^2*H2/(t*(1-H2)))
    } else {
      # alpha>0, OU in equilibrium
      sqrt(2*alpha*H2*sigmae^2/(1-H2))
    }
  } else {
    # t is finite
    if(alpha==0) {
      # BM
      sqrt(sigmae^2*H2/(t*(1-H2)))
    } else {
      # alpha>0, OU not in equilibrium
      sqrt(2*alpha*H2*sigmae^2/((1-exp(-2*alpha*t))*(1-H2)))
    } 
  }
}

#' Genotypic variance at time
#' @param alpha,sigma numeric, parameters of the OU process acting on the genetic contributions 
#' @param sigmae numeric, environmental standard deviation
#' @param t time from the beginning of the process at which heritability should be calculated, i.e.
#' epidemiologic time
#' @export
sigmaG.poumm <- function(alpha, sigma, sigmae, t) {
  lenoutput <- max(sapply(list(alpha, sigma, sigmae, t), length))
  if(lenoutput>1) {
    alpha <- rep(alpha, lenoutput/length(alpha))
    sigma <- rep(sigma, lenoutput/length(sigma))
    sigmae <- rep(sigmae, lenoutput/length(sigmae))
    t <- rep(t, lenoutput/length(t))
  }
  ifelse(alpha>0, 0.5*(1-exp(-2*alpha*t))/alpha*sigma^2, t*sigma^2)  
}

#' Calculate sigmae given alpha, sigma, and H2.
#' @details Uses the formula 
#' H2=varStOU(alpha, sigma)/(varStOU(alpha, sigma)+sigmae^2)
#' @export
sigmae.poumm <- function(alpha, sigma, H2) {
  sigmaStOU2 <- varStOU(alpha, sigma)
  sqrt(sigmaStOU2/H2 - sigmaStOU2)
}

#' Broad-sense heritability estimated from the empirical variance of
#' the observed phenotypes and sigmae
#' @param z numerical vector of observed phenotypes
#' @param sigmae numerical standard deviation of the environmental deviation
#' @return numerical between 0 and 1 
#' @export
H2e.poumm <- function(z, sigmae, tree=NULL, tFrom=0, tTo=Inf) {
  if(!is.null(tree)) {
    tipTimes <- nodeTimes(tree)[1:length(tree$tip.label)]
    z <- z[which(tipTimes>=tFrom & tipTimes<=tTo)]
  }
  1-(sigmae^2/var(z))
}

#' Broad-sense heritability estimated at time t
#' @param alpha,sigma numeric, parameters of the OU process acting on the genetic contributions 
#' @param sigmae numeric, environmental standard deviation
#' @param t time from the beginning of the process at which heritability should be calculated, i.e.
#' epidemiologic time
#' @param tm average time for within host evolution from getting infected until getting measured or 
#' passing on the infection to another host
#' @export
H2.poumm <- function(alpha, sigma, sigmae, t, tm=0
                         #, as=c('vector', 'mean', 'median')
                         ) {
  lenoutput <- max(sapply(list(alpha, sigma, sigmae, t, tm), length))
  if(lenoutput>1) {
    alpha <- rep(alpha, lenoutput/length(alpha))
    sigma <- rep(sigma, lenoutput/length(sigma))
    sigmae <- rep(sigmae, lenoutput/length(sigmae))
    t <- rep(t, lenoutput/length(t))
    tm <- rep(tm, lenoutput/length(tm))
  }
  sigmag2 <- ifelse(alpha>0, 0.5*(1-exp(-2*alpha*t))/alpha*sigma^2, t*sigma^2)  
  sigmagm2 <- ifelse(alpha>0, 0.5*(1-exp(-2*alpha*tm))/alpha*sigma^2, tm*sigma^2)
  
  ifelse(t>tm, (sigmag2-sigmagm2)/(sigmag2+sigmae^2), rep(0, length(t)))
}



