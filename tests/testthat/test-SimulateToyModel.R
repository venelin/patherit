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

library(testthat)
library(patherit)

context("Setup toy model")

# Here we setup the toy-model, specifying the pathogen strains (here, called genotypes), 
# the general host-type x carried-strain effects (GEVs), 
# the SIR parameters (birth/death/recovery/transmission/mutation rate), 
# the total population sizes and the time for simulation


# make results reproducible (using same seed as for the discrete generation case)
set.seed(5)

numsAllelesAtSites <- c(4, 4, 4, 4, 4)
# genotype encodings with their allele contents
genotypes <- generateGenotypes(numsAllelesAtSites)

# number of host-types
n <- 24
pe <- runif(n)
pe <- pe/sum(pe)

# General host-type x strain effects (expected phenotypes for genotype by environment combinations)
GEVs <- matrix(NA, nrow=n, ncol=nrow(genotypes))

for(g in 1:nrow(genotypes)) {
  GEVs[, g] <- rnorm(n=n, mean=2+2/nrow(genotypes)*(g), sd=0.4)
}

matplot(GEVs, type='l', lty=1, col=1:6)

sigmae <- .6

timeStep <- 0.05

pg.init <- rep(0, nrow(genotypes))
pg.init[1] <- 1

sde <- rep(sigmae, n)

# death-rate as a function of viral load and natural death rate mu
rateDie <- function(z, mu) {
  V <- 10^z
  Dmin <- 2
  Dmax <- 25*12
  D50 <- 10^3
  Dk <- 1.4
  (V^Dk+D50^Dk)/(Dmin*(V^Dk+D50^Dk)+((Dmax-Dmin)*D50^Dk)) + mu
}

meanRateDie <- function(z, mu) {
  rep(0.01, length(z))
}

# infectionrate as a function of viral load and rate of risky contacts
rateInfect <- function(z, rateContact) {
  V <- 10^z
  Emin <- .3
  Emax <- .6
  E50 <- 10^3
  Ek <- 1.4
  E <- Emin+(Emax-Emin)*V^Ek/(V^Ek+E50^Ek)
  E*rateContact
}

meanRateInfect <- function(z, rateContact) {
  rep(0.45*rateContact, length(z))
}

# per locus mutation rates
rateMutate <- function(GEValues, es, envs, genes) {
  z <- es+GEValues[cbind(envs, genes)]
  V <- 10^z
  Mmin <- 0.00
  Mmax <- 0.2
  M50 <- 10^3
  Mk <- 1.4
  Mmin+(Mmax-Mmin)*V^Mk/(V^Mk+M50^Mk)
}

meanRateMutate <- function(GEValues, es, envs, genes) {
  rep(0.01, length(genes))
}

rateSampleNeutral <- 1/(4*12)

# maximum time before starting graceful fadeout of the epidemic
maxTime <- 2400

# continue the epidemic outbreak until reaching ...
maxNTips <- 10000

# is the special environmental effect unique for each pathogen genoetype in an individual
eUniqForEachG <- TRUE

# are only beneficial (i.e. increasing the trait-value) mutations allowed
selectWithinHost <- TRUE

# time to continue the simulation of transmission after reaching maxNTips (this was introduced
# in order to study post-outbreak dynamics, i.e. epidemic waves after exhaustion of the 
# susceptible pool)
expandTimeAfterMaxNTips <- 4

#initial population size (equilibrium in the absence of disease)
N <- 1e5

# setting the between-host and within-host dynamics of the simulation:
# natural death rate
mu <- 1/850

# constant birth rate that maintains this equilibrium
nu <- ifelse(is.finite(N), mu*N, 0)

rateContact <- 1

rateSample <- rateSampleNeutral

selectWithinHost <- TRUE

epidemic <- simulateEpidemic(
  Ninit=N, nu=nu, mu=mu, pe=pe, sde=sde, pg.init=pg.init, GEValues=GEVs, 
  rateContact=1/6, rateInfect=rateInfect, rateDie=rateDie, 
  rateSample=rateSampleNeutral, 
  rateMutate=meanRateMutate, 
  numsAllelesAtSites=numsAllelesAtSites, eUniqForEachG=eUniqForEachG, 
  selectWithinHost=TRUE,
  timeStep=timeStep, maxTime=maxTime, maxNTips=maxNTips, 
  expandTimeAfterMaxNTips=expandTimeAfterMaxNTips, 
  process="select/select")



