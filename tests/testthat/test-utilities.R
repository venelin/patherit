library(testthat)
library(patherit)

context("Test utilities")

m <- matrix(rnorm(20), 4, 5)

cumul <- function(m) {
  for(i in 2:ncol(m)) {
    m[, i] = m[, i-1]+m[, i]
  }
  m
}

print(m)
m2 <- cumul(m)
m3 <- cumulativeRows(m)

all.equal(m3, m)
all.equal(m3, m2)
print(m2)
print(m3)
test_that(all.equal(cumul(m) == cumulativeRows(m)))
