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
