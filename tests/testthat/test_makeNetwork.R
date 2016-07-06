library(ergm)
context("Test makeNetwork")

nNodes = sample(5:100, 1)

test_that("makeNetwork makes a network object",
          expect_is(makeNetwork(nNodes, runif(1, .01, .5)), "network"))

test_that("makeNetwork makes an empty graph if density is zero",
          expect_equal(network.edgecount(makeNetwork(nNodes, 0)), 0))

test_that("makeNetwork makes a complete graph if density is one",
          expect_equal(network.edgecount(makeNetwork(nNodes, 1)),
                       (nNodes ^ 2 - nNodes) / 2))
