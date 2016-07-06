library(ergm)
context("clustering coefficient function")

nNodes = sample(5:100, 1)

test_that("A complete graph of random size has a clustering coefficient of 1",
          expect_equal(gwdegree:::clusteringCoef(network(nNodes, density = 1, directed = FALSE)), 1))

test_that("A graph with a single two-path has a clustering coefficient of 0",
          expect_equal(gwdegree:::clusteringCoef(network(3, numedges = 2, directed = FALSE)), 0))

test_that("Clustering coefficient can't be calculated for an empty graph of random size",
          expect_true(is.nan(gwdegree:::clusteringCoef(network(nNodes, density = 0, directed = FALSE)))))


