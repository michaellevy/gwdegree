library(sna)
context("clustering coefficient function")

networks = list(
  nEmpty = gwdegree:::makeNetwork(100, 0),
  nFull = gwdegree:::makeNetwork(100, 1),
  nOther = gwdegree:::makeNetwork(100, .2)
)

test_that("Degree centralization is the same for my function as package:sna",
          expect_equal(
            sapply(networks, gwdegree:::degCent),
            sapply(networks, sna::centralization, sna::degree, "graph")
          ))

