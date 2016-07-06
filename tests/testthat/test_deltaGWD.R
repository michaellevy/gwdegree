library(ergm)
context("delta gwdegree")

test_that("deltaGWD(), which is based on my analytical solution, produces the same change statistics as ergm::summary() for two random nodes from a random network with a random decay-parameter value in [0, 3].",
          {
            decay = runif(1, 0, 3)
            n = network(sample(2:50, 1), density = runif(1, .01, .5), directed = FALSE)
            nodesToAdd = sample(network.size(n), 2, replace = FALSE)
            nodesDegrees = degree(n, nodes = nodesToAdd, gmode = "graph")
            analyticalResult = sum(gwdegree:::deltaGWD(nodesDegrees, decay))

            originalGWD = summary(n ~ gwdegree(decay, fixed = TRUE))
            add.edge(n, nodesToAdd[1], nodesToAdd[2])  # modifies `n` in place
            newGWD = summary(n ~ gwdegree(decay, fixed = TRUE))
            simResult = unname(newGWD - originalGWD)

            expect_equal(analyticalResult, simResult)
          })
