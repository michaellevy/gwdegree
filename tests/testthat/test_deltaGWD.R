library(ergm)
context("Testing analytical delta-gwdegree solution")

test_that("Deleting an edge: deltaGWD() matches simulation-based calculation for a random dyad in random network with a random decay-parameter value in [0, 3].",
          {
            decay = runif(1, 0, 3)
            repeat ({  # Make sure we have an edge that is present
              n = network(sample(5:50, 1), density = runif(1, .01, .5), directed = FALSE)
              nodesToAdd = sample(network.size(n), 2, replace = FALSE)
              if (n[nodesToAdd[1], nodesToAdd[2]])
                break
            })
            nodesDegrees = degree(n, nodes = nodesToAdd, gmode = "graph")
            analyticalResult = -sum(gwdegree:::deltaGWD(nodesDegrees - 1, decay))

            originalGWD = summary(n ~ gwdegree(decay, fixed = TRUE))
            delete.edges(n, eid = get.edgeIDs(n, nodesToAdd[1], nodesToAdd[2]))  # modifies `n` in place
            newGWD = summary(n ~ gwdegree(decay, fixed = TRUE))
            simResult = unname(newGWD - originalGWD)

            expect_equal(analyticalResult, simResult)
          })

test_that("Adding an edge: deltaGWD() matches simulation-based calculation for a random dyad in random network with a random decay-parameter value in [0, 3].",
          {
            decay = runif(1, 0, 3)
            repeat ({  # Make sure we do not have an edge that is present
              n = network(sample(5:50, 1), density = runif(1, .01, .5), directed = FALSE)
              nodesToAdd = sample(network.size(n), 2, replace = FALSE)
              if (!n[nodesToAdd[1], nodesToAdd[2]])
                break
            })
            nodesDegrees = degree(n, nodes = nodesToAdd, gmode = "graph")
            analyticalResult = sum(gwdegree:::deltaGWD(nodesDegrees, decay))

            originalGWD = summary(n ~ gwdegree(decay, fixed = TRUE))
            add.edge(n, nodesToAdd[1], nodesToAdd[2])  # modifies `n` in place
            newGWD = summary(n ~ gwdegree(decay, fixed = TRUE))
            simResult = unname(newGWD - originalGWD)

            expect_equal(analyticalResult, simResult)
          })
