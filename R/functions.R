utils::globalVariables(".")  # For the magrittr pipe to circumvent cran-check note

#' Run GW-Degree Shiny Application
#'
#' @return NULL Called for side-effect of launching shiny application
#' @export
#'
#' @import ergm sna ggplot2 magrittr scales dplyr gtools shiny shinydashboard
#' @importFrom grDevices topo.colors
#' @importFrom graphics plot
#' @importFrom network network.size
#' @importFrom network network
#' @importFrom stats rbinom
#' @importFrom tidyr gather
#'
#' @examples
#' \dontrun{
#' gwdegree()
#' }
gwdegree <- function() {
  shiny::runApp(system.file('app.R', package='gwdegree'))
  return(NULL)
}

plotDeltaGWD = function(d, cols) {
  ggplot(d, aes_(x = ~degree, y = ~delta_GWD, color = ~as.factor(theta_s))) +
    geom_line(size = 1) +
    ylab(expression(paste(delta, ' GWD'))) +
    xlab('Node degree') +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = cols, name = expression(theta[s]),
                       guide = guide_legend(reverse = TRUE))

}

makeDeltaGWDdF = function(kmin = 0, kmax = 20, decay = 1) {
  d = expand.grid(degree = kmin:kmax, theta_s = decay)
  d$delta_GWD = deltaGWD(d$degree, d$theta_s)
  d
}

deltaGWD = function(k, theta_s)  (1 - exp(-theta_s))^k

plotNet = function(net, vCol, mtext) {
  plot(net,
       vertex.col = vCol,
       pad = 0,
       edge.col = 'gray',
       vertex.border = 'white',
       vertex.lwd = .35,
       vertex.cex = 4 - log10(network::network.size(net))
  )
  mtext(mtext, cex = 1.5)
}


simCCCent = function(gwdRange = c(-2, 2), gwespRange = c(-.5, .5),
                     theta_s = 1.5, theta_t = .5, gridSize = 3, nsim = 1,
                     netSize = 100, density = .02) {
  # Simulate GWDxGWESP-grid of networks and measure centralization and clustering
  # Returns a data.frame with one row for each cell in the grid (ie gridSize^2 rows). Four columns: gwd and gwesp parameter values from which the networks were simulated, and measured centralization and clustering coefficient (mean if nsim > 1) of each network.
  dfForSim = expand.grid(
    gwd = seq(gwdRange[1], gwdRange[2], len = gridSize),
    gwesp = seq(gwespRange[1], gwespRange[2], len = gridSize))

  N = makeNetwork(netSize, density)

  cbind(dfForSim,
        do.call(rbind,
                lapply(1:nrow(dfForSim), function(i) {
                  n = simulate.formula(N ~ gwdegree(theta_s, TRUE) + gwesp(theta_t, TRUE),
                                       coef = unlist(dfForSim[i, ]),
                                       constraints = ~ edges,
                                       nsim = nsim)
                  data.frame(
                    Centralization = mean(degCent(n)),
                    ClusteringCoef = mean(clusteringCoef(n))
                  )
                })
        )
  )
}

degCent = function(n) {
  degs = sna::degree(n)
  sum(abs(max(degs) - degs)) / sna::degree(n, tmaxdev = TRUE)
}

clusteringCoef = function(net)
  unname(3 * summary(net ~ triangles) / summary(net ~ twopath))

makeNetwork = function(nSize, nDensity) {
  # Same as ergm::as.network.numeric, but doesn't create note on R CMD check
  m = matrix(rep(0, nSize^2), nrow = nSize)
  m[upper.tri(m)] = rbinom((nSize ^ 2 - nSize)/ 2 , 1, nDensity)
  network::network(m, directed = FALSE)
}

plotHeatmaps = function(df) {
  # Plot grid heat maps of centralization and clustering. df is of the structure produced by simCCCent.
  list(
    cent =
      ggplot(df, aes_(x = ~gwd, y = ~gwesp, fill = ~Centralization)) +
      geom_raster(alpha = .7) +
      scale_fill_gradientn(colours = topo.colors(1e3), name = "") +
      xlab(expression(theta[GWD])) +
      ylab(expression(theta[GWESP])) +
      theme(aspect.ratio = 1,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()),

    cc =
      ggplot(df, aes_(x = ~gwd, y = ~gwesp, fill = ~ClusteringCoef)) +
      geom_raster(alpha = .7) +
      scale_fill_gradientn(colours = topo.colors(1e3), name = "") +
      xlab(expression(theta[GWD])) +
      ylab(expression(theta[GWESP])) +
      theme(aspect.ratio = 1,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  )
}
