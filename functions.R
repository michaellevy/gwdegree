plotDeltaGWD = function(d, cols) {
    ggplot(d, aes(x = degree, y = delta_GWD, color = as.factor(theta_s))) + 
        geom_line(size = 1) +
        ylab(expression(paste(delta, ' GWD'))) +
        xlab('Node degree') +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme_bw(base_size = 20) +
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

makeRandomGraph = function(nodes, density, directed = FALSE)
{
    require(statnet)
    
    edges = density * nodes * (nodes - 1)
    if(!directed) 
        edges = edges / 2
    
    maxEdges = nodes * (nodes - 1) / abs(directed - 2)
    if(edges >= maxEdges)
        stop("Too many edges. I won't make you a complete graph!")
    n = network.initialize(nodes, directed = directed)
    ed = as.data.frame(t(replicate(edges, sort(sample(network.size(n), 2)))))
    while(any(duplicated(ed))) {
        for(i in which(duplicated(ed)))
            ed[i, ] = sort(sample(network.size(n), 2))
    }
    add.edges(n, ed$V1, ed$V2)  
    n
}

plotNet = function(net, vCol, mtext) {
    plot(net, 
         vertex.col = vCol, 
         pad = 0,
         edge.col = 'gray', 
         vertex.border = 'white',
         vertex.lwd = .35,
         vertex.cex = 4 - log10(network.size(net))
         ) 
    mtext(mtext)
}

simCCCent = function(gwdRange = c(-2, 2), gwespRange = c(-.5, .5),
                     theta_s = 1.5, theta_t = .5, gridSize = 3, nsim = 1,
                     netSize = 100, meanDegree = 3) {
    
    dfForSim = expand.grid(
            gwd = seq(gwdRange[1], gwdRange[2], len = gridSize),
            gwesp = seq(gwespRange[1], gwespRange[2], len = gridSize))
    
    N = intergraph::asNetwork(igraph::erdos.renyi.game(netSize, netSize * meanDegree, "gnm"))
    
    dfForSim = 
        # Could parallelize this! Not sure how many cores shinyapps.io gives...xs
        lapply(1:nrow(dfForSim), function(i) {
            n = simulate.formula(N ~ gwdegree(theta_s, TRUE) + gwesp(theta_t, TRUE), 
                                 coef = unlist(dfForSim[i, ]),
                                 constraints = ~ edges,
                                 nsim = nsim)
            data.frame(
                Centralization = mean(centralization(n, 'degree', mode = 'graph')),
                ClusteringCoef = mean(clusteringCoef(n))
            )
        }) %>%
        do.call(rbind, .) %>%
        cbind(dfForSim, .)
    
    plotHeatmaps(dfForSim)
    
}

clusteringCoef = function(net)
    unname(3 * summary(net ~ triangles) / summary(net ~ twopath))

plotHeatmaps = function(df) {
    list(
        cent = 
            ggplot(df, aes(x = gwd, y = gwesp, fill = Centralization)) +
            geom_raster(alpha = .7) +
            scale_fill_gradientn(colours = topo.colors(1e3), name = 'Degree\ncentralization')+
            # scale_fill_gradient2(name = 'Degree\ncentralization', mid = 'lightgray', 
            #                      midpoint = strSims$Centralization[strSims$GWD == 0 & strSims$GWESP == 0]) +
            xlab(expression(theta[GWD])) +
            ylab(expression(theta[GWESP])) +
            theme(aspect.ratio = 1),
        
        cc = 
            ggplot(df, aes(x = gwd, y = gwesp, fill = ClusteringCoef)) +
            geom_raster(alpha = .7) +
            scale_fill_gradientn(colours = topo.colors(1e3), name = 'Clustering\nCoefficient')+
            # scale_fill_gradient2(name = 'Clustering\nCoefficient', mid = 'lightgray', 
            #                      midpoint = strSims$ClusteringCoef[strSims$GWD == 0 & strSims$GWESP == 0]) +
            xlab(expression(theta[GWD])) +
            ylab(expression(theta[GWESP])) +
            theme(aspect.ratio = 1)
    )
# cowplot::plot_grid(cent, cc, labels = letters[1:2], align = 'vh')
}
