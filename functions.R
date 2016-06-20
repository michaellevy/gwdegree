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