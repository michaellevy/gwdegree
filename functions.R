plotDeltaGWD = function(d, cols) {
    p = ggplot(d, aes(x = degree, y = delta_GWD, color = as.factor(theta_s))) + 
        geom_line(size = 1) +
        ylab(expression(paste(delta, ' GWD'))) +
        # xlab('Degree of node receiving half-edge') +
        xlab('Node degree') +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme_bw(base_size = 20) +
        theme(aspect.ratio = 1) +
        scale_color_manual(values = cols, name = expression(theta[s]),
                           guide = guide_legend(reverse = TRUE))
    #     
    # if(length(unique(d$theta_s)) == 1)
    #     p = p + scale_color_manual(values = cols, guide = 'none') else
    #         p = p + scale_color_manual(values = cols, name = expression(theta[s]),
    #                                    guide = guide_legend(reverse = TRUE))
        p
        
}

makeDeltaGWDdF = function(kmin = 0, kmax = 20, decay = 1) {
    d = expand.grid(degree = kmin:kmax, theta_s = decay)
    d$delta_GWD = deltaGWD(d$degree, d$theta_s)
    d
}

deltaGWD = function(k, theta_s)  (1 - exp(-theta_s))^k
