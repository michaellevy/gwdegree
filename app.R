require(shiny)
require(shinydashboard)

# colors
gwdDecayOptions = list(3, 2, 1, .5, 0, -.5, -1)
cols = structure(rev(RColorBrewer::brewer.pal(length(gwdDecayOptions) + 2, 'BuPu')[-1:-2]),
                 names = unlist(gwdDecayOptions))

# Tabs ####

### Change stat  ####
changeStatTab = 
    tabItem(tabName = "stat",
            h2("GW Degree Statistic Behavior"),
            HTML("Plot the change to the GWD statistic for adding an edge to a node of given degree, 
            at various decay-parameter (&theta;<sub>S</sub>) values.<br><br>
            Three take-home points:"),
            tags$ol(
                tags$li(HTML("The statistic increases with any edge, and it always increases more for edge added to a 
                        lower-degree node than a higher-degree node. The implication is that for positive values of the corresponding 
                        parameter (&theta;<sub>GWD</sub>), the likelihood of an edge on a lower-degree node is greater than an edge on a high-degree node.
                        This is discussed in more detail below the charts.")), 
                tags$li(HTML("&theta;<sub>S</sub> should never be less than zero.")), 
                tags$li(HTML("GWD always weights changes to lower-degree nodes more than changes to higher-degree nodes.
                        The closer &theta;<sub>S</sub> is to zero, the more dramatic the difference. Thus, if you want
                        &theta;<sub>GWD</sub> to capture a popularity effect that operates <em>among popular actors</em>, choose a large
                        value for &theta;<sub>S</sub>. The particular value will depend on the size and density of the network -- the
                        \"Parameter and Degree Distribution\" tab is meant to facilitate this choice."))
            ),
            tags$br(), tags$br(),
            
            fluidRow(
                column(width = 5,
                       box(width = NULL, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                           status = "primary", solidHeader = TRUE,
                           checkboxGroupInput("theta_s", label = NULL, # h4(HTML("&theta;<sub>S</sub> (aka \"decay\")")),
                                              choices = gwdDecayOptions,
                                              selected = 1)  # gwdDecayOptions[unlist(gwdDecayOptions) >= 0])
                       ),
                       box(width = NULL, title = "Degree Range",
                           status = "primary", solidHeader = TRUE,
                           sliderInput("degreeRange", step = 1,
                                       label = NULL, ticks = FALSE,
                                       min = 0, max = 100, 
                                       value = c(0, 20))
                       )
                       
                ),
                
                column(width = 7,
                       box(width = NULL, 
                           status = "primary", solidHeader = TRUE,
                           plotOutput("statPlot")
                       )
                )
            ),
            tags$br(), tags$br(),
            h3("Analytical Solution"),
            "The GWD statistic is given by:",
            withMathJax("$$GWD(\\mathbf{y}, \\theta_s) = e^{\\theta_s} \\sum_{k=1}^{n-1}[1-(1-e^{-\\theta_s})^{k}]D_{k}(\\mathbf{y})$$"),
            withMathJax("Where \\(\\mathbf{y}\\) is the network of \\(n\\) nodes, \\(D_{k}\\) of which have degree-\\(k\\), and \\(\\theta_s\\)
                    is a decay parameter that controls the severity of geometric weighting. 
                    The change statistic for adding a half-edge to a node of degree-k, then, is:"),
            withMathJax("$$\\delta GWD = (1 - e^{-\\theta_s})^k$$"),
            HTML("From the change statistic, we see why &theta;<sub>S</sub> < 0 makes interpretation of &theta;<sub>GWD</sub>
                 difficult: the change statistic is positive or negative depending on whether the degree of the node
                 receiving an edge has even- or odd-degree. Furthermore, we see that &theta;<sub>S</sub> < -ln(2) is deeply problematic:
                 the exponentiated term is less than -1, so as k increases, the statistic explodes."),
            
            h3("Implication"),
            HTML("In an ERGM, the change in the log-odds of an edge is given by the 
            product of the change statistic and the parameter value (akin to a logistic regression). 
            So, for positive GWD parameters,
            edges that increase the change statistic more are more likely; conversely,
            for negative GWD parameter values, edges that
            increase the change statistic less are more likely. Note in the plot above that edges to nodes of low degree increase the statistic more.
            This implies that positive GWD parameters reward edges to low-degree nodes, and negative GWD parameters reward edges to higher-degree nodes.
            Therefore, <b>negative GWD parameter estimates are consistent with networks more centralized (i.e. with higher-variance degree distributions)
            than expected by chance, <i>ceteris paribus.</i></b> <br><br>
            You may notice that for &theta;<sub>GWD</sub> < 0, all edges are 
            penalized and for &theta;<sub>GWD</sub> > 0, all edges are rewarded. This fact is offset by the inclusion of a 
            density or edges term in the model which establishes the baseline likelihood of any edge, or fixing the number
            of edges in the network: It is the <i>relative</i> likelihood of edges on various nodes that is of interest.")
    )


### Degree distribution  ####
degDistTab = 
    tabItem(
        tabName = "degdist",
        h2("Parameter & Degree Distribution"),
        HTML("The two networks plotted below have the same number of nodes and edges,
             but the purple one is a random network (all edges equiprobable), 
            while the blue one is simulated with the 
             selected GWD parameter and decay values. Adjust the sliders and click
            \"Simulate!\" to see how the blue
             network changes relative to the baseline of the purple one.<br><br>
             The histograms at right are average degree distributions from multiple 
             networks simulated the same way as the plotted networks. You can examine
            degree centralization directly as a response to parameter values in the
            \"Centralization, Clustering, & GWESP\" tab at left.<br><br>
             Note that when &theta;<sub>GWD</sub> < 0, there are more high- and
             low-degree nodes, and when &theta;<sub>GWD</sub> > 0, there are more
             middle-degree nodes. That is, negative values of &theta;<sub>GWD</sub>
             are consistent with a network more centralized than expected by chance. 
             The particular shape of the degree distribution
             is a function of &theta;<sub>S</sub> (aka \"decay\"), as well as network
             size and density. As &theta;<sub>S</sub> approaches zero, changes
             to low-degree nodes have more of an effect on the statistic than changes
             to high-degree nodes. See the \"Statistic Behavior\" tab for the intuition on this.
             Increasing the number of simulated networks reduces 
             noise from stochastic simulation, but simulations can be computationally
             expensive, especially for large and dense networks."),
        tags$br(), tags$br(),
        
        # Inputs
        fluidRow(
            box(width = 4, title = "Network Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("netSize", 
                            label = "Number of Nodes", ticks = FALSE,
                            min = 5, max = 500, 
                            value = 100, step = 5),
                sliderInput("meanDegree", 
                            label = "Average Degree", ticks = FALSE,
                            min = .04, max = 8, 
                            value = 2, 
                            step = .1)
            ),
            
            box(width = 4, title = "GWDegree Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("gwd", label = HTML("&theta;<sub>GWD</sub>"), 
                            ticks = FALSE,
                            min = -5, max = 5, value = -3, step = 0.1),
                # box(width = 6, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                # status = "primary", solidHeader = TRUE,
                sliderInput("decay", label = HTML("&theta;<sub>S</sub> (aka \"decay\")"), 
                            ticks = FALSE,
                            min = 0, max = 5, value = 2, step = 0.1)
            ),
            
            box(width = 4, title = "Simulation Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("reps", 
                            label = "Number of Simulated Networks", ticks = FALSE,
                            min = 10, max = 30, 
                            value = 15, step = 1),
                actionButton("goDD", label = strong("Simulate!"))
            )
        ),
        
        
        # Outputs
        
        fluidRow(
            ## Sample graphs (visNet?)
            box(status = "primary", plotOutput("bothGraphs")),
            
            ## Degree Distribution
            box(status = "primary", plotOutput("degDistPlot"))
        )
    )


### Confoundedness  ####
gwespTab = 
    tabItem(
        tabName = "gwesp",
        h2("Confoundedness of Centralization & Clustering"),
        HTML("Here you can examine how network centralization (the tendency for edges to
        be concentrated on popular actors) and clustering (the tendency to form complete triangles)
        respond across a range of not only GWD-, but also GWESP-parameter values, and
        their respective decay parameters, as well as network size and density. To examine
        just &theta;<sub>GWD</sub>'s effects, look across the horizontal 
        &theta;<sub>GWESP</sub> = 0 line. Note, as in the other tabs, that greater
        values of &theta;<sub>GWD</sub> correspond to <i>less</i> centralized networks.<br><br>
        Note also the confoundedness of clustering and centralization: Negative values
        of &theta;<sub>GWD</sub> create centralized networks, where there are many
        edges in close proximity, which leads to 
        complete triangles by chance. Likewise, positive &theta;<sub>GWESP</sub> values
        create complete triangles, which puts edges near each other, incidentally causing
        centralized networks. Thankfully, ERGMs that include
        both GWD and GWESP parameters can distinguish the two effects to the extent
        that they leave a distinctive structural signature in the network. 
        For this reason, modeling
        efforts that include one GW-term should at least begin with both terms (and 
        probably GWDSP as well to distinguish tendency for two-paths from the tendency
        for two-paths that become closed triangles.)
        <br><br>
        Note that these simulations are computationally expensive. Especially for larger
             and denser networks, the calculations will several dozen seconds. If you would
             like to run these simulations for larger networks, consider cloning the app
             and deploying it locally. Instructions and code are "),
        tags$a(href="https://github.com/michaellevy/GWDegree-Shiny", "on GitHub"), ".",
        
        
        tags$br(), tags$br(),
        
        # Inputs
        fluidRow(
            
            box(width = 3, title = "Network Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("netSize3", 
                            label = "Number of Nodes", ticks = FALSE,
                            min = 5, max = 250, 
                            value = 50, step = 5),
                sliderInput("meanDegree3", 
                            label = "Average Degree", ticks = FALSE,
                            min = .04, max = 8, 
                            value = 2, step = .1)
                
            ),

            box(width = 3, title = "GWDegree Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("gwd3", label = HTML("&theta;<sub>GWD</sub>"), 
                            ticks = FALSE,
                            min = -5, max = 5, value = c(-3, 3), step = 0.1),
                # box(width = 6, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                # status = "primary", solidHeader = TRUE,
                sliderInput("theta_s3", label = HTML("&theta;<sub>S</sub> (aka \"decay\")"), 
                            ticks = FALSE,
                            min = 0, max = 5, value = 2, step = 0.1)
            ),
            
            box(width = 3, title = "GWESP Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("gwesp3", label = HTML("&theta;<sub>GWESP</sub>"), 
                            ticks = FALSE,
                            min = -2, max = 2, value = c(-1, 1), step = 0.1),
                # box(width = 6, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                # status = "primary", solidHeader = TRUE,
                sliderInput("theta_t3", label = HTML("&theta;<sub>T</sub> (aka \"alpha\")"), 
                            ticks = FALSE,
                            min = 0, max = 1, value = .25, step = 0.05)
            ),
            
            box(width = 3, title = "Simulation Parameters",
                status = "primary", solidHeader = TRUE,
                selectInput("reps3", 
                            label = "Number of Simulated Networks", 
                            choices = as.list(1:3), selected = 3),
                selectInput("gridSize",
                            label = "Grid size",
                            choices = list("5 x 5" = 5
                                           , "9 x 9" = 9
                                           # , "17 x 17" = 17
                                           # , "31 x 31" = 31
                                           ),
                            selected = 5),
                actionButton("goGWESP", label = strong("Simulate!"))
            )
        ),
        
        fluidRow(
            box(status = "primary", solidHeader = TRUE, 
                plotOutput("centHeatmap"), 
                title = "Degree Centralization"),
            box(status = "primary", solidHeader = TRUE,
                plotOutput("ccHeatmap"),
                title = "Clustering Coefficient")
        )
        
        
        
    )




# UI  ####
header = dashboardHeader(title = "GW Degree")

sidebar = 
    dashboardSidebar(sidebarMenu(
        menuItem("Statistic Behavior", tabName = "stat"),
        menuItem("Parameter & Degree Distribution", tabName = "degdist", selected = TRUE),
        menuItem("Centralization, Clustering, & GWESP", tabName = "gwesp")
    ))


body =
    dashboardBody(
        tabItems(
            changeStatTab,
            degDistTab,
            gwespTab
        )
    )

ui = dashboardPage(header, sidebar, body)

# Server  ####
server =
    shinyServer(function(input, output, session) {
        
        lapply(c('ergm', 'network', 'sna', 'ggplot2', 'dplyr', 'scales', 'magrittr'), 
               require, char = TRUE)
        theme_set(theme_bw())
        source("functions.R")
        
        
        ### Change-statistic plot  ####
        deltaGWDdf = reactive({
            makeDeltaGWDdF(kmin = input$degreeRange[1],
                           kmax = input$degreeRange[2],
                           decay = as.numeric(input$theta_s))
        })
        output$statPlot = renderPlot({
            plotDeltaGWD(d = deltaGWDdf(), cols = cols)
        })
        
        
        ### Degree Distribution plots ####
        
        #### Update the meanDegree slider to avoid having too many edges:
        observe({
            updateSliderInput(session, "meanDegree", 
                              value = log10(input$netSize),
                              min = 4 / input$netSize,   # Need at least one edge
                              max = min(10, input$netSize / 20 + 3, (input$netSize - 2) / 2)
                              )  
            # (N - 1) / 2 is actual max for undir'd net, but it trips up simulate() 
            # and nearly-complete graphs are expensive and not of much interest.
        })

        netCols = structure(cols[c(2, 5)], names = c("random", "gwd"))
        
        #### Simulate many and plot degree distributions
        degDists =
            reactive({
                
                input$goDD  # This button causes the reactive to update, but
                            # not all the isolated code that directly takes
                            # slider inputs below.

                
                isolate({
                    networks = vector('list', 2L)
                    networks[[1]] = replicate(input$reps, 
                                              intergraph::asNetwork(igraph::erdos.renyi.game(input$netSize, input$netSize * input$meanDegree / 2, "gnm")), 
                                              simplify = FALSE)
                    
                    networks[[2]] = simulate.formula(networks[[1]][[1]] ~ gwdegree(input$decay, TRUE), coef = input$gwd,
                                                     constraints = ~ edges,
                                                     nsim = input$reps)
                    names(networks) = c('random', 'gwd')
                    # Calculate the degrees of each graph and make a freq table of each
                    # The factor levels is to include 0-counts
                    # Then take the median count of each degree
                    
                    dd = lapply(seq_along(networks), function(type) {
                        # type is stack of simulated graphs of a particular type
                        # Split into each simulation, tabulate degrees, and use `smartbind` to match on names (ie, degree)
                        deg = 
                            degree(networks[[type]], g = 1:length(networks[[type]]), gmode = 'graph') %>% 
                            split(., 1:ncol(.)) %>%
                            lapply(table) %>%
                            do.call(gtools::smartbind, .) %>% 
                            tidyr::gather(degree, count) %>%
                            # Need number of sims to join with correct number of 0s for degrees that never occur
                            mutate(., 
                                   sim = rep(1:input$reps, length(unique(.$degree))),
                                   degree = as.integer(degree)) %>% 
                            arrange(degree)
                        
                        # Where smartbind fills with NA, there were no nodes of that degree, so replace with zero:
                        deg[is.na(deg)] = 0
                        deg$type = names(networks[type])
                        
                        list(deg = deg, 
                             var = data.frame(
                                 type = names(networks)[type],
                                 var = mean(sapply(networks[[type]], function(n) var(degree(n, gmode = 'graph'))))
                                 ))
                    }) 
                    
                    vari = do.call(rbind, lapply(dd, "[[", "var"))
                    deg = do.call(rbind, lapply(dd, "[[", "deg"))
                    
                    # Fill in the degrees-never-found to keep box-widths constant
                    filled = expand.grid(
                        degree = 0:max(deg$degree),
                        type = names(networks),
                        sim = 1:input$reps,
                        stringsAsFactors = FALSE
                    )
                    filled = dplyr::left_join(filled, deg, by = c("degree", "type", "sim"))
                    filled$count[is.na(filled$count)] = 0
                    filled$type = factor(filled$type, levels = sort(unique(filled$type), decreasing = TRUE))
                })
                
                list(dd = filled, variance = vari, random = networks[[1]][[1]], gwd = networks[[2]][[1]])
                
                
            })
        
        #### Plot degree distributions ####
        
        # Pretty the x-axis labels and get outlier points colored
        # And add average variance
        output$degDistPlot = renderPlot({
            
            ggplot(degDists()[["dd"]], aes(x = factor(degree), y = count, fill = type)) +
                geom_boxplot(position = position_dodge(width = .5),
                             alpha = .75, color = 'black') +
                ylab('Node Count') + 
                # xlab('Degree') +
                scale_x_discrete(name = "Degree",
                                 breaks = seq(0, max(degDists()[["dd"]]$degree), by = 5)) +
                scale_fill_manual(values = netCols, name = 'Network',
                                  labels = c('Random Network', 'GWD Network')) +
                theme(legend.justification = c(1, 1), legend.position = c(1, 1))
        })            
        
        #### Plot sample networks ####
        output$bothGraphs = renderPlot({
            
            par(mfrow = c(1, 2), mar = c(0, 0, 3.5, 0))
            plotNet(degDists()[["random"]], netCols["random"], "Random Network\n(Null Model)")
            plotNet(degDists()[["gwd"]], netCols["gwd"], "GWD Network\n")
        })
        
        
        ### GWD & GWESP ####
        
        ##### Reactive degree slider:
        #### Update the meanDegree slider to avoid having too many edges, as in deg-dist tab
        observe({
            updateSliderInput(session, "meanDegree3", 
                              value = log10(input$netSize),
                              min = 4 / input$netSize,   # Need at least one edge
                              max = min(10, input$netSize / 20 + 3, (input$netSize - 2) / 2)
            )
        })
        
        heatmaps = reactive({
            input$goGWESP
            isolate({
                simCCCent(gwdRange = input$gwd3,
                          gwespRange = input$gwesp3,
                          theta_s = input$theta_s3,
                          theta_t = input$theta_t3,
                          gridSize = as.integer(input$gridSize),
                          nsim = as.integer(input$reps3),
                          netSize = input$netSize3,
                          meanDegree = input$meanDegree3)
            })
        })
        output$centHeatmap = renderPlot(heatmaps()$cent)
        output$ccHeatmap = renderPlot(heatmaps()$cc)
        
        
    })

# Run the app  ####
shinyApp(ui = ui, server = server)

