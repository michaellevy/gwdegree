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
            "Plot the change to the GWD statistic for adding an edge to a node of given degree, at various decay-parameter values.",
            tags$br(), tags$br(),
            "Note that for any positive value of ", HTML("&theta;<sub>S</sub>"),
            " the statistic is increasing for any edge, and an edge on a lower-degree 
            node increases the statistic more than an edge on a higher-degree node.
            The implications of this are discussed below.",
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
                    The change statistic for adding a half-edge to a node of degree-k is:"),
            withMathJax("$$\\delta GWD = (1 - e^{-\\theta_s})^k$$"),
            
            h3("Implication"),
            "In an ERGM, the change in the log-odds of an edge is given by the 
            product of the change statistic and the parameter value. 
            So for positive GWD parameters,
            edges that increase the change statistic more are more likely; conversely,
            for negative GWD parameter values, edges that
            increase the change statistic less are more likely. Note in the plot above that edges to nodes of low degree increase the statistic more.
            This implies that positive GWD parameters reward edges to low-degree nodes, and negative GWD parameters reward edges to higher-degree nodes.
            Therefore, ", strong("negative GWD parameter estimates are consistent with networks more centralized (i.e. with higher-variance degree distributions)
            than expected by chance, ", em("ceteris paribus."))
    )


### Degree distribution  ####
degDistTab = 
    tabItem(
        tabName = "degdist",
        h2("Parameter & Degree Distribution"),
        HTML("Simulate networks with various values of
        the GWD-parameter and GWD-decay-parameter and examine the affects on 
        network structure (left) and degree distribution (right), compared to
        random graphs of the same size and density. <br><br>
        The degree distribution histograms are averaged over the chosen number of 
        simulated networks; the networks at left are sample realizations from
        those simulations.
        &theta;<sub>GWD</sub> is an ERGM parameter associated with the GW-Degree statistic.
        &theta;<sub>S</sub> (aka \"decay\") is a shape parameter that controls 
        how severely edges after the first are discounted. See the \"Statistic Behavior\"
        tab for the intuition on how &theta;<sub>S</sub> works."),
        tags$br(), tags$br(),
        
        # Inputs
        fluidRow(
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
            
            box(width = 4, title = "Network Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("netSize", 
                            label = "Number of Nodes", ticks = FALSE,
                            min = 5, max = 500, 
                            value = 100, step = 5),
                sliderInput("meanDegree", 
                            label = "Average Degree", ticks = FALSE,
                            min = 0, max = 10, 
                            value = 3, step = .1)
            ),
            
            box(width = 4, title = "Simulation Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("reps", 
                            label = "Number of Simulated Networks", ticks = FALSE,
                            min = 10, max = 100, 
                            value = 30, step = 1))
        ),
        
        
        # Outputs
        
        fluidRow(
            ## Sample graphs (visNet?)
            box(status = "primary", plotOutput("bothGraphs")),

            ## Degree Distribution
            box(status = "primary", plotOutput("degDistPlot"))
        )
    )


### GWESP  ####
gwespTab = 
    tabItem(
        tabName = "gwesp",
        h2("Confoundedness of Centralization & Clustering"),
        "The two affect each other ... ",
        
        # Inputs
        fluidRow(
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
            
            box(width = 3, title = "Network Parameters",
                status = "primary", solidHeader = TRUE,
                sliderInput("netSize3", 
                            label = "Number of Nodes", ticks = FALSE,
                            min = 5, max = 500, 
                            value = 50, step = 5),
                sliderInput("meanDegree3", 
                            label = "Average Degree", ticks = FALSE,
                            min = 0, max = 10, 
                            value = 3, step = .1)
            ),
            
            box(width = 3, title = "Simulation Parameters",
                status = "primary", solidHeader = TRUE,
                selectInput("reps3", 
                            label = "Number of Simulated Networks", 
                            choices = as.list(1:10), selected = 3),
                selectInput("gridSize",
                            label = "Grid size",
                            choices = list("5 x 5" = 5,
                                           "9 x 9" = 9,
                                           "17 x 17" = 17,
                                           "31 x 31" = 31),
                            selected = 5)
                )
        ),
        
                   
                   fluidRow(
                       box(status = "primary", plotOutput("centHeatmap")),
                       box(status = "primary", plotOutput("ccHeatmap"))
                   )
                   
                   
                   
)




# UI  ####
header = dashboardHeader(title = "GW Degree")

sidebar = 
    dashboardSidebar(sidebarMenu(
        menuItem("Statistic Behavior", tabName = "stat"),
        menuItem("Parameter & Degree Distribution", tabName = "degdist", selected = TRUE),
        menuItem("Confounded with GWESP", tabName = "gwesp")
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
    shinyServer(function(input, output) {
        
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
        netCols = structure(cols[c(2, 5)], names = c("random", "gwd"))

        #### Simulate many and plot degree distributions
        degDists = 
            reactive({
                networks = vector('list', 2L)
                networks[[1]] = replicate(input$reps, 
                                          intergraph::asNetwork(igraph::erdos.renyi.game(input$netSize, input$netSize * input$meanDegree, "gnm")), 
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
                        do.call(gtools::smartbind, .)
                    # Where smartbind fills with NA, there were no nodes of that degree, so replace with zero:
                    deg[is.na(deg)] = 0
                    
                    # Calculate median count for each degree and mean prob of each degree and organize in df
                    data.frame(degree = as.integer(colnames(deg)),
                               meanFreq = apply(deg, MARGIN = 2, mean),
                               type = names(networks[type]))
                    
                }) %>%
                    do.call(rbind, .)
                
                # Fill in the degrees-never-found to keep box-widths constant
                filled = expand.grid(
                    degree = 0:max(dd$degree),
                    type = names(networks)
                )
                filled = dplyr::left_join(filled, dd, by = c("degree", "type"))
                filled$meanFreq[is.na(filled$meanFreq)] = 0
                
                list(dd = filled, random = networks[[1]][[1]], gwd = networks[[2]][[1]])
                
            })
        
        #### Plot degree distributions ####
        output$degDistPlot = renderPlot({
            ggplot(degDists()[["dd"]], aes(x = degree, y = meanFreq, fill = type)) +
                geom_bar(stat = 'identity', position = position_dodge(width = .5),
                         alpha = .75, width = 1.2, color = 'black') +
                ylab('Mean frequency') + xlab('Degree') +
                xlim(c(0, NA)) + 
                scale_fill_manual(values = netCols, name = 'Network', 
                                  labels = c('Random Graph', 'Chosen Parameters')) +
                theme(legend.justification = c(1, 1), legend.position = c(1, 1))
        })            
        
        #### Plot sample networks ####
        output$bothGraphs = renderPlot({
            par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
            plotNet(degDists()[["random"]], netCols["random"], "Random Graph")
            plotNet(degDists()[["gwd"]], netCols["gwd"], "Chosen Parameters")
        })
        
        
        ### GWD & GWESP ####
        heatmaps = reactive({
            simCCCent(gwdRange = input$gwd3,
                      gwespRange = input$gwesp3,
                      theta_s = input$theta_s3,
                      theta_t = input$theta_t3,
                      gridSize = as.integer(input$gridSize),
                      nsim = as.integer(input$reps3),
                      netSize = input$netSize3,
                      meanDegree = input$meanDegree3)
        })
        output$centHeatmap = renderPlot(heatmaps()$cent)
        output$ccHeatmap = renderPlot(heatmaps()$cc)
        
        
    })

# Run the app  ####
shinyApp(ui = ui, server = server)

