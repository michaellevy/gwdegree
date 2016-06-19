lapply(c('shiny', 'shinydashboard', 'statnet', 'ggplot2', 'dplyr', 'scales'), 
       library, char = TRUE)
theme_set(theme_bw())
source("functions.R")


# Tabs ####

### Change stat  ####
gwdDecayOptions = list(3, 2, 1, .5, 0, -.5, -1)
decayColors = structure(rev(RColorBrewer::brewer.pal(length(gwdDecayOptions) + 2, 'BuPu')[-1:-2]),
                        names = unlist(gwdDecayOptions))
changeStatTab = 
    tabItem(tabName = "stat",
            h2("GWDegree Statistic Behavior"),
            "Plot the change to the GWD statistic for adding an edge to a node of given degree, at various decay-parameter values.",
            tags$br(), tags$br(),
            
            fluidRow(
                column(width = 5,
                       box(width = NULL, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                           status = "primary", solidHeader = TRUE,
                           # Slider for single decay value:
                           # sliderInput("theta_s", step = 0.1,
                           #             label = h4(HTML("&theta;<sub>S</sub> (aka \"decay\")")), 
                           #             min = -1, max = 5, value = 1.5)
                           
                           # Checkboxes for multiple decay values:
                           checkboxGroupInput("theta_s", label = NULL, # h4(HTML("&theta;<sub>S</sub> (aka \"decay\")")),
                                              choices = gwdDecayOptions,
                                              selected = gwdDecayOptions[unlist(gwdDecayOptions) >= 0])
                       ),
                       box(width = NULL, title = "Degree Range",
                           status = "primary", solidHeader = TRUE,
                           sliderInput("degreeRange", step = 1,
                                       label = NULL, 
                                       min = 0, max = 50, 
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
            "The change in the log-odds of an edge is given by the 
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
degDistTab = tabItem(tabName = "degdist",
                     h2("Network plot and degree distribution histogram here."),
                     "Some text to orient the user.
                     Simulations are costly. E.g. 50 replicates of a 500-node network
                     takes about a minute.",
                     tags$br(), tags$br(),
                     
                     # Inputs
                     fluidRow(
                         box(width = 6, title = HTML("&theta;<sub>GWD</sub>"),
                             status = "primary", solidHeader = TRUE,
                             sliderInput("gwd", label = NULL, ticks = TRUE,
                                         min = -5, max = 5, value = -2, step = 0.1)),
                         box(width = 6, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                             status = "primary", solidHeader = TRUE,
                             sliderInput("decay", label = NULL, ticks = TRUE,
                                         min = 0, max = 5, value = 1.5, step = 0.1))
                     ),
                     
                     fluidRow(
                         box(width = 4, title = "Network Size",
                             status = "primary", solidHeader = TRUE,
                             sliderInput("netSize", 
                                         label = NULL, ticks = FALSE,
                                         min = 3, max = 1000, 
                                         value = 50, step = 1)),
                         box(width = 4, title = "Network Density",
                             status = "primary", solidHeader = TRUE,
                             sliderInput("netDensity", 
                                         label = NULL, ticks = FALSE,
                                         min = 0, max = 1, 
                                         value = .05, step = .01)),
                         box(width = 4, title = "Number Simulations",
                             status = "primary", solidHeader = TRUE,
                             sliderInput("reps", 
                                         label = NULL, ticks = FALSE,
                                         min = 1, max = 500, 
                                         value = 10, step = 1))
                     ),
                     
                     fluidRow(
                         box(title = "Run New Simulations", status = "primary",
                             actionButton("simulate", label = "Go!"))
                     ),
                     
                     
                     # Outputs
                     
                     ## Sample graph (visNet?)
                     
                     ## Degree Distribution
                     fluidRow(
                         box(status = "primary", plotOutput("degDistPlot"))
                     )
)


### GWESP  ####
gwespTab = tabItem(tabName = "gwesp",
                   h2("Centralization and clustering heatmaps here.")
                   
)




# UI  ####
header = dashboardHeader(title = "GW Degree")

sidebar = 
    dashboardSidebar(sidebarMenu(
        menuItem("Statistic Behavior", tabName = "stat"),
        menuItem("Parameter & Degree Distribution", tabName = "degdist"),
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
        
        ### Change stat plot  ####
        deltaGWDdf = reactive({
            makeDeltaGWDdF(kmin = input$degreeRange[1],
                           kmax = input$degreeRange[2],
                           decay = as.numeric(input$theta_s))
        })
        output$statPlot = renderPlot({
            plotDeltaGWD(d = deltaGWDdf(), cols = decayColors)
        })
        
        
        ### Degree Distribution ####
        degDists = 
            eventReactive(input$simulate, {
                
                
                networks = vector('list', 2L)
                networks[[1]] = replicate(input$reps, makeRandomGraph(input$netSize, input$netDensity), simplify = FALSE)
                
                # Simulate networks
                # Don't care about the ERGM estimates; so turn down the control parameters
                # The edges term is there only to get the model to converge. Estimate is not used (edges constrained)
                m = ergm(networks[[1]][[1]] ~ gwdegree(0, FALSE) + edges,
                         control = control.ergm(MCMC.samplesize = 1e1, MCMLE.maxit = 1, 
                                                loglik.control = control.logLik.ergm(nsteps = 1)))
                coefs = coef(m)
                coefs[1] = input$gwd
                coefs[2] = input$decay
                networks[[2]] = simulate(m, coef = coefs, constraints = ~edges, nsim = reps
                                         # , control = control.simulate.ergm(MCMC.burnin = 1e4, 
                                         #                               MCMC.interval = 1e4)
                )
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
                filled = left_join(filled, dd, by = c("degree", "type"))
                filled$meanFreq[is.na(filled$meanFreq)] = 0
                filled
                
            })
        
        output$degDistPlot = renderPlot({
            ggplot(degDists(), aes(x = degree, y = meanFreq, fill = type)) +
                geom_bar(stat = 'identity', position = position_dodge(width = .5),
                         alpha = .75, width = 1.2, color = 'black') +
                ylab('Mean frequency') + xlab('Degree') +
                xlim(c(0, NA)) + 
                scale_fill_brewer(palette = 'BuPu', name = 'Network', 
                                  labels = c('Random Graph', 'Chosen Parameters')) +
                theme(legend.justification = c(1, 1), legend.position = c(1, 1))
            
        })            
        
    })

# Run the app  ####
shinyApp(ui = ui, server = server)

