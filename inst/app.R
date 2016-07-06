utils::globalVariables(".")  # For the magrittr pipe to circumvent cran-check note

# colors
gwdDecayOptions = list(3, 2, 1, .5, 0, -.5, -1)
cols = structure(rev(RColorBrewer::brewer.pal(length(gwdDecayOptions) + 2, 'BuPu')[-1:-2]),
                 names = unlist(gwdDecayOptions))

# Tabs ####

### Change stat  ####
changeStatTab =
  tabItem(tabName = "stat",
          h2("GW Degree Statistic Behavior"),
          HTML("Plot the change to the GWD statistic for adding a half-edge to a node of given degree,
            at various decay-parameter (&theta;<sub>S</sub>) values.<br><br>
               The GWD statistic is given by:"),
          withMathJax("$$GWD(\\mathbf{y}, \\theta_s) = e^{\\theta_s} \\sum_{k=1}^{n-1}[1-(1-e^{-\\theta_s})^{k}]D_{k}(\\mathbf{y})$$"),
          withMathJax("Where \\(\\mathbf{y}\\) is a network of \\(n\\) nodes, \\(D_{k}\\) of which have degree-\\(k\\), and \\(\\theta_s\\)
                      is a decay parameter that controls the severity of geometric weighting.
                      Noting that adding a half-edge to a node of degree-k increments \\(D_{k+1}\\) and decrements \\(D_{k}\\), the change statistic is:"),
          withMathJax("$$\\delta GWD = (1 - e^{-\\theta_s})^k$$"),
          withMathJax("Note that \\(\\theta_s < 0\\) makes GWD hard to interpret, and when \\(\\theta_s < -ln(2)\\),
                      GWD becomes poorly behaved."),
            tags$br(), tags$br(),

          fluidRow(
            column(width = 5,
                   box(width = NULL, title = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                       status = "primary", solidHeader = TRUE,
                       checkboxGroupInput("theta_s", label = NULL, # h4(HTML("&theta;<sub>S</sub> (aka \"decay\")")),
                                          choices = gwdDecayOptions,
                                          selected = gwdDecayOptions[unlist(gwdDecayOptions) >= 0])
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
          )
  )


### Degree distribution  ####
degDistTab =
  tabItem(
    tabName = "degdist",
    h2("Parameter & Degree Distribution"),
    HTML("The two networks plotted below are the same size and density. The purple
             one is a random network (all edges equiprobable),
            and the blue one is simulated with the
             selected GWD parameter and decay values. Adjust the sliders and click
            \"Simulate!\" to see how the blue
             network changes relative to the baseline of the purple one.<br><br>
             The boxplots at right show degree distributions from multiple
             networks simulated the same way as the plotted networks. You can examine
            degree centralization directly as a response to parameter values in the
            \"Centralization, Clustering, & GWESP\" tab at left."),
    tags$br(), tags$br(),

    # Inputs
    fluidRow(
      box(width = 4, title = "Network Parameters",
          status = "primary", solidHeader = TRUE,
          sliderInput("netSize",
                      label = "Number of Nodes", ticks = FALSE,
                      min = 20, max = 1e3,
                      value = 50, step = 5),
          sliderInput("density",
                      label = "Network Density", ticks = FALSE,
                      min = .001, max = .2,
                      value = .08,
                      step = .001)
      ),

      box(width = 4, title = "GWDegree Parameters",
          status = "primary", solidHeader = TRUE,
          sliderInput("gwd", label = HTML("&theta;<sub>GWD</sub>"),
                      ticks = FALSE,
                      min = -5, max = 5, value = -3, step = 0.1),
          sliderInput("decay", label = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                      ticks = FALSE,
                      min = 0, max = 5, value = 2, step = 0.1)
      ),

      box(width = 4, title = "Simulation Parameters",
          status = "primary", solidHeader = TRUE,
          sliderInput("reps",
                      label = "Number of Simulated Networks", ticks = FALSE,
                      min = 3, max = 50,
                      value = 20, step = 1),
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
        be concentrated on popular actors) and clustering (the tendency for open two-paths to become closed triangles)
        respond across a range of both GWD- and GWESP-parameter values, their respective decay parameters, and network size and density.
        <br><br>
        Note that these simulations are computationally expensive. Larger and denser networks or simulations on larger grids may take minutes for each replicate network."),
    tags$br(), tags$br(),

    # Inputs
    fluidRow(

      box(width = 3, title = "Network Parameters",
          status = "primary", solidHeader = TRUE,
          sliderInput("netSize3",
                      label = "Number of Nodes", ticks = FALSE,
                      min = 20, max = 1e3,
                      value = 50, step = 5),
          sliderInput("density3",
                      label = "Network Density", ticks = FALSE,
                      min = .001, max = .2,
                      value = .08,
                      step = .001)
      ),

      box(width = 3, title = "GWDegree Parameters",
          status = "primary", solidHeader = TRUE,
          sliderInput("gwd3", label = HTML("&theta;<sub>GWD</sub>"),
                      ticks = FALSE,
                      min = -5, max = 5, value = c(-3, 3), step = 0.1),
          sliderInput("theta_s3", label = HTML("&theta;<sub>S</sub> (aka \"decay\")"),
                      ticks = FALSE,
                      min = 0, max = 5, value = 2, step = 0.1)
      ),

      box(width = 3, title = "GWESP Parameters",
          status = "primary", solidHeader = TRUE,
          sliderInput("gwesp3", label = HTML("&theta;<sub>GWESP</sub>"),
                      ticks = FALSE,
                      min = -2, max = 2, value = c(-1, 1), step = 0.1),
          sliderInput("theta_t3", label = HTML("&theta;<sub>T</sub> (aka \"alpha\")"),
                      ticks = FALSE,
                      min = 0, max = 1, value = .25, step = 0.05)
      ),

      box(width = 3, title = "Simulation Parameters",
          status = "primary", solidHeader = TRUE,
          selectInput("reps3",
                      label = "Number of Replicate Networks",
                      choices = as.list(1:3), selected = 3),
          selectInput("gridSize",
                      label = "Grid size",
                      choices = list("5 x 5" = 5
                                     , "9 x 9" = 9
                                     , "17 x 17" = 17
                                     , "31 x 31" = 31
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
    menuItem("Parameter & Degree Distribution", tabName = "degdist"),
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

    # lapply(c('ergm', 'network', 'sna', 'ggplot2', 'dplyr', 'scales', 'magrittr'),
    #        require, char = TRUE)
    # source("../R/functions.R")
    theme_set(theme_bw(base_size = 20))


    ### Change-statistic plot  ####
    deltaGWDdf = reactive({
      gwdegree:::makeDeltaGWDdF(kmin = input$degreeRange[1],
                                kmax = input$degreeRange[2],
                                decay = as.numeric(input$theta_s))
    })
    output$statPlot = renderPlot({
      gwdegree:::plotDeltaGWD(d = deltaGWDdf(), cols = cols)
    })


    ### Degree Distribution plots ####

    #### Update the density slider to avoid having too many edges:
    observe({
      updateSliderInput(session, "density",
                        value = 4 / input$netSize,
                        min = .001,
                        max = round(10 / input$netSize, 2)
      )
    })

    netCols = structure(cols[c(2, 5)], names = c("random", "gwd"))

    #### Simulate many and plot degree distributions
    degDists =
      reactive({

        input$goDD  # This button causes the reactive to update, but
        # not all the isolated code that directly takes slider inputs below.

        isolate({
          networks = vector('list', 2L)
          networks[[1]] = replicate(input$reps,
                                    gwdegree:::makeNetwork(input$netSize, density = input$density),
                                    simplify = FALSE)
          networks[[2]] =
            lapply(networks[[1]], function(basis) {  # To keep number edges same across types
              if (!network.edgecount(basis)) return(basis)
              simulate.formula(basis ~ gwdegree(input$decay, TRUE), coef = input$gwd,
                               constraints = ~ edges)
            })
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
              dplyr::mutate(.,
                     sim = rep(1:input$reps, length(unique(.$degree))),
                     degree = as.integer(degree)) %>%
              dplyr::arrange(degree)

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
        scale_x_discrete(name = "Degree",
                         breaks = seq(0, max(degDists()[["dd"]]$degree), by = 5)) +
        scale_fill_manual(values = netCols, name = 'Network',
                          labels = c('Random Network', 'GWD Network')) +
        theme(legend.justification = c(1, 1), legend.position = c(1, 1))
    })

    #### Plot sample networks ####
    output$bothGraphs = renderPlot({

      par(mfrow = c(1, 2), mar = c(0, 0, 3.5, 0))
      gwdegree:::plotNet(degDists()[["random"]], netCols["random"], "Random Network\n(Null Model)")
      gwdegree:::plotNet(degDists()[["gwd"]], netCols["gwd"], "GWD Network\n")
    })


    ### GWD & GWESP ####

    ##### Reactive density slider:
    observe({
      updateSliderInput(session, "density3",
                        value = 4 / input$netSize3,
                        min = .001,
                        max = round(10 / input$netSize3, 2)
      )
    })

    heatmaps = reactive({
      input$goGWESP
      isolate({
        gwdegree:::simCCCent(gwdRange = input$gwd3,
                             gwespRange = input$gwesp3,
                             theta_s = input$theta_s3,
                             theta_t = input$theta_t3,
                             gridSize = as.integer(input$gridSize),
                             nsim = as.integer(input$reps3),
                             netSize = input$netSize3,
                             density = input$density3) %>%
          gwdegree:::plotHeatmaps()
      })
    })
    output$centHeatmap = renderPlot(heatmaps()$cent)
    output$ccHeatmap = renderPlot(heatmaps()$cc)


  })

# Run the app  ####
shinyApp(ui = ui, server = server)
