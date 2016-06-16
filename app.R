lapply(c('shiny', 'shinydashboard', 'statnet', 'ggplot2', 'dplyr', 'scales'), 
       library, char = TRUE)
theme_set(theme_bw())
source("functions.R")

# UI  ####
header = dashboardHeader(title = "GW Degree")

sidebar = 
    dashboardSidebar(sidebarMenu(
        menuItem("Statistic Behavior", tabName = "stat"),
        menuItem("Parameter & Degree Distribution", tabName = "degdist"),
        menuItem("Confounded with GWESP", tabName = "gwesp")
    ))

gwdDecayOptions = list(3, 2, 1, .5, 0, -.5, -1)
decayColors = structure(RColorBrewer::brewer.pal(length(gwdDecayOptions), 'BuPu'),
                        names = unlist(gwdDecayOptions))
body =
    dashboardBody(
        tabItems(
            
            ### Change stat  ####
            tabItem(tabName = "stat",
                    h2("GWDegree Statistic Behavior"),
                    "Plot the change to the GWD statistic for adding an 
                    edge to a node of given degree, at various decay parameter values.",
                    tags$br(), tags$br(),
                    "ANALYTICAL SOLUTION HERE?",
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
                                                      selected = 1) # gwdDecayOptions[unlist(gwdDecayOptions) >= 0])
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
                    "The change in the log-odds of an edge is given by the 
                    product of the change statistic and the parameter value. 
                    So for positive GWD parameters,
                    edges that increase the change statistic more are more likely,
                    and conversely for negative GWD parameter values, edges that
                    increase the change statistic less are more likely."
            ),
            
            ### Degree distribution  ####
            tabItem(tabName = "degdist",
                    h2("Network plot and degree distribution histogram here.")
            ),
            
            ### GWESP  ####
            tabItem(tabName = "gwesp",
                    h2("Centralization and clustering heatmaps here.")
            )
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
    })

# Run the app  ####
shinyApp(ui = ui, server = server)

