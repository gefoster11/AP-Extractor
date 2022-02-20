#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# This application is meant to load an action potential locations file following extraction with matlab with columns in this order ID, condition, time, sample_loc
# This application is mean to load a signals file with columns ID, condition, time, pressure, ecg, MSNA_raw, MSNA_int following processing with AP_File_combine.R
# The detected action potentials can be viewed and accepted or rejected
# Once finished click the download button to save a .zip file containing the outputted action potential files

source("./functions/helper_functions.R")

# Sets max file size to 500 Mb
options(shiny.maxRequestSize = 500*1024^2, scipen = 999)

# Package Dependency
packages = c("remotes",
             "shiny",
             "shinyBS",
             "tidyverse",
             "MESS",
             "zip",
             "thematic",
             "shinythemes",
             "signal",
             "features")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Package Dependency - Install plotly from github
   if (!require("plotly", character.only = TRUE)) {
      remotes::install_github("ropensci/plotly")
      library("plotly", character.only = TRUE)
    }

# ---- Define UI for data upload app ----
ui <- fluidPage(
    theme = shinytheme("united"),
    # ---- App title ----
    titlePanel("Action Potential Extractor"),
    
    # ---- Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # ---- Sidebar panel for inputs ----
        sidebarPanel(width = 3,
            
            # ---- Input: Select Raw MSNA file ----
            fileInput("Signals_10KHz", "Raw MSNA at 10 KHz",
                      accept = ".txt"), #Choose .txt file with Time in column one 
            #and the raw MSNA signal with appropriate bandpass filter in column two. Data at 10 KHz


            # ---- Input: Select AP Locations file ----
            fileInput("Locations", "AP Locations",
                      accept = ".txt"),

            # ---- Input: MSNA Signals @ 25 Hz ----
            fileInput("Signals_25KHz", "Neurogram File",
                      accept = ".csv"),
                        
            # ---- Input: MSNA burst/beat file ----
            fileInput("MSNAapp_file", "MSNA Burst/Beat File",
                      accept = ".csv"),
            
            tags$hr(),
            
            # ---- Input: Analysis conditions ----
            selectizeInput('conditions', "Choose Analysis Conditions", choices = "", multiple = TRUE, width = "90%"),
            bsTooltip("conditions", "Choose Conditions to Analyze", 
                      "right", options = list(container = "body")),
            
            # ---- Input: Analysis conditions ----
            selectizeInput('baseline_condition', "Which condition is baseline?", choices = "", width = "90%"),
            bsTooltip("baseline_condition", "Choose baseline condition", 
                      "right", options = list(container = "body")),
        
            # ---- Output sample frequency ----
            textOutput("fs"),
            
            # ---- Horizontal line ----
            tags$hr(),
            
            actionButton("reset", "Reset", width = "90%"),
            bsTooltip("reset", "click to reset selected bursts",
                      "right", options = list(container = "body")),
            
            # ---- Horizontal line ----
            tags$hr(),

            downloadButton("save", "Download"),
            bsTooltip("save", "saves the selected data to a csv file", "right", 
                      options = list(container = "body")),
            
            # ---- Horizontal line ----
            tags$hr(),
            
            # ---- Reference ----
            tags$div(
              HTML("<p>2022; Created by Glen Foster</p><br>")),

        ),
        
        # ---- Main panel ----
        mainPanel(
            
            # ---- Output: Tabset w/ parsing table, selected data, and plot ----
            tabsetPanel(type = "tabs",
                        tabPanel("Data", 
                                 dataTableOutput("data")),
                        tabPanel("Action Potential Selection",

                                 fluidRow(
                                   plotlyOutput("ap_characteristics", height = "200px")
                                 ),
                                 # ---- Horizontal line ----
                                 tags$hr(),
                                 
                                 fluidRow(
                                   # ---- Input: RangeSlider Scale ----
                                   column(width = 3, align = "center",
                                          actionButton("toggle", "Remove AP", width = "75%"),
                                          bsTooltip("toggle", "Removed selected AP",
                                                    "right", options = list(container = "body")),
                                          tags$hr()
                                   ),
                                   
                                 ),
                                 # ---- Main Plot ----
                                 fluidRow(
                                   plotlyOutput("ap_plot", height = "600px")  
                                 )
                                 
                                 ),
                        )
            )
    )
)

# ---- Define server logic ----
server <- function(input, output, session) {
    
    # ---- initialize reactive values ----
    values <- reactiveValues(
        fs = NULL, # sample frequency
        df = NULL, # data
        plot.click.results = NULL, # store plot clicks
        )
    
    # ---- Listen: delimiter inputs ----
    loadlisten <- reactive({
        list(input$Signals_10KHz$datapath, input$Locations$datapath)
    })
    
    # ---- Observe: data load - df ----
    observeEvent(loadlisten(), {
        
        req(input$Signals_10KHz$datapath, input$Locations$datapath)
      
      #browser()  
      
        tryCatch( {
            # load locations
            locations <- read_tsv(file = input$Locations$datapath) 
            
            # load signals
            signals <- read_tsv(file = input$Signals_10KHz$datapath)
            
            # determine and set fs
            fs <- round(1/(signals$time[2] - signals$time[1]), 1)
            values$fs <- fs

            withProgress(message = "Processing Data - be patient!", {
            df <- AP_extract(locations, signals) %>%
              mutate(ap_keep = TRUE)
            
            # determine clusters
            breaks <- find_cluster_bins(df) # Helper function
            
            h <- df$amplitude[df$include == TRUE] %>% 
              hist(breaks = breaks, plot = FALSE)  
            
            df <- df %>% 
              ungroup() %>% 
              mutate(breaks = cut(.$amplitude, breaks = h$breaks), 
                     cluster = cut(.$amplitude, breaks = h$breaks, labels = FALSE, ordered_result = TRUE)) %>%
              group_by(cluster) %>% 
              nest()
            
            df$mid <- h$mids[df$cluster]
            df$lower_bound <- h$breaks[df$cluster]
            
            df <- df %>% ungroup() %>% unnest(data)
            values$df <- df
            
            myChoices <- unique(df$condition)
            updateSelectizeInput(session, "conditions", choices = myChoices)
            updateSelectizeInput(session, "baseline_condition", choices = myChoices)
            
            }
            )
        }, error = function(e) {
            stop(safeError(e))
        }
        )
    })

    # ---- Output: sample frequency ----
    output$fs <- renderText({
        fs <- values$fs
        print(paste0("Sample Frequency = ", fs))
    })
    
    
    # ---- Output: Data Table ----
    output$data <- renderDataTable({
      req(values$df, input$conditions)
        
        return(values$df %>% 
                 plotly::filter(condition == input$conditions) %>%
                 select(!data)
               )
      
    }, options = list(pageLength = 10))
    
    
    # ---- Re-cluster data when conditions change ----
    observe({
      
      #### Unsure if this part is working yet!
      
      req(values$df, input$conditions, input$baseline_condition)
      
      df <- values$df %>% 
        select(!c(mid, lower_bound))
        
      breaks <- find_cluster_bins(df) # Helper function
      
      h <- df$amplitude[df$include == TRUE & df$ap_keep == TRUE] %>% 
        hist(breaks = breaks, plot = FALSE)  
      
      df <- df %>% 
        ungroup() %>% 
        mutate(breaks = cut(.$amplitude, breaks = h$breaks), 
               cluster = cut(.$amplitude, breaks = h$breaks, labels = FALSE, ordered_result = TRUE)) %>%
        group_by(cluster) %>% 
        nest()
      
        df$mid <- h$mids[df$cluster]
        df$lower_bound <- h$breaks[df$cluster]
      
      values$df <- df %>% ungroup() %>% unnest(data)
      
    })
    

    # ---- Output: Create AP Plot by amplitude cluster ----
    output$ap_characteristics <- renderPlotly({
      req(values$df)
      
      df <- values$df %>% 
        filter(include == TRUE & ap_keep == TRUE)
      
      h1 <- df$inter_spike_int[df$inter_spike_int < 600] %>% 
        hist(breaks = "Scott", plot = FALSE)
      h1$probs <- h1$counts/sum(h1$counts)*100
      p1 <- plot_ly(x = h1$mids, y = h1$probs) %>% 
        add_bars(name = "Inter Spike Interval") %>% 
        layout(xaxis = list(title = "Inter-Spike Interval (ms)"), yaxis = list(title = "Probability (%)"))
      h2 <- df$amplitude %>% 
        hist(breaks = "Scott", plot = FALSE)  
      h2$probs <- h2$counts/sum(h2$counts)*100  
      p2 <- plot_ly(x = h2$mids, y = h2$probs) %>% 
        add_bars(name = "Amplitude") %>% 
        layout(xaxis = list(title = "Amplitude (V)"), yaxis = list(title = "Probability (%)"))
      
      plot <- subplot(p1, p2, titleY = TRUE, titleX = TRUE, shareY = TRUE, margin = 0.07) %>% 
        layout(showlegend = FALSE)

      plot
    })

    
    # ---- Output: Create AP Plot ----
    output$ap_plot <- renderPlotly({
      req(values$df)
        #browser() 
      df <- values$df
      
      #browser()
      
      temp <- df %>% 
        unnest(data) %>%
        filter(include == TRUE & ap_keep == TRUE)
      
      df_key <- highlight_key(temp %>% group_by(AP_no), ~AP_no)
      
      # Plot overlayed APs
      FigA <- plot_ly(source = "A") %>%
        add_lines(data = df_key,
                  color = ~factor(breaks),
                  colors = c("blue", "red"),
                  x = ~sample,
                  y = ~AP_v,
                  line = list(width = 2),
                  legendgroup = ~ap_keep,
                  text = ~paste("AP no: ",
                                round(AP_no, 2), '<br>time:',
                                round(time, 2)
                  ),
                  customdata = ~AP_no,
        ) %>% layout(xaxis = list(title = "", showticklabels = FALSE))
      
      
      # Calculate Average AP
      # Mean AP by condition and cluster
      mean <- temp %>%
        group_by(ID, sample) %>% 
        summarise(n = n(), mean = mean(AP_v), sd = sd(AP_v), se = sd/sqrt(n)) %>% 
        mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)
      
      # Plot Mean AP by cluster by condition
      FigB <- plot_ly(source = "A") %>%
        add_lines(data = mean,
                  x = ~sample,
                  y = ~mean,
                  name = "Mean",
                  showlegend = FALSE,
                  line = list(color = "black")
        ) %>%
        add_ribbons(data = mean,
                    x = ~sample,
                    ymin = ~lower95,
                    ymax = ~upper95,
                    name = "95%CI",
                    showlegend = FALSE,
                    line = list(color = "grey", dash = 'dot', width = 1),
                    fillcolor = 'rgba(7, 164, 181, 0.2)'
        ) %>% 
        layout(xaxis = list(title = "", showticklabels = FALSE))
      
      Fig <-subplot(FigA, FigB, shareY = TRUE) %>%
        layout(title = paste0(mean$n[1], " Action Potentials Detected"),
               yaxis = list(title = "MSNA (mV)", fixedrange = FALSE)
               ) 
      
      #### Plot mean plots by cluster and plot below All APs and Mean AP.
      mean_cluster <- temp %>%
        group_by(ID, breaks, sample) %>% 
        summarise(n = n(), mean = mean(AP_v), sd = sd(AP_v), se = sd/sqrt(n)) %>% 
        mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)
      
      FigC <- mean_cluster %>% ungroup() %>%
        group_by(breaks) %>%
        group_map(~ plot_ly(data=., 
                            x = ~sample, 
                            y = ~mean, 
                            showlegend = FALSE,
                            #color = ~breaks, 
                            type = "scatter", mode="lines") %>%
                    add_ribbons(  x = ~sample,
                                  ymin = ~lower95,
                                  ymax = ~upper95,
                                  name = "95%CI",
                                  showlegend = FALSE,
                                  line = list(color = "grey", dash = 'dot', width = 1),
                                  fillcolor = 'rgba(7, 164, 181, 0.2)'
                    ) %>%
                    layout(yaxis = list(title = "MSNA (mV)", fixedrange = FALSE),
                           xaxis = list(title = "", showticklabels = FALSE))
        ) %>%
        subplot(nrows = 1, shareX = TRUE, shareY=TRUE)
      
      subplot(Fig, FigC, nrows = 2) %>% 
        highlight(on = "plotly_click", off = "plotly_relayout", debounce = 250)
      
      
    })
  
    # # Click Register
    # observeEvent(event_data("plotly_click", source = "A"), { 
    #   click_results <- values$plot.click.results
    #   
    #   event <- event_data("plotly_click", source = "A")
    #   
    #   browser()
    #   
    #   if (event$customdata %in% click_results$customdata) {
    #     click_results <- click_results %>% filter(customdata != event$customdata)
    #   } else {
    #     click_results <- click_results %>% rbind(., event)  
    #   }
    #   
    # values$plot.click.results <- click_results
    #   
    # })
    # 
    # # Clear Click Register
    # observeEvent(event_data("plotly_doubleclick", source = "A"), { 
    #   
    #   values$plot.click.results <- NULL
    # })

    # ---- Observe: plot click ----
    observeEvent(input$toggle, {
      
      eventData <- event_data("plotly_click", source = "A")
      
      if ("customdata" %in% names(eventData)) {
        
        
        df <- isolate(values$df)
        ap_no <- eventData$customdata
        
        df$ap_keep[df$AP_no == ap_no] <- !df$ap_keep[df$AP_no == ap_no]
        
        values$df <- df
        
      }
      
    })
    

    # ---- Observe: reset ----
    observeEvent(input$reset, {
        df <- values$df
        df$ap_keep <- TRUE
        values$df <- df
    })
    
    
    # ---- Downloadable csv of selected dataset ----
    output$save <- downloadHandler(

        filename = function() {
            paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-ape.zip", sep = "")
        },
        content = function(file){


          #go to a temp dir to avoid permission issues
          owd <- setwd(tempdir())
          on.exit(setwd(owd))

            fileName <- c(paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-all_aps.csv", sep = ""),
                          paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-mean_ap.csv", sep = ""),
                          paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-mean_cluster_ap.csv", sep = ""))

            df <- values$df %>% unnest(data) %>% filter(include == TRUE & ap_keep == TRUE)
            mean <- df %>%
              group_by(ID, condition, sample) %>% 
              summarise(n = n(), mean = mean(AP_v), sd = sd(AP_v), se = sd/sqrt(n)) %>% 
              mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)
            
            mean_cluster <- df %>%
              group_by(ID, condition, cluster, sample) %>% 
              summarise(n = n(), mean = mean(AP_v), sd = sd(AP_v), se = sd/sqrt(n)) %>% 
              mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)
            
              write_csv(df, fileName[[1]])  
              write_csv(mean, fileName[[2]])
              write_csv(mean_cluster, fileName[[3]])
          
          #create the zip file
          zip(file,fileName)

        }
    )
    
  
}

# Create Shiny app ----
thematic::thematic_shiny()
shinyApp(ui, server)
