#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# This application is meant to load an action potential locations file following extraction with matlab with columns in this order ID, condition, time, sample_loc
# This application is mean to load a signals file with columns ID, condition, time, pressure, ecg, MSNA_raw, MSNA_int following processing with AP_File_combine.R
# The detected action potentials can be viewed and accepted or rejected
# Once finished click the download button to save a .zip file containing the outputted action potential files

source("./functions/helper_functions.R")

# Sets max file size to 800 Mb
options(shiny.maxRequestSize = 800*1024^2,
        shiny.launch.browser = .rs.invokeShinyWindowExternal,
        scipen = 999)

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
             "features",
             "kableExtra",
             "DT",
             "ggpubr",
             "stringr"
             )

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
   # theme = shinytheme("united"),
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
            selectizeInput('base', "Choose Baseline Conditions", choices = "", multiple = FALSE, width = "90%"),
            bsTooltip("base", "Choose baseline condition", 
                      "right", options = list(container = "body")),
            
            # ---- Input: Analysis conditions ----
            selectizeInput('conditions', "Choose Analysis Conditions", choices = "", multiple = TRUE, width = "90%"),
            bsTooltip("conditions", "Choose Conditions to Analyze", 
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
              HTML("<p>2023; Created by Glen Foster</p><br>")),

        ),
        
        # ---- Main panel ----
        mainPanel(
            
            # ---- Output: Tabset w/ parsing table, selected data, and plot ----
            tabsetPanel(type = "tabs",
                        tabPanel("Data", 
                                 fluidRow(
                                   tags$div(
                                     HTML("<p>Click on table rows to toggle action potential 
                                          keep status (i.e. True/False)</p><br>")),
                        
                                 ),
                                 fluidRow(
                                   DT::dataTableOutput("data")  
                                 ),
                                 fluidRow(
                                   DT::dataTableOutput("absolute_summary")  
                                 ),
                                 fluidRow(
                                   DT::dataTableOutput("nCluster_summary")  
                                 ),
                                 ),
                        tabPanel("Action Potential Selection",

                                 fluidRow(
                                   plotlyOutput("ap_characteristics", height = "200px")
                                 ),
                                 # ---- Horizontal line ----
                                 tags$hr(),
                                 
                                 # ---- Input: Choose Cluster ----
                                 selectizeInput('Cluster', "Choose Cluster", choices = "", width = "20%"),
                                 bsTooltip("Cluster", "Choose cluster to plot", 
                                           "right", options = list(container = "body")),
                                 
                                 tags$hr(),
                                 
                                 # ---- Main Plot ----
                                 fluidRow(
                                   plotlyOutput("ap_plot", height = "600px")
                                 )
                                 
                                 ),
                        tabPanel("AP Characteristics & Signal Check",
                                 # ---- Main Plot ----
                                 fluidRow(
                                   #### SNR Check ####
                                   tableOutput("SNRCheck")
                                 ),
                                 fluidRow(
                                   #### AP Characteristics Plot ####
                                   plotlyOutput("ap_characteristics2", height = "600px")
                                 ),
                                 
                                 ),
                        
                        tabPanel("Time Series",
                                 selectizeInput('conditions2', "Choose Conditions to Plot", choices = "", multiple = FALSE, width = "20%"),
                                 bsTooltip("conditions2", "Choose Conditions to Plot", 
                                           "right", options = list(container = "body")),
                                 # ---- Main Plot ----
                                 fluidRow(
                                   #### time series plot ####
                                   plotlyOutput("TimeSeries", height = "600px")
                                 )),
                        
                        tabPanel("Correlograms",
                                 # ---- Main Plot ----
                                 fluidRow(
                                   #### Correlograms ####
                                   plotOutput("correlograms", height = "600px")
                                 )),

                        )
            )
    )
)

# ---- Define server logic ----
server <- function(input, output, session) {
  
    # ---- initialize reactive values ----
    values <- reactiveValues(
        # fs = NULL, # sample frequency
        # df = NULL, # data
        # beat = NULL, #beat/burst data
        # signals25 = NULL, #neurogram
        )
    
    # ---- Listen: AP data inputs ----
    loadlisten <- reactive({
        list(input$Signals_10KHz$datapath, input$Locations$datapath)
    })
    
    # ---- Listen: beat/burst data inputs ----
    loadlisten_beat <- reactive({
      list(input$Signals_25KHz$datapath, input$MSNAapp_file$datapath)
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
            
            withProgress(message = "Extracting APs", {
            df <- AP_extract(locations, signals) %>%
              mutate(ap_keep = TRUE)
            
            #browser()
            
            myChoices <- unique(df$condition)
            updateSelectizeInput(session, "base", choices = myChoices, selected = myChoices[[1]])
            
            # determine clusters
            breaks <- find_cluster_bins2(x = df, base_cond = myChoices[[1]]) # Helper function
            
            h <- df$ap_amplitude[df$ap_include == TRUE] %>% 
              hist(breaks = breaks, plot = FALSE)  
            
            df <- df %>% 
              ungroup() %>% 
              mutate(breaks = cut(.$ap_amplitude, breaks = h$breaks), 
                     cluster = cut(.$ap_amplitude, breaks = h$breaks, labels = FALSE, ordered_result = TRUE)) %>%
              group_by(cluster) %>% 
              nest()
            
            df$mid <- h$mids[df$cluster]
            df$lower_bound <- h$breaks[df$cluster]
            
            #browser()
            
            df <- df %>% 
              ungroup() %>% 
              unnest(data) %>%
              relocate(ID, condition) %>%
              arrange(ap_no)
            
            df <- df %>%
              mutate(ncluster = (cluster/max(cluster, na.rm = TRUE) * 10) %>% ceiling()) %>%
              relocate(ncluster, .after = cluster)
            
            values$df <- df
            values$signals10khz <- signals
            values$times <- signals %>% group_by(ID, condition) %>% summarise(t_min = min(time), t_max = max(time))
  
            
            myClusters <- c("ALL", unique(df$cluster) %>% sort(.))
            
            
            updateSelectizeInput(session, "conditions", choices = myChoices, selected = myChoices)
            updateSelectizeInput(session, "Cluster", choices = myClusters, selected = "ALL")
            updateSelectizeInput(session, "conditions2", choices = myChoices, selected = "")
            
            }
            )
        }, error = function(e) {
            stop(safeError(e))
        }
        )
    })

    # ---- Observe: data load - df ----
    observeEvent(loadlisten_beat(), {

      req(input$Signals_10KHz$datapath,
          input$Locations$datapath,
          input$Signals_25KHz$datapath,
          input$MSNAapp_file$datapath,
          values$df,
          input$conditions)

      tryCatch( {
        # load locations
        
        beat <- read_csv(file = input$MSNAapp_file$datapath)

        # load signals
        signals25 <- read_csv(file = input$Signals_25KHz$datapath)

        values$beat <- beat
        values$signals25 <- signals25
        
        withProgress(message = "Aligning APs to beats", {
          times <- values$times
          ap_beat <- values$beat
          APs <- values$df
          signal <- values$signals10khz
          
          df <- link_AP(times, ap_beat, APs, signal)
          values$df <- df %>%
            relocate(ID, condition) %>%
            arrange(beat_no, ap_no)
        })

      }, error = function(e) {stop(safeError(e))}
      )
    })
    
    # ---- Observe: summarize data - summary_df ----
    observe(priority = -1, {
      
      req(input$Signals_10KHz$datapath,
          input$Locations$datapath,
          input$Signals_25KHz$datapath,
          input$MSNAapp_file$datapath,
          values$df,
          input$conditions)

      df <- values$df
      
      if("beat_no" %in% colnames(df)) {

              summary_df <- Absolute_Summary(df)
              summary_df <- nCluster_Summary(summary_df %>% unnest(data))

              values$summary_df <- summary_df
      }

     

    })
    
    
    # ---- Output: sample frequency ----
    output$fs <- renderText({
        fs <- values$fs
        print(paste0("Sample Frequency = ", fs))
    })
    
    
    # ---- Output: Data Table ----
    output$data <- DT::renderDataTable({
      req(values$df, input$conditions)

        df <- values$df

        return(df %>% 
                 plotly::filter(condition %in% input$conditions) %>%
                 select(!data)
               )
      
    }, options = list(pageLength = 10), filter = "top", caption = "All Data")
    
    output$SNRCheck <- function() {
      req(values$df, input$conditions, values$beat, values$signals25)
      
      df <- values$df %>% 
        plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>%
        plotly::filter(condition %in% input$conditions)
      
      df %>% select(ap_noise, SNR, condition) %>% group_by(condition) %>% summarize(n = n(), 
                                                                                    "Noise" = round(mean(ap_noise), 3), 
                                                                                    "sd SNR" = round(sd(SNR), 2),
                                                                                    "SNR" = round(mean(SNR), 2)
                                                                                    ) %>% 
        relocate("SNR", .before = "sd SNR") %>%
        mutate("SNR > 3.75" = "SNR" > 3.75) %>%
        kable(., caption = "SNR Quality Check - Is SNR sufficiently high (> 3.75) to minimize false positive detection?") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        column_spec(1, bold = T) %>%
        add_footnote("See Salmanpour et al 2010") %>%
        print()
      
    }
    
    # ---- Output: Absolute Summary Data Table ----
    output$absolute_summary <- DT::renderDataTable({
      req(values$summary_df, input$conditions)
      
      summary_df <- values$summary_df

      #browser()
      
      temp <- summary_df %>% ungroup() %>%
        select(ID, condition, Absolute_Summary) %>%
        unnest(Absolute_Summary) %>%
        group_by(ID, condition) %>%
        mutate_all(unlist) %>%
        summarize_all(unique) %>%
        rename_all(~stringr::str_replace(., "Absolute_Summary_", "")) %>%
        pivot_longer(3:last_col()) %>%
        arrange(ID, name, condition) %>%
        mutate(name = factor(name))
        

      return(temp %>% 
               plotly::filter(condition %in% input$conditions)
      )
      
    }, options = list(pageLength = 10, autoWidth = TRUE), filter = "top", caption = "Summary Data") 
    
    
    # ---- Output: nCluster Summary Data Table ----
    output$nCluster_summary <- DT::renderDataTable({
      req(values$summary_df, input$conditions)
      
      summary_df <- values$summary_df
      
      #browser()
      
      temp <- summary_df %>% select(ID, condition, ncluster, nCluster_Summary) %>%
        unnest(nCluster_Summary) %>%
        group_by(ID, condition, ncluster) %>%
        mutate_all(unlist) %>%
        summarize_all(unique) %>%
        pivot_longer(4:last_col()) %>%
        arrange(ID, name, ncluster, condition) %>% 
        mutate(name = factor(name))
      
      return(temp %>% 
               plotly::filter(condition %in% input$conditions)
               
      )
      
    }, options = list(pageLength = 10, autoWidth = TRUE), filter = "top", caption = "Summary Data of Normalized Clusters")
    
    
    # ---- Re-cluster data when conditions change ----
    observe(priority = 1, {
      
      req(values$df, input$conditions, input$base)
    
      df <- values$df %>% 
        select(!c(mid, lower_bound))

      #browser()
      
      base_cond <- if(input$base %in% input$conditions) {
        input$base
      } else {
        input$conditions[[1]]
      }
      
      breaks <- find_cluster_bins2(x = df %>% 
                                    plotly::filter(condition %in% input$conditions) %>%
                                    plotly::filter(ap_include == TRUE & ap_keep == TRUE),
                                   base_cond = base_cond) # Helper function
      
      h <- df$ap_amplitude[df$ap_include == TRUE & df$ap_keep == TRUE & df$cluster %in% input$conditions] %>% 
        hist(breaks = breaks, plot = FALSE)  
      
      ap_base_max <- df %>% 
        plotly::filter(condition == base_cond) %>%
        plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>% select(ap_amplitude) %>% max(., na.rm = TRUE)
      
      df <- df %>% 
        ungroup() %>% 
        mutate(breaks = cut(.$ap_amplitude, breaks = h$breaks), 
               cluster = cut(.$ap_amplitude, breaks = h$breaks, labels = FALSE, ordered_result = TRUE),
               ap_namplitude = (ap_amplitude/ap_base_max)*100) %>%
        relocate(ap_namplitude, .after = ap_amplitude) %>%
        group_by(cluster) %>% 
        nest()
      
        df$mid <- h$mids[df$cluster]
        df$lower_bound <- h$breaks[df$cluster]
      
       df <- df %>% 
          ungroup() %>% 
          unnest(data) %>%
          relocate(ID, condition) %>%
          arrange(ap_no)
       
       df <- df %>%
         mutate(ncluster = (cluster/max(cluster, na.rm = TRUE) * 10) %>% ceiling()) %>%
         relocate(ncluster, .after = cluster)
       
        if("beat_no" %in% colnames(df)) {
          burstamp_base_max <- df %>% 
            plotly::filter(condition == base_cond) %>%
            plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>% select(burst_amplitude) %>% max(., na.rm = TRUE)
          burstarea_base_max <- df %>% 
            plotly::filter(condition == base_cond) %>%
            plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>% select(burst_area) %>% max(., na.rm = TRUE)
          
          
          df <- df %>% arrange(beat_no, ap_no) %>% 
            mutate(burst_namplitude = (burst_amplitude/burstamp_base_max)*100,
                   burst_narea = (burst_area/burstarea_base_max)*100) %>%
            relocate(mid, lower_bound, .after = breaks) %>%
            relocate(burst_namplitude, burst_narea, .after = burst_area)
        }

        myClusters <- c("ALL", unique(df$cluster) %>% sort(.))
        updateSelectizeInput(session, "Cluster", 
                             selected = isolate(input$Cluster), 
                             choices = myClusters)    
          
      values$df <- df
      
    })
    

    # ---- Output: Create AP Plot by amplitude cluster ----
    output$ap_characteristics <- renderPlotly({
      req(values$df)
      
      df <- values$df %>% 
        plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>%
        plotly::filter(condition %in% input$conditions)
      
      h1 <- df$inter_ap_int[df$inter_ap_int < 600] %>% 
        hist(breaks = "Scott", plot = FALSE)
      h1$probs <- h1$counts/sum(h1$counts)*100
      p1 <- plot_ly(x = h1$mids, y = h1$probs) %>% 
        add_bars(name = "Inter AP Interval") %>% 
        layout(xaxis = list(title = "Inter-AP Interval (ms)"), yaxis = list(title = "Probability (%)"))
      h2 <- df$ap_amplitude %>% 
        hist(breaks = "Scott", plot = FALSE)  
      h2$probs <- h2$counts/sum(h2$counts)*100  
      p2 <- plot_ly(x = h2$mids, y = h2$probs) %>% 
        add_bars(name = "AP Amplitude") %>% 
        layout(xaxis = list(title = "AP Amplitude (V)"), yaxis = list(title = "Probability (%)"))
      
      plot <- subplot(p1, p2, titleY = TRUE, titleX = TRUE, shareY = TRUE, margin = 0.07) %>% 
        layout(showlegend = FALSE)

      plot
    })

    
    # ---- Output: Create AP Plot ----
    output$ap_plot <- renderPlotly({
      req(values$df)
       # browser() 
      df <- values$df
      
      temp <- df %>% 
        unnest(data) %>%
        plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>%
        plotly::filter(condition %in% input$conditions)
      
      SelectedCluster <- input$Cluster
      
      #browser() 
    
      if (SelectedCluster == "ALL") {
        temp2 <- temp
      } else {
        temp2 <- temp %>% plotly::filter(cluster == as.numeric(SelectedCluster))
      }

    numAP <- length(unique(temp2$ap_no))  
        
    df_key <- highlight_key(temp2 %>% group_by(ap_no), ~ap_no)
      
      # Plot overlayed APs
      FigA <- plot_ly(source = "A") %>%
        add_lines(data = df_key,
                  x = ~sample,
                  y = ~ap_v,
                  color = ~factor(breaks),
                  colors = c("blue", "red"),

                  line = list(width = 2),
                  legendgroup = ~ap_keep,
                  text = ~paste("AP no: ",
                                round(ap_no, 2), '<br>time:',
                                round(ap_time, 2)
                  ),
                  customdata = ~ap_no,
        ) %>%
        layout(xaxis = list(title = "", showticklabels = FALSE))
        
      
      # Calculate Average AP
      # Mean AP by condition and cluster
      mean <- temp2 %>%
        group_by(ID, sample) %>%
        summarise(n = n(), mean = mean(ap_v), sd = sd(ap_v), se = sd/sqrt(n)) %>%
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
        layout(title = paste0(numAP, " Action Potentials Detected"),
               yaxis = list(title = "MSNA (mV)", fixedrange = FALSE)
               )

      #### Plot mean plots by cluster and plot below All APs and Mean AP.
      mean_cluster <- temp %>%
        group_by(ID, breaks, sample) %>%
        summarise(n = n(), mean = mean(ap_v), sd = sd(ap_v), se = sd/sqrt(n)) %>%
        mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)

      # FigC <- mean_cluster %>% ungroup() %>%
      #   group_by(breaks) %>%
      #   group_map(~ plot_ly(data=.,
      #                       x = ~sample,
      #                       y = ~mean,
      #                       showlegend = FALSE,
      #                       #color = ~breaks,
      #                       type = "scatter", mode="lines") %>%
      #               add_ribbons(  x = ~sample,
      #                             ymin = ~lower95,
      #                             ymax = ~upper95,
      #                             name = "95%CI",
      #                             showlegend = FALSE,
      #                             line = list(color = "grey", dash = 'dot', width = 1),
      #                             fillcolor = 'rgba(7, 164, 181, 0.2)'
      #               ) %>%
      #               layout(yaxis = list(title = "MSNA (mV)", fixedrange = FALSE),
      #                      xaxis = list(title = "", showticklabels = FALSE))
      #   ) %>%
      #   subplot(nrows = 1, shareX = TRUE, shareY=TRUE)
      
      #browser()
      
      panel <- . %>% plot_ly(data=.,
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
        add_annotations(
          text = ~unique(breaks),
          x = 0.5,
          y = 1,
          yref = "paper",
          xref = "paper",
          yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 10)
        ) %>%
        layout(
          showlegend = FALSE,
          shapes = list(
            type = "rect",
            x0 = 0,
            x1 = 1,
            xref = "paper",
            y0 = 0,
            y1 = 16,
            yanchor = 1,
            yref = "paper",
            ysizemode = "pixel",
            fillcolor = toRGB("gray80"),
            line = list(color = "transparent")),
          yaxis = list(title = "MSNA (mV)", 
                       fixedrange = FALSE),
          xaxis = list(title = "", 
                       showticklabels = FALSE))
        
      
      FigC <- mean_cluster %>% ungroup() %>%
           group_by(breaks) %>%
            do(p = panel(.)) %>%
             subplot(nrows = ceiling(NROW(.)/5), shareX = TRUE, shareY=TRUE)

      subplot(Fig, FigC, nrows = 2) %>%
        layout(clickmode = "event+select") %>%
        highlight(on = "plotly_hover", off = "plotly_relayout",  
                  selected = attrs_selected(showlegend = FALSE), debounce = 250) %>%
        event_register('plotly_click')

      
    })
  

    # ---- Output: Create AP Plot by amplitude cluster ----
    output$ap_characteristics2 <- renderPlotly({
      req(values$df, input$conditions, values$beat, values$signals25)
      
      df <- values$df %>% 
        plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>%
        plotly::filter(condition %in% input$conditions)
      
      df <- df %>% select(condition, ap_no, ap_amplitude, ap_latency, inter_ap_int, SNR) %>%
        pivot_longer(cols = c(ap_amplitude:SNR), names_to = "variable")
      
      p <- df %>% ggplot(aes(x = value, fill = condition)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~variable, scales = "free", labeller = as_labeller(c(ap_amplitude = "Amplitude (mV)", 
                                                                        ap_latency = "Latency (s)", 
                                                                        inter_ap_int = "Inter-AP Interval (ms)", SNR = "Signal-to-Noise Ratio"))) + 
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        ggtitle("Action Potential Characteristics") +
        theme_bw() +
        theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x = element_blank())
      ggplotly(p)
    })
    
    # ---- Output: Create time series plot ----
    output$TimeSeries <- renderPlotly({
      req(values$df, input$conditions2, values$beat, values$signals25)
      
      #browser()
      
      df <- values$df %>% 
        plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>%
        plotly::filter(condition %in% input$conditions2)
      
      t_min <- df$ap_time %>% min() 
      t_max <- df$ap_time %>% max()
      
      signals25 <- values$signals25 %>% plotly::filter(MSNA_time >= t_min & MSNA_time <= t_max)
      ap_linerange <- df$ap_time %>% cbind(ap_time = ., lower = 0, upper = 1) %>% as.data.frame() %>% pivot_longer(c("upper", "lower"))
      
      signals10 <- values$signals10khz %>% plotly::filter(time >= t_min & time <= t_max)
      
      figA <- plot_ly() %>%
        add_lines(data = signals25,
                  x = ~MSNA_time,
                  y = ~BP,
                  name = "BP",
                  line = list(color = "black")
        ) %>%
        layout(yaxis = list(title = "BP (mm Hg)"))
      
      figB <- plot_ly() %>%
        add_lines(data = signals25,
                  x = ~MSNA_time,
                  y = ~ECG,
                  name = "ECG",
                  line = list(color = "black")
        ) %>%
        layout(yaxis = list(title = "ECG (V)"))
      
      figC <- plot_ly() %>%
        add_lines(data = signals25,
                  x = ~MSNA_time,
                  y = ~MSNA_vjust,
                  name = "Neurogram",
                  line = list(color = "black")
        ) %>%
        layout(yaxis = list(title = "MSNA (V)"))
 
      figD <- plot_ly(ap_linerange) %>% group_by(ap_time) %>%
        add_lines(
                  x = ~ap_time,
                  y = ~value,
                  name = "AP",
                  line = list(color = "black")
        ) %>%
        layout(yaxis = list(title = "AP", fixedrange = FALSE),
               xaxis = list(title = "Time (s)")) %>%
        rangeslider()
      
      subplot(figA, figB, figC, figD, nrows = 4, shareX =TRUE, titleY = TRUE, heights = c(0.25,0.25,0.4,0.1))
      
    })
    
    
  # ---- Output: Correlograms ----
    output$correlograms <- renderPlot({
      
      req(values$df, input$conditions, values$beat, values$signals25)
      
     # browser()
      
      df <- values$df 
      
      ap_time <- df %>% plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>% plotly::filter(condition %in% input$conditions) %>% 
        group_by(condition, ap_no) %>% 
        summarise(TIME = unique(ap_time)) %>% select(!ap_no) %>% group_by(condition) %>% nest(spike = TIME)
      R_time <- df %>% plotly::filter(condition %in% input$conditions) %>% group_by(condition, beat_no) %>% summarise(TIME = unique(beat_time)) %>% select(!beat_no) %>% 
        group_by(condition) %>% nest(R = TIME)
      burst_time <- df %>% plotly::filter(condition %in% input$conditions) %>% select(condition, beat_no, burst_time) %>% plotly::filter(!is.na(burst_time)) %>% group_by(condition, beat_no) %>%
        summarise(TIME = unique(burst_time)) %>% select(!beat_no) %>% group_by(condition) %>% nest(burst = TIME)
      
    cor_df <- full_join(R_time, burst_time, by = "condition") %>% full_join(., ap_time, by = "condition")
      
    cor_df <- cor_df %>% mutate(burst_ECG = map2(.x = R, .y = burst, .f = correlation_histogram),
                                ap_ECG = map2(.x = R, .y = spike, .f = correlation_histogram))
      
    FigA <- cor_df %>% select(condition, burst_ECG) %>% unnest(burst_ECG) %>% 
      ggplot(aes(x = burst_ECG)) +
      geom_histogram(binwidth = 0.05) + #50ms binwidth
      geom_vline(xintercept = 0, linetype = "dotted") +
      scale_x_continuous(expand = c(0,0), limits = c(-3, +3), breaks = c(-3,-2,-1,0,1,2,3)) +
      labs(x = "Time (s)", y = "Counts/bin") +
      facet_grid(condition ~ .) +
      ggtitle("MSNA Burst vs ECG R") +
      theme_bw()
    
    FigB <- cor_df %>% select(condition, ap_ECG) %>% unnest(ap_ECG) %>% 
      ggplot(aes(x = ap_ECG)) +
      geom_histogram(binwidth = 0.05) + #50ms binwidth
      geom_vline(xintercept = 0, linetype = "dotted") +
      scale_x_continuous(expand = c(0,0), limits = c(-3, +3), breaks = c(-3,-2,-1,0,1,2,3)) +
      labs(x = "Time (s)", y = "Counts/bin") +
      facet_grid(condition ~ .) +
      ggtitle("MSNA AP vs ECG R") +
      theme_bw()
    
    ggarrange(FigA, FigB)

    })
    
    
    
    
    
    # ---- Observe: plot click ----
    observeEvent(event_data("plotly_click", source = "A"), { #input$toggle, {
      
      #browser()
      
      eventData <- event_data("plotly_click", source = "A")
      
      if ("customdata" %in% names(eventData)) {
        
        df <- values$df
        ap_no <- eventData$customdata
        
        df$ap_keep[which(df$ap_no %in% ap_no)] <- !df$ap_keep[which(df$ap_no %in% ap_no)]
        
        values$df <- df
        
      }
      
    })
    
    # ---- Observe: table click ----
    observeEvent(input$data_rows_selected, {

      df_selected <- input$data_rows_selected
      df <- values$df
      df$ap_keep[df_selected] <- !df$ap_keep[df_selected]
      values$df <- df
      
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
            
            df <- values$df %>% unnest(data) %>% plotly::filter(ap_include == TRUE & ap_keep == TRUE) %>%
              plotly::filter(condition %in% input$conditions) %>% select(ID, condition, ap_time, ap_no, sample, ap_v)
            df2 <- values$df %>% select(!data) %>% plotly::filter(condition %in% input$conditions)
            df2 <- df2[(df2$ap_include != FALSE | is.na(df2$ap_include)) & (df2$ap_keep != FALSE | is.na(df2$ap_keep)),]

            if("beat_no" %in% colnames(df)) {            
              fileName <- c(paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-all_aps.csv", sep = ""),
                            paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-all_aps_summary.csv", sep = ""),
                           paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-mean_ap.csv", sep = ""),
                           paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-mean_cluster_ap.csv", sep = ""),
                          paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-summary.csv", sep = ""),
                          paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-summary_ncluster.csv", sep = "")) 
              } else {
              fileName <- c(paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-all_aps.csv", sep = ""),
                            paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-all_aps_summary.csv", sep = ""),
                            paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-mean_ap.csv", sep = ""),
                            paste(tools::file_path_sans_ext(input$Signals_10KHz$name), "-mean_cluster_ap.csv", sep = ""))              
              }
            
            
            mean <- df %>%
              group_by(ID, condition, sample) %>%
              summarise(n = n(), mean = mean(ap_v), sd = sd(ap_v), se = sd/sqrt(n)) %>%
              mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)

            mean_cluster <- df %>%
              group_by(ID, condition, cluster, sample) %>%
              summarise(n = n(), mean = mean(ap_v), sd = sd(ap_v), se = sd/sqrt(n)) %>%
              mutate(upper95 = mean + 1.96*se, lower95 = mean - 1.96*se)

            if("beat_no" %in% colnames(df)) {
              summary_df <- values$summary_df %>% ungroup() %>%
                select(ID, condition, Absolute_Summary) %>%
                unnest(Absolute_Summary) %>%
                group_by(ID, condition) %>%
                mutate_all(unlist) %>%
                summarize_all(unique) %>%
                rename_all(~stringr::str_replace(., "Absolute_Summary_", "")) %>%
                pivot_longer(3:last_col()) %>%
                arrange(ID, name, condition) %>%
                mutate(name = factor(name))
              
              ncluster_df <- values$summary_df %>% select(ID, condition, ncluster, nCluster_Summary) %>%
                unnest(nCluster_Summary) %>%
                group_by(ID, condition, ncluster) %>%
                mutate_all(unlist) %>%
                summarize_all(unique) %>%
                pivot_longer(4:last_col()) %>%
                arrange(ID, name, ncluster, condition) %>% 
                mutate(name = factor(name))
            }

              write_csv(df, fileName[[1]])
              write_csv(df2, fileName[[2]])
              write_csv(mean, fileName[[3]])
              write_csv(mean_cluster, fileName[[4]])
              
              if("beat_no" %in% colnames(df)) {
              write_csv(summary_df, fileName[[5]])
              write_csv(ncluster_df, fileName[[6]])
              }

          #create the zip file
          zip(file,fileName)
        }
     )
     
 }

# Create Shiny app ----
thematic::thematic_shiny()
shinyApp(ui, server)
