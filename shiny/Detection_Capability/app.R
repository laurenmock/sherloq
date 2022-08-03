
#----- Shiny App for Detection Capability -----#

# this app is based on the sherloq package

# source UI functions
source("ui_LoB.R")
source("ui_LoD_classical.R")
source("ui_LoD_precision_profile.R")
source("ui_LoD_probit.R")
source("ui_LoQ_functional.R")
source("ui_LoQ_total_error.R")

# source server functions
# source("server_LoB.R")
# source("server_LoD_classical.R")

library(sherloq)
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(DT)
library(shinydashboard)


ui <- dashboardPage(

  dashboardHeader(title = "Detection Capability"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("LoB",
               tabName = "LoB",
               icon = icon("1")),
      menuItem("LoD (Classical)",
               tabName = "LoD_classical",
               icon = icon("2")),
      menuItem("LoD (Precision Profile)",
               tabName = "LoD_precision_profile",
               icon = icon("3")),
      menuItem("LoD (Probit)",
               tabName = "LoD_probit",
               icon = icon("4")),
      menuItem("LoQ (Functional)",
               tabName = "LoQ_functional",
               icon = icon("5")),
      menuItem("LoQ (Total Error)",
               tabName = "LoQ_total_error",
               icon = icon("6"))

    ) # end sidebar menu
  ), # end dashboard sidebar


  dashboardBody(
    tabItems(

      ui_LoB,
      ui_LoD_classical,
      ui_LoD_precision_profile,
      ui_LoD_probit,

      tabItem(tabName = "LoQ_functional",
              titlePanel("LoQ (Functional) Calculator")
      ), # end of LoQ_functional tab

      tabItem(tabName = "LoQ_total_error",
              titlePanel("LoQ (Total Error) Calculator")
      ) # end of LoQ_total_error tab

    )
  ) # end dashboard body
) # end dashboard page/UI



server <- function(input, output) {

  #----- LoB -----#

  observeEvent(input$file_LoB_IN,{

    # blank reactive to store loaded data
    reactiveDF_LoB <- reactiveValues(df = NULL)

    req(input$file_LoB_IN)

    # store loaded data in reactive
    reactiveDF_LoB$df <- read.csv(file = input$file_LoB_IN$datapath)

    # Update select input
    updateSelectInput(inputId = "col_lot_LoB_IN",
                      label = "Column with reagent lot number:",
                      choices  = c(colnames(reactiveDF_LoB$df), "N/A"))
    updateSelectInput(inputId = "col_sample_LoB_IN",
                      label = "Column with unique values for each sample:",
                      choices  = colnames(reactiveDF_LoB$df))
    updateSelectInput(inputId = "col_value_LoB_IN",
                      label = "Column with measurement values:",
                      choices  = colnames(reactiveDF_LoB$df))

    results <- eventReactive(input$go_LoB, {

      LoB(df = reactiveDF_LoB$df,
          col_lot = if(input$col_lot_LoB_IN == "N/A" | input$seplots_LoB_IN == "pool")
            {NULL}else{input$col_lot_LoB_IN},
          col_sample = input$col_sample_LoB_IN,
          col_value = input$col_value_LoB_IN,
          parametric = ifelse(input$parametric_LoB_IN == "nonpar",
                              FALSE, TRUE),
          alpha = input$alpha_LoB_IN,
          plot = input$plot_LoB_IN
      )
    })

    output$results_LoB_OUT <- renderTable({
      results()$LoB_values
    })

    output$plot_LoB_OUT <- renderPlot(
      results()$LoB_plot
    )

  }) # end LoB server



  #----- LoD (classical) -----#

  observeEvent(input$file_LoDC_IN,{

    # blank reactive to store loaded data
    reactiveDF_LoDC <- reactiveValues(df = NULL)

    req(input$file_LoDC_IN)

    # store loaded data in reactive
    reactiveDF_LoDC$df <- read.csv(file = input$file_LoDC_IN$datapath)

    # Update select input
    updateSelectInput(inputId = "col_lot_LoDC_IN",
                      label = "Column with reagent lot number:",
                      choices  = c(colnames(reactiveDF_LoDC$df), "N/A"))
    updateSelectInput(inputId = "col_sample_LoDC_IN",
                      label = "Column with unique values for each sample:",
                      choices  = colnames(reactiveDF_LoDC$df))
    updateSelectInput(inputId = "col_value_LoDC_IN",
                      label = "Column with measurement values:",
                      choices  = colnames(reactiveDF_LoDC$df))

    results <- eventReactive(input$go_LoDC, {

      LoD_classical(df = reactiveDF_LoDC$df,
                    col_lot = if(input$col_lot_LoDC_IN == "N/A" | input$seplots_LoDC_IN == "pool")
                      {NULL}else{input$col_lot_LoDC_IN},
                    col_sample = input$col_sample_LoDC_IN,
                    col_value = input$col_value_LoDC_IN,
                    LoB = input$LoB_LoDC_IN,
                    beta = input$beta_LoDC_IN,
                    plot = input$plot_LoDC_IN
      )
    })

    output$results_LoDC_OUT <- renderTable({
      results()$LoD_values
    })

    output$plot_LoDC_OUT <- renderPlot(
      results()$LoD_plot
    )
  }) # end LoD_classical server


  #----- LoD (precision profile) -----#

  observeEvent(input$file_LoDPP_IN,{

    # blank reactive to store loaded data
    reactiveDF_LoDPP <- reactiveValues(df = NULL)

    req(input$file_LoDPP_IN)

    # store loaded data in reactive
    reactiveDF_LoDPP$df <- read.csv(file = input$file_LoDPP_IN$datapath)

    # update input with column names
    updateSelectInput(inputId = "col_lot_LoDPP_IN",
                      label = "Column with reagent lot number:",
                      choices  = c(colnames(reactiveDF_LoDPP$df), "N/A"))
    updateSelectInput(inputId = "col_sample_LoDPP_IN",
                      label = "Column with unique values for each sample:",
                      choices  = c(colnames(reactiveDF_LoDPP$df), "N/A"))
    updateSelectInput(inputId = "col_avg_LoDPP_IN",
                      label = "Column with mean measurement values:",
                      choices  = colnames(reactiveDF_LoDPP$df))
    updateSelectInput(inputId = "col_sd_wl_LoDPP_IN",
                      label = "Column with within-lab precision:",
                      choices  = colnames(reactiveDF_LoDPP$df))

    results <- eventReactive(input$go_LoDPP, {

      LoD_precision_profile(df = reactiveDF_LoDPP$df,
                            col_lot = if(input$col_lot_LoDPP_IN == "N/A" | input$seplots_LoDPP_IN == "pool")
                              {NULL}else{input$col_lot_LoDPP_IN},
                            col_sample = if(input$col_sample_LoDPP_IN == "N/A")
                              {NULL}else{input$col_sample_LoDPP_IN},
                            col_avg = input$col_avg_LoDPP_IN,
                            col_sd_wl = input$col_sd_wl_LoDPP_IN,
                            n_measures = input$n_meas_LoDPP_IN,
                            n_samples = if(input$col_sample_LoDPP_IN == "N/A")
                              {input$n_samples_LoDPP_IN}else{NULL},
                            LoB = input$LoB_LoDPP_IN,
                            beta = input$beta_LoDPP_IN,
                            model = input$model_LoDPP_IN,
                            sadler_start = c(input$sadler_start1_LoDPP_IN,
                                             input$sadler_start2_LoDPP_IN,
                                             input$sadler_start3_LoDPP_IN)
      )
    })

    output$results_LoDPP_OUT <- renderTable({
      results()$LoD_values
    })

    output$plot_LoDPP_OUT <- renderPlot(
      results()$model_plot
    )

    output$mod_coeffs_LoDPP_OUT <- renderTable({
      lapply(results()$model, function(x) summary(x)$coef)
    })

    output$pred_tab_LoDPP_OUT <- renderDataTable({
      results()$model_predictions |>
        setNames(c("Reagent", "Mean", "SD_Pred"))
    })


  }) # end LoD_precision_profile server


  #----- LoD (probit) -----#

  observeEvent(input$file_LoDP_IN,{

    # blank reactive to store loaded data
    reactiveDF_LoDP <- reactiveValues(df = NULL)

    req(input$file_LoDP_IN)

    # store loaded data in reactive
    reactiveDF_LoDP$df <- read.csv(file = input$file_LoDP_IN$datapath)

    # update input with column names
    updateSelectInput(inputId = "col_lot_LoDP_IN",
                      label = "Column with reagent lot number:",
                      choices  = c(colnames(reactiveDF_LoDP$df), "N/A"))
    updateSelectInput(inputId = "col_conc_LoDP_IN",
                      label = "Column with concentration:",
                      choices  = colnames(reactiveDF_LoDP$df))
    updateSelectInput(inputId = "col_01_LoDP_IN",
                      label = "Column with sample call (0 if negative, 1 if positive):",
                      choices  = colnames(reactiveDF_LoDP$df))
    updateSelectInput(inputId = "col_obs_pos_LoDP_IN",
                      label = "Column with the number of positive calls per concentration level:",
                      choices  = colnames(reactiveDF_LoDP$df))
    updateSelectInput(inputId = "col_tot_LoDP_IN",
                      label = "Column with the total number of calls per concentration level:",
                      choices  = colnames(reactiveDF_LoDP$df))

    results <- eventReactive(input$go_LoDP, {

      LoD_probit(df = reactiveDF_LoDP$df,
                 col_lot = if(input$col_lot_LoDP_IN == "N/A" | input$seplots_LoDP_IN == "pool")
                   {NULL}else{input$col_lot_LoDP_IN},
                 col_conc = input$col_conc_LoDP_IN,
                 col_01 = if(input$format_LoDP_IN == "long"){input$col_01_LoDP_IN}else{NULL},
                 col_obs_pos = if(input$format_LoDP_IN == "wide"){input$col_obs_pos_LoDP_IN}else{NULL},
                 col_tot = if(input$format_LoDP_IN == "wide"){input$col_tot_LoDP_IN}else{NULL},
                 LoB = input$LoB_LoDP_IN,
                 beta = input$beta_LoDP_IN,
                 log10_trans = ifelse(input$log10_trans_LoDP_IN, TRUE, FALSE)
      )
    })

    output$results_LoDP_OUT <- renderTable({
      results()$LoD_values
    })

    output$plot_LoDP_OUT <- renderPlot(
      results()$hit_rate_plot
    )

    output$mod_coeffs_LoDP_OUT <- renderTable({
      lapply(results()$probit_model, function(x) summary(x)$coef)
    })

    output$pred_tab_LoDP_OUT <- renderDataTable({
      results()$model_predictions
    })

  }) # end LoD_probit server


} # end server

# Run the application
shinyApp(ui = ui, server = server)
