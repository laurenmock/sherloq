
library(sherloq)
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(DT)
library(shinydashboard)

ui_LoD_precision_profile <-
  tabItem(tabName = "LoD_precision_profile",

          titlePanel("LoD (Precision Profile) Calculator"),

          sidebarLayout(

            sidebarPanel(

              fileInput("file_LoDPP_IN",
                        "Choose CSV File",
                        multiple = FALSE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),

              selectInput("col_lot_LoDPP_IN", "Column with reagent lot number:", ""),
              selectInput("col_sample_LoDPP_IN", "Column with unique values for each sample:", ""),
              conditionalPanel(
                condition = "input.col_sample_LoDPP_IN == 'N/A'",
                numericInput("n_samples_LoDPP_IN",
                             label = "Number of samples:",
                             value = "")),
              selectInput("col_avg_LoDPP_IN", "Column with mean measurement values:", ""),
              selectInput("col_sd_wl_LoDPP_IN", "Column with within-lab precision:", ""),
              numericInput("LoB_LoDPP_IN",
                           label = "LoB:",
                           value = ""),
              numericInput("n_meas_LoDPP_IN",
                           label = "Total number of measurement values (per reagent lot, if separated):",
                           value = ""),
              radioButtons("seplots_LoDPP_IN",
                           label = "Separate or pool reagent lots?",
                           choices = list("Follow CLSI guidelines (pool if 4+ lots)" = "CLSI",
                                          "Separate lots" = "separate",
                                          "Pool lots" = "pool"),
                           selected = "CLSI"),
              numericInput("beta_LoDPP_IN",
                           label = "Beta:",
                           value = "0.05"),
              radioButtons("model_LoDPP_IN",
                           label = "Precision profile model:",
                           choices = list("Lowest AIC" = "lowestAIC",
                                          "Linear" = "linear",
                                          "Quadratic" = "quadratic",
                                          "Sadler" = "sadler"),
                           selected = "lowestAIC"),
              conditionalPanel(
                condition = "input.model_LoDPP_IN == 'sadler'",
                withMathJax(),
                h4(helpText("Sadler model structure: \\(y = (\U03B2_0 + \U03B2_1x)^{\U03B2_2}\\)")),
                h5("The Sadler model fit relies on the `nls` function and requires inital estimates
                   for each coefficient."),
                numericInput("sadler_start1_LoDPP_IN",
                             label = helpText("Initial estimate for \\(\U03B2_0\\):"),
                             value = ""),
                numericInput("sadler_start2_LoDPP_IN",
                             label = helpText("Initial estimate for \\(\U03B2_1\\):"),
                             value = ""),
                numericInput("sadler_start3_LoDPP_IN",
                             label = helpText("Initial estimate for \\(\U03B2_2\\):"),
                             value = "")),
              actionButton("go_LoDPP", "Calculate")

            ), # end sidebar panel

            mainPanel(

              tableOutput("results_LoDPP_OUT"),
              plotOutput("plot_LoDPP_OUT"),
              tableOutput("mod_coeffs_LoDPP_OUT"),
              dataTableOutput("pred_tab_LoDPP_OUT")

            ) # end main panel

          ) # end sidebar layout

  ) # end of LoD_precision_profile tab
