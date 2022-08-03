
library(sherloq)
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(DT)
library(shinydashboard)

ui_LoD_probit <-
  tabItem(tabName = "LoD_probit",

          titlePanel("LoD (Probit) Calculator"),

          sidebarLayout(

            sidebarPanel(

              fileInput("file_LoDP_IN",
                        "Choose CSV File",
                        multiple = FALSE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),

              selectInput("col_lot_LoDP_IN", "Column with reagent lot number:", ""),
              selectInput("col_conc_LoDP_IN", "Column with concentration:", ""),
              checkboxInput("log10_trans_LoDP_IN",
                            "Perform a log10 transformation of concentration values before fitting probit models"),
              radioButtons("format_LoDP_IN",
                           label = "This data...",
                           choices = list("has a row for each measurement (long format)" = "long",
                                          "is already grouped by concentration (aggregated)" = "wide"),
                           selected = "long"),
              conditionalPanel(
                condition = "input.format_LoDP_IN == 'long'",
                selectInput("col_01_LoDP_IN",
                            "Column with sample call (0 if negative, 1 if positive):", "")),
              conditionalPanel(
                condition = "input.format_LoDP_IN == 'wide'",
                selectInput("col_obs_pos_LoDP_IN",
                            "Column with the number of positive calls per concentration level:", ""),
                selectInput("col_tot_LoDP_IN",
                            "Column with the total number of calls per concentration level:", "")),

              numericInput("LoB_LoDP_IN",
                           label = "LoB:",
                           value = ""),
              radioButtons("seplots_LoDP_IN",
                           label = "Separate or pool reagent lots?",
                           choices = list("Follow CLSI guidelines (pool if 4+ lots)" = "CLSI",
                                          "Separate lots" = "separate",
                                          "Pool lots" = "pool"),
                           selected = "CLSI"),
              numericInput("beta_LoDP_IN",
                           label = "Beta:",
                           value = "0.05"),
              actionButton("go_LoDP", "Calculate")

            ), # end sidebar panel

            mainPanel(

              tableOutput("results_LoDP_OUT"),
              plotOutput("plot_LoDP_OUT"),
              tableOutput("mod_coeffs_LoDP_OUT"),
              dataTableOutput("pred_tab_LoDP_OUT")

            ) # end main panel

          ) # end sidebar layout

  ) # end of LoD_probit tab
