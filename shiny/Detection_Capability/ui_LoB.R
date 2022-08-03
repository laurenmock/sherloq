
library(sherloq)
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(DT)
library(shinydashboard)

ui_LoB <-
  tabItem(tabName = "LoB",

          titlePanel("LoB Calculator"),

          sidebarLayout(

            sidebarPanel(

              fileInput("file_LoB_IN",
                        "Choose CSV File",
                        multiple = FALSE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),

              selectInput("col_lot_LoB_IN", "Column with reagent lot number:", ""),
              selectInput("col_sample_LoB_IN", "Column with unique values for each sample:", ""),
              selectInput("col_value_LoB_IN", "Column with measurement values:", ""),
              radioButtons("seplots_LoB_IN",
                           label = "Separate or pool reagent lots?",
                           choices = list("Follow CLSI guidelines (pool if 4+ lots)" = "CLSI",
                                          "Separate lots" = "separate",
                                          "Pool lots" = "pool"),
                           selected = "CLSI"),
              radioButtons("parametric_LoB_IN",
                           label = "Approach",
                           choices = list("Non-parametric (recommended)" = "nonpar",
                                          "Parametric" = "par"),
                           selected = "nonpar"),
              numericInput("alpha_LoB_IN",
                           label = "Alpha:",
                           value = "0.05"),
              radioButtons("plot_LoB_IN",
                           label = "Display measurement results with...",
                           choices = list("Box plot(s)" = "boxplot",
                                          "Histogram(s)" = "histogram"),
                           selected = "boxplot"),
              actionButton("go_LoB", "Calculate")

            ), # end sidebar panel

            mainPanel(

              tableOutput("results_LoB_OUT"),
              hr(),
              plotOutput("plot_LoB_OUT"),

            ) # end main panel

          ) # end sidebar layout

  ) # end of LoB tab
