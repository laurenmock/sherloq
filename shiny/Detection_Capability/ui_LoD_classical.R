
library(sherloq)
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(DT)
library(shinydashboard)

ui_LoD_classical <-
  tabItem(tabName = "LoD_classical",

          titlePanel("LoD (Classical) Calculator"),

          sidebarLayout(

            sidebarPanel(

              fileInput("file_LoDC_IN",
                        "Choose CSV File",
                        multiple = FALSE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),

              selectInput("col_lot_LoDC_IN", "Column with reagent lot number:", ""),
              selectInput("col_sample_LoDC_IN", "Column with unique values for each sample:", ""),
              selectInput("col_value_LoDC_IN", "Column with measurement values:", ""),
              numericInput("LoB_LoDC_IN",
                           label = "LoB:",
                           value = ""),
              radioButtons("seplots_LoDC_IN",
                           label = "Separate or pool reagent lots?",
                           choices = list("Follow CLSI guidelines (pool if 4+ lots)" = "CLSI",
                                          "Separate lots" = "separate",
                                          "Pool lots" = "pool"),
                           selected = "CLSI"),
              numericInput("beta_LoDC_IN",
                           label = "Beta:",
                           value = "0.05"),
              radioButtons("plot_LoDC_IN",
                           label = "Display measurement results with...",
                           choices = list("Box plot(s)" = "boxplot",
                                          "Histogram(s)" = "histogram"),
                           selected = "boxplot"),
              actionButton("go_LoDC", "Calculate")

            ), # end sidebar panel

            mainPanel(

              tableOutput("results_LoDC_OUT"),
              plotOutput("plot_LoDC_OUT"),

            ) # end main panel

          ) # end sidebar layout

  ) # end of LoD_classical tab
