library(shiny)
library(shinyjs)
library(plotly)
column_size <- 6
ui_plot_width <- "100%"
ui_plot_height <- "auto"

shinyUI(fluidPage(
    useShinyjs(),
    titlePanel('Shiny Single Cell Browser'),
    sidebarLayout(
        sidebarPanel(
            selectizeInput('dataset_1', h4('Datasets'),
                        selected = NULL,
                        choices = NULL,
                        options = list(placeholder = 'Select a dataset')
                    ),
            selectizeInput('dataset_2', NULL,
                        selected = NULL,
                        choices = NULL,
                        options = list(placeholder = 'Select a dataset')
                    ),
            selectizeInput('gene_symbol', h4('Gene (Ensembl)'),
                        selected = "",
                        choices = "",
                        options = list(placeholder = 'Select a gene')
                    ),
            textAreaInput("gene_list", NULL, "", placeholder = "Add gene list here",
                        height = '80px'),
            fluidRow(column(12, actionButton("gene_list_submit", "Plot Gene list"), align="center")),
            radioButtons("plot_type", "Plot type:", choices = c("Static" = "ggplot2", "Interactive" = "plotly"),
                            selected = "ggplot2", inline = TRUE),
            checkboxInput("featureplot_check", "Always show 2D gene plot", value = FALSE),
            checkboxInput("figure_scaling_check", "Auto scale figure size", value = TRUE),
            # checkboxGroupInput("checkboxes", label = NULL,
            #     choices = list("Always show 2D gene plot" = "featureplot_check",
            #                     "Auto scale figure size" = "figure_scaling_check"),
            #     selected = c("figure_scaling_check")),
            sliderInput("plot_width", "Manual scale figure size",
                  min = 200, max = 800, value = 500, step = 10, ticks = FALSE),
            hr(),
            radioButtons("figure_type", h4(HTML("Download Figure:")),
                  choiceNames = c("2D Cluster Plot", "Dot Plot", "2D Gene Plot"),
                  choiceValues = c("clusterplot", "dotplot", "featureplot"), selected = "dotplot", inline = FALSE),
            fluidRow(
                column(5, radioButtons("dataset", "Dataset:",
                      choices = c(1, 2), selected = 1, inline = FALSE)),
                column(7, radioButtons("filetype", "File Type:",
                       choices = c("pdf", "png"), selected = "pdf", inline = FALSE))
            ),
            fluidRow(column(12, downloadButton('download', 'Download'), align="center")),
            verbatimTextOutput("value"),  # for debugging purpose only.

            style = "position:fixed;width:inherit;",
            width=2
        ),
        mainPanel(
            fluidRow(
                column(column_size, uiOutput("title1"), align = 'center'),
                column(column_size, uiOutput("title2"), align = 'center')
            ),
            fluidRow(
                column(column_size, uiOutput("description1"), align = 'center'),
                column(column_size, uiOutput("description2"), align = 'center')
            ),
            fluidRow(
                column(column_size, uiOutput("clusterplot1_ui"), align="center"),
                column(column_size, uiOutput("clusterplot2_ui"), align="center")
            ),
            fluidRow(
                column(column_size, uiOutput("featureplot1_ui"), align="center"),
                column(column_size, uiOutput("featureplot2_ui"), align="center")
            ),
            fluidRow(
                column(column_size, uiOutput("dotplot1_ui"), align="center"),
                column(column_size, uiOutput("dotplot2_ui"), align="center")
            ), width = 9
        ), position = c('left')
    )
)
)
