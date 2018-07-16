library(shiny)
library(shinyjs)
library(plotly)

ui_plot_width <- "100%"
ui_plot_height <- "auto"

window_title <- 'Shiny Single Cell Browser'

shinyUI(fluidPage(
    useShinyjs(),
    titlePanel(
        fluidRow(
            column(6, window_title)  # we can potentially add more info here
        ), windowTitle = window_title  # this is to display on the headers of web browser.
    ),
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

            fluidRow(class = "loading",
                column(12,
                    conditionalPanel(
                                condition="($('html').hasClass('shiny-busy'))",
                                img(height=16, src="horizontal_spinner.gif")
                    ), align="center", style="margin-top:8px;"
                ), style="height:32px;"
            ),

            radioButtons("plot_type", "Visualization Mode:", choices = c("Basic" = "ggplot2", "Interactive" = "plotly"),
                            selected = "ggplot2", inline = TRUE),
            radioButtons("layout_type", "Layout:", choices = c("Vertical" = "vertical", "Horizontal" = "horizontal"),
                            selected = "vertical", inline = TRUE),
            #checkboxInput("featureplot_check", "Always show 2D gene plot", value = FALSE),
            checkboxInput("auto_scaling_check", "Auto scale figure size", value = TRUE),
            # checkboxGroupInput("checkboxes", label = NULL,
            #     choices = list("Always show 2D gene plot" = "featureplot_check",
            #                     "Auto scale figure size" = "auto_scaling_check"),
            #     selected = c("auto_scaling_check")),
            sliderInput("plot_width", "Manual scale figure size",
                  min = 200, max = 800, value = 500, step = 10, ticks = FALSE),
            div(),  # or hr()
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

            # style = "position:fixed;width:inherit;",  # floating sidebar prevents reaching to the bottom when it's longer than window height.
            # tags$head(tags$style(".loading{height:30px;}"  # two ways to add styles: 1. here, 2. use style="xxx" in each box.
            #   )),
            width=2
        ),
        mainPanel(
            # fluidRow(
            #     column(6, uiOutput("title1"), align = 'center'),
            #     column(6, uiOutput("title2"), align = 'center')
            # ),
            # fluidRow(
            #     column(6, uiOutput("description1"), align = 'center'),
            #     column(6, uiOutput("description2"), align = 'center')
            # ),
            uiOutput("main_panel")
            # fluidRow(
            #     column(4, uiOutput("plot_ui_r1c1"), align="center"),
            #     column(4, uiOutput("plot_ui_r1c2"), align="center"),
            #     column(4, uiOutput("plot_ui_r1c3"), align="center")
            # ),
            # fluidRow(
            #     column(4, uiOutput("plot_ui_r2c1"), align="center"),
            #     column(4, uiOutput("plot_ui_r2c2"), align="center"),
            #     column(4, uiOutput("plot_ui_r2c3"), align="center")
            # ),
            # fluidRow(
            #     column(4, uiOutput("plot_ui_r3c1"), align="center"),
            #     column(4, uiOutput("plot_ui_r3c2"), align="center"),
            #     column(4, uiOutput("plot_ui_r3c3"), align="center")
            # )
            , width = 10
        ), position = c('left')
    )
)
)
