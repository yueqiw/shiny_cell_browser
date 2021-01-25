library(shiny)
library(shinyjs)
library(plotly)
library(DT)
library(shinythemes)

json_file <- rjson::fromJSON(file = './data/config.json')
json_config <- json_file$config
window_title <- json_config$ui_title
title_link_text <- json_config$title_link_text
title_link_url <- json_config$title_link_url

ui <- fluidPage(
  tags$script(HTML(
   '$(document).ready(function () {
    $.getJSON("https://ipapi.co/json/", function (client) {
        Shiny.onInputChange("client", client);
    });
});'
  )),
  theme = shinytheme("yeti"),
  titlePanel(
    fluidRow(
            column(9, 
              div(
                h3(window_title, style="display:inline;"),
                tags$a(
                  h4(title_link_text, style="display:inline;"),
                  href=title_link_url, target="_blank"
                )
              )
            ),  
            column(3, div(
                        tags$a(h5("Browser App", style="display:inline;color:black;vertical-align:middle;"),
                        href="https://github.com/yueqiw/shiny_cell_browser", target="_blank"),
                        tags$a(img(height = 24, width = 24, src = "GitHub-Mark-64px.png", style="display:inline;vertical-align:middle;"),
                        href="https://github.com/yueqiw/shiny_cell_browser", target="_blank")
                    ),
                    align="right")
        ), windowTitle = window_title  # this is to display on the headers of web browser.
  ),
  shinyjs::useShinyjs(),
  sidebarLayout(fluid = TRUE,
              sidebarPanel(width = 2,
                            selectInput(inputId = "selected_dataset", label = "Dataset", choices = NULL),
                            selectizeInput(inputId = "selected_gene", label = "Gene", choices = NULL, options = list(placeholder = 'Select a gene')),
                            textAreaInput(inputId = "gene_list", label = "Gene List (up to 10 genes)", placeholder = "Up to 10,\n1 per line"),
                            fluidRow(column(12, actionButton(inputId = "gene_list_submit", label = "Dot plot genes"), align = "center")),
                            h5(strong("Table controls")),
                            fluidRow(column(12, actionButton(inputId = "reset_table", label = "Reset table"), align = "center")),
                            shinyjs::hidden(
                              textInput(inputId = "hidden_selected_gene", label = "Hidden Selected Gene", value = NULL)
                            ),
                            shinyjs::hidden(
                              textAreaInput(inputId = "hidden_gene_list", label = NULL, value = NULL)
                            ),
                            shinyjs::hidden(
                              textInput(inputId = "hidden_selected_cluster", label = "Hidden Selected Cluster", value = NULL)
                            )
              ),
              mainPanel(fluid = TRUE, width = 10,
                        fluidRow(width = 12,
                                  column(width = 4,
                                        plotlyOutput(outputId = "cluster_plot", width = "100%", height = "auto")
                                  ),
                                  column(width = 4,
                                        plotlyOutput(outputId = "expression_plot", width = "100%", height = "auto")
                                  ),
                                  column(width = 4,
                                        plotlyOutput(outputId = "dot_plot", width = "100%", height = "auto")
                                  )
                        ),
                        fluidRow(width = 12,
                                  textOutput('cluster_gene_table_title'),
                                  br(),
                                  DTOutput('cluster_gene_table')
                        )
              )
  ),
  windowTitle = window_title
)