library(shiny)
library(shinyjs)
library(plotly)
library(DT)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("U Brain Browser"),
  sidebarLayout(fluid=TRUE,
                sidebarPanel(width=2,
                             selectInput(inputId = "selected_dataset", label = "Dataset", choices = NULL),
                             selectizeInput(inputId = "selected_gene", label = "Gene", choices = NULL, options = list(placeholder = 'Select a gene')),
                             textAreaInput(inputId = "gene_list", label = "Gene List", placeholder = "Up to 10,\n1 per line"),
                             fluidRow(column(12, actionButton(inputId = "gene_list_submit", label = "Dot plot genes"), align="center")),
                             h5(strong("Table controls")),
                             fluidRow(column(12,actionButton(inputId = "reset_table", label = "Reset table"), align="center")),
                             shinyjs::hidden(
                               textInput(inputId = "hidden_selected_gene", label="Hidden Selected Gene", value=NULL)
                             ),
                             shinyjs::hidden(
                               textAreaInput(inputId = "hidden_gene_list", label = NULL, value=NULL)
                             ),
                             shinyjs::hidden(
                               textInput(inputId = "hidden_selected_cluster", label="Hidden Selected Cluster", value=NULL)
                             )
                ),
                mainPanel(fluid=TRUE,width=10,
                          fluidRow(width=12,
                                   column(width=4,
                                          plotlyOutput(outputId = "cluster_plot",width="100%")
                                   ),
                                   column(width=4,
                                          plotlyOutput(outputId = "expression_plot",width="100%")
                                   ),
                                   column(width=4,
                                          plotlyOutput(outputId = "dot_plot",width="100%")
                                   )
                          ),
                          fluidRow(width=12,
                                   textOutput('cluster_gene_table_title'),
                                   br(),
                                   DTOutput('cluster_gene_table')
                          )
                )
  )
)