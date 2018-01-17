column_size <- 6
plot_size <- 500

pageWithSidebar(
    headerPanel('Gene Plot for organoid scSeq'),
    sidebarPanel(
        textInput('gene_symbol', h3(HTML('Gene Name<br/>(Ensembl, upper-case)')), value = 'TBR1'),
        actionButton("submit", "Submit"),
        width=2
    ),
    mainPanel(
        verbatimTextOutput("value_1"),
        tableOutput("exp_level_1"),
        fluidRow(
            column(column_size, plotOutput("plot1", height = plot_size)),
            column(column_size, plotOutput("plot2", height = plot_size))
        ),
        verbatimTextOutput("value_2"),
        tableOutput("exp_level_2"),
        fluidRow(
            column(column_size, plotOutput("plot3", height = plot_size)),
            column(column_size, plotOutput("plot4", height = plot_size))
        )
    )
)
