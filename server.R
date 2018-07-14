library(tidyverse)
library(Seurat)
library(rjson)
source("./utils.R")

calc_pt_size <- function(n) {30 / n^0.5}
plot_inch <- 4
dpi <- 125  # for web display. The saved plots have higher dpi
png.arguments <- c(4,4,300)
rasterize_ncells <- 1500

json_data <- fromJSON(file = './config.json')
json_data <- json_data$data
datasets <- 1:length(json_data)
dataset_names <- sapply(json_data, function(x) x$name)
dataset_selector <- as.list(c(datasets, "no_data"))
names(dataset_selector) <- c(dataset_names, "None")

# fixed_plot_size <- 500

function(input, output, session) {
    #dataset_selector <- c(1:4)
    print(dataset_selector)

    updateSelectizeInput(session, 'dataset_1', choices = dataset_selector, selected = dataset_selector[[1]], server = F)
    updateSelectizeInput(session, 'dataset_2', choices = dataset_selector, selected = dataset_selector[[2]], server = F)


    data_1 <- reactiveValues()
    data_2 <- reactiveValues()
    seurat_data_1 <- reactiveValues()
    seurat_data_2 <- reactiveValues()
    ncells_1 <- reactiveValues()
    ncells_2 <- reactiveValues()
    pt_size_1 <- reactiveValues()
    pt_size_2 <- reactiveValues()
    colors_1 <- reactiveValues()
    colors_2 <- reactiveValues()
    all_genes <- reactiveValues()

    dataset_1 <- reactive(input$dataset_1)
    dataset_2 <- reactive(input$dataset_2)

    observeEvent({
        input$dataset_1
        input$dataset_2
    }, {
        #print(dataset_1())
        if (dataset_1() %in% datasets) {
            data1 <<- json_data[[as.numeric(dataset_1())]]
            #print(data1)
            seurat_data_1 <<- readRDS(data1$file)
            seurat_data_1 <<- SetAllIdent(seurat_data_1,  data1$clusters)
            ncells_1 <<- length(seurat_data_1@cell.names)
            pt_size_1 <<- calc_pt_size(ncells_1)
            data1$embedding <<- data1$embedding
            colors_1 <<- seurat_data_1@misc[[sprintf("%s_colors",  data1$clusters)]]
            if (is.null(colors_1)) colors_1 <<- rainbow(n_distinct(seurat_data_1@ident))
            genes_1 <- rownames(seurat_data_1@data)
        } else {
            seurat_data_1 <<- NULL
            genes_1 <- NULL
        }

        if (dataset_2() %in% datasets) {
            data2 <<- json_data[[as.numeric(dataset_2())]]
            seurat_data_2 <<- readRDS(data2$file)
            seurat_data_2 <<- SetAllIdent(seurat_data_2, data2$clusters)
            ncells_2 <<- length(seurat_data_2@cell.names)
            pt_size_2 <<- calc_pt_size(ncells_2)
            data2$embedding <<- data2$embedding
            colors_2 <<- seurat_data_2@misc[[sprintf("%s_colors", data2$clusters)]]
            if (is.null(colors_2)) colors_2 <<- rainbow(n_distinct(seurat_data_2@ident))
            genes_2 <- rownames(seurat_data_2@data)
        } else {
            seurat_data_2 <<- NULL
            genes_2 <- NULL
        }

        all_genes <<- as.list(c("", union(genes_1, genes_2)))
        #print(all_genes[1:5])

        updateSelectizeInput(session, 'gene_symbol', choices = all_genes, server = F)
    })

    gene_names <- reactiveValues()
    observeEvent(input$gene_symbol, {
        gene_names <<- input$gene_symbol
    })
    observeEvent(input$gene_list_submit, {
        gene_names <<- trimws(strsplit(input$gene_list, '\n')[[1]])
    })

    #output$value <- renderText({ gene_names })
    #checked_boxes <- reactive(input$checkboxes)  # is none selected %in% doesn't work
    always_show_featureplot <- reactive(input$featureplot_check)
    auto_scaling <- reactive(input$figure_scaling_check)
    fixed_plot_size <- reactive(input$plot_width)
    plot_size <- reactiveValues()
    plot_width <- reactiveValues()

    observeEvent({
        input$figure_scaling_check
        input$plot_width
        session$clientData$output_clusterplot1_width
    }, {
        if(auto_scaling()) {
            plot_size <<- "auto"
            plot_width <<- function() {session$clientData$output_clusterplot1_width}
        } else {
            plot_size <<- function() {fixed_plot_size()}
            plot_width <<- function() {fixed_plot_size()}
        }

    }, ignoreNULL=FALSE)  # if using checkboxGroupInput, need this option when no choice is selected

    plot_1_ptsize <- function() {(plot_width() / dpi / plot_inch)^2 * pt_size_1}
    plot_2_ptsize <- function() {(plot_width() / dpi / plot_inch)^2 * pt_size_2}

    scale_factor <- function() {plot_width() / dpi/ 6.33}

    dotplot <- function(data, genes, scale_factor) {
            dotplot_width <- function() {scale_factor * dpi  * (1.33 + n_distinct(data@ident)/3)}
            dotplot_height <- function() {scale_factor * dpi * (1 + length(genes)/4)}
            pdf(file=NULL)
            p <- DotPlot2(data, color_scaling = "zero-one", size_scaling = "area",
                        genes.plot = genes, plot.legend = F, cols.use = c("white", "magenta"), do.return = T,
                        dot.min=0.01, dot.scale = 6 * scale_factor) +
                        theme(axis.text.y = element_text(face="bold", size=15 * scale_factor),
                              axis.text.x = element_text(size=15 * scale_factor)) + coord_flip()
            a <- dev.off()
            list(p = p, width = dotplot_width, height = dotplot_height)
        }


    SingleFeaturePlot <- function(data, gene, pt_size, embedding) {
            # only one gene
            gene_exp <- FetchData(data, gene)[,1]
            cell_order <- names(gene_exp)[gene_exp %>% order()]
            p <- FeaturePlot(object = data, features.plot = gene, cells.use = cell_order,
                        cols.use = c("grey90", "red"), pt.size = pt_size, pch.use = 19, no.axes = TRUE,
                        reduction.use = embedding, nCol=1, do.return = TRUE)
            p[[1]]
        }

    clusterplot1 <- reactiveValues()
    clusterplot2 <- reactiveValues()
    featureplot1 <- reactiveValues()
    featureplot2 <- reactiveValues()
    dotplot1 <- reactiveValues()
    dotplot2 <- reactiveValues()

    observeEvent({
        input$dataset_1
        input$dataset_2
        input$figure_scaling_check
        input$plot_width
        scale_factor()
    }, {
        # print(seurat_data_1)
        if (!is.null(seurat_data_1)) {
            output$title1 <- renderUI({h4(data1$name)})
            output$description1 <- renderUI({ sprintf("%s cells", ncells_1) })
            clusterplot1 <<- list(p = DimPlot(seurat_data_1, reduction.use = data1$embedding, do.label = T, pt.size = plot_1_ptsize(),
                                                label.size = 4, cols.use = colors_1, no.legend = T, no.axes = TRUE, do.return = TRUE,
                                                vector.friendly = (ncells_1 > rasterize_ncells), png.arguments=png.arguments),
                                    width = plot_width, height = plot_width)
            output$clusterplot1 <- renderPlot({clusterplot1$p}, width = plot_size, height = clusterplot1$height, res = dpi)
        } else {
            output$title1 <- NULL
            output$description1 <- NULL
            output$clusterplot1 <- NULL
        }

        if (!is.null(seurat_data_2)) {
            output$title2 <- renderUI({h4(data2$name)})
            output$description2 <- renderUI({ sprintf("%s cells", ncells_2) })
            clusterplot2 <<- list(p = DimPlot(seurat_data_2, reduction.use = data2$embedding, do.label = T, pt.size = plot_2_ptsize(),
                                            label.size = 4, cols.use = colors_2, no.legend = T, no.axes = TRUE, do.return = TRUE,
                                            vector.friendly = (ncells_2 > rasterize_ncells), png.arguments=png.arguments),
                                width = plot_width, height = plot_width)
            output$clusterplot2 <- renderPlot({clusterplot2$p}, width = plot_size, height = clusterplot2$height, res = dpi)
        } else {
            output$title2 <- NULL
            output$description2 <- NULL
            output$clusterplot2 <- NULL
        }
    })

    observeEvent({
        input$dataset_1
        input$dataset_2
        input$gene_symbol
        input$gene_list_submit
        input$featureplot_check
        input$figure_scaling_check
        input$plot_width
        scale_factor()
    }, {

        if (length(gene_names) > 1 || (length(gene_names) == 1 && gene_names != "")) {
            if (!is.null(seurat_data_1)) {
                dotplot1 <<- dotplot(seurat_data_1, gene_names, scale_factor())
                output$dotplot1 <- renderPlot({dotplot1$p}, width = dotplot1$width, height = dotplot1$height, res = dpi)
            }
            if (!is.null(seurat_data_2)) {
                dotplot2 <<- dotplot(seurat_data_2, gene_names, scale_factor())
                output$dotplot2 <- renderPlot({dotplot2$p}, width = dotplot2$width, height = dotplot2$height, res = dpi)
            }

            if (length(gene_names) == 1) {
                if (!is.null(seurat_data_1) && (gene_names %in% rownames(seurat_data_1@data))) {
                    featureplot1 <<- list(p=SingleFeaturePlot(seurat_data_1, gene_names, plot_1_ptsize(), data1$embedding),
                                        width=plot_width, height=plot_width)
                    output$featureplot1 <- renderPlot({featureplot1$p}, width = featureplot1$width, height = featureplot1$height, res = dpi)
                }
                if (!is.null(seurat_data_2) && (gene_names %in% rownames(seurat_data_2@data))) {
                    featureplot2 <<- list(p=SingleFeaturePlot(seurat_data_2, gene_names, plot_2_ptsize(), data2$embedding),
                                        width=plot_width, height=plot_width)
                    output$featureplot2 <- renderPlot({featureplot2$p}, width = featureplot2$width, height = featureplot2$height, res = dpi)
                }
                shinyjs::show('featureplot1')
                shinyjs::show('featureplot2')
            } else if (!always_show_featureplot()) {
                shinyjs::hide('featureplot1')
                shinyjs::hide('featureplot2')
            } else {
                shinyjs::show('featureplot1')
                shinyjs::show('featureplot2')
            }
        }

        if (is.null(seurat_data_1)) {
            output$dotplot1 <- NULL
            output$featureplot1 <- NULL
        }
        if (is.null(seurat_data_2)) {
            output$dotplot2 <- NULL
            output$featureplot2 <- NULL
        }
    })

    output$download <- downloadHandler(
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            sprintf("%s_%s.%s", input$figure_type, input$dataset, input$filetype)},
        # This function should write data to a file given to it by the argument 'file'.
        content = function(file) {
          # Write to a file specified by the 'file' argument
          plot_use <- get(sprintf("%s%s", input$figure_type, input$dataset))
          #print(plot_use())
          if (input$filetype == 'pdf') {
              ggsave(plot = plot_use$p, filename = file, width = plot_use$width() / dpi, height = plot_use$height() / dpi)
          } else if (input$filetype == 'png') {
              ggsave(plot = plot_use$p, filename = file, width = plot_use$width() / dpi, height = plot_use$height() / dpi, units = 'in', dpi = 300)
          }
        }
    )

}
