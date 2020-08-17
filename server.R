library(tidyverse)
library(Seurat)
library(plotly)  # dev branch
#library(rjson)
source("./utils.R")

# -------------------------
# Shiny Single Cell Browser
# -------------------------
# Interactive visualization of single cell RNAseq datasets

ui_plot_width <- "100%"
ui_plot_height <- "auto"

plot_inch <- 4
dpi <- 125  # for web display. The saved plots have higher dpi
fixed_plot_size <- dpi * plot_inch
dotplot_height <- 500  # height in vertical layout (can expand to make it longer. )
# calc_pt_size <- function(n) {30 / n^0.5}
calc_pt_size <- function(n) {100 / n^0.63}

# png.arguments <- c(4,4,300)
# rasterize_ncells <- 150000  # This might nor work for plotly. usually use ~2000. But now vector.friendly does not work for ggplot 3.0, so I disabled this.

json_file <- rjson::fromJSON(file = './data/config.json')
json_data <- json_file$data
datasets <- 1:length(json_data)
dataset_names <- sapply(json_data, function(x) x$name)
dataset_selector <- as.list(c(datasets, "none"))
names(dataset_selector) <- c(dataset_names, "None")
# print(dataset_selector)

config <- json_file$config

startup_datasets <- list(config$startup_data1, config$startup_data2, config$startup_data3)
startup_datasets[sapply(startup_datasets, is.null)] <- "none"

if (!all(startup_datasets %in% c(datasets, "none"))) {
    warning(sprintf("`startup_data1/2/3` needs to be one of %s", paste(dataset_selector, collapse = ", ")))
    startup_datasets[!startup_datasets %in% dataset_selector] <- "none"
}

default_single_gene <- config$default_single_gene


# TODO:
# add figure legend -- done.
# Add precomputed marker genes into Seurat data
# Add Plotly interactive plots  -- done
# Click on cluster labels and show significant marker genes
# Click on dotplot gene and show tsne/umap  -- done
# Calculate gene module score with gene lists and plot on tsne/umap
# package


# -------------------------
# Read data into global environment. Share between different sessions.

read_data <- function(x) {
    # load data and metadata specified by the JSON string.
    # x: individual json string, with [name, file, clusters embedding]
    seurat_data <- readRDS(x$file)
    seurat_data <- SetAllIdent(seurat_data,  x$clusters)
    ncells <- length(seurat_data@cell.names)
    pt_size <- calc_pt_size(ncells)
    colors <- seurat_data@misc[[sprintf("%s_colors",  x$clusters)]]
    if (is.null(colors)) {
        set.seed(2)
        colors <- sample(rainbow(n_distinct(seurat_data@ident)))
    }
    genes <- rownames(seurat_data@data)
    list(
        name = x$name,
        seurat_data = seurat_data,
        ncells = ncells,
        pt_size = pt_size,
        embedding = x$embedding,
        colors = colors,
        genes = genes
    )
}

# data_list <- lapply(json_data, read_data)

logging::loginfo("loading data...")
data_list <- rep(list(NULL), length(json_data))
data_list[[1]] <- read_data(json_data[[1]])
logging::loginfo("loaded dataset #1.")


# --------------------------
# Sever sessions

function(input, output, session) {

    #print(dataset_selector)

    updateSelectizeInput(session, 'dataset_1', choices = dataset_selector, selected = startup_datasets[[1]], server = F)
    updateSelectizeInput(session, 'dataset_2', choices = dataset_selector, selected = startup_datasets[[2]], server = F)
    updateSelectizeInput(session, 'dataset_3', choices = dataset_selector, selected = startup_datasets[[3]], server = F)

    # -----------------------
    # Main UI
    output$main_panel <- renderUI({
        if (input$layout_type == 'horizontal') {
            list(
                fluidRow(
                    column(4,
                        switch(input$plot_type,
                            "plotly" = plotlyOutput("clusterplot1_plotly", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height),
                            "ggplot2" = plotOutput("clusterplot1", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height)
                        ), align="center"
                    ),
                    column(4,
                        switch(input$plot_type,
                            "plotly" = plotlyOutput("featureplot1_plotly", height = ui_plot_height, width = ui_plot_width),
                            "ggplot2" = plotOutput("featureplot1", height = ui_plot_height, width = ui_plot_width)
                        ), align="center"
                    ),
                    column(4,
                        shinydashboard::box(width = 12,
                            switch(input$plot_type,
                                "plotly" = plotlyOutput("dotplot1_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("dotplot1", height = ui_plot_height, width = ui_plot_width)
                            ), style = sprintf("height:%spx; overflow-y:scroll;",
                                            if (!is.null(session$clientData$output_dotplot1_width)) {
                                                session$clientData$output_dotplot1_width
                                            } else {
                                                session$clientData$output_dotplot1_plotly_width
                                            }
                                        )
                        ), align="center"
                    )
                ),
                fluidRow(
                    column(4,
                        switch(input$plot_type,
                            "plotly" = plotlyOutput("clusterplot2_plotly", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height),
                            "ggplot2" = plotOutput("clusterplot2", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height)
                        ), align="center"
                    ),
                    column(4,
                        switch(input$plot_type,
                            "plotly" = plotlyOutput("featureplot2_plotly", height = ui_plot_height, width = ui_plot_width),
                            "ggplot2" = plotOutput("featureplot2", height = ui_plot_height, width = ui_plot_width)
                        ), align="center"
                    ),
                    column(4,
                        shinydashboard::box(width = 12,
                            switch(input$plot_type,
                                "plotly" = plotlyOutput("dotplot2_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("dotplot2", height = ui_plot_height, width = ui_plot_width)
                            ), style = sprintf("height:%spx; overflow-y:scroll;",
                                            if (!is.null(session$clientData$output_dotplot2_width)) {
                                                session$clientData$output_dotplot2_width
                                            } else {
                                                session$clientData$output_dotplot2_plotly_width
                                            }
                                       ) #
                        ), align="center"
                    )
                ),
                fluidRow(
                    column(4,
                        switch(input$plot_type,
                            "plotly" = plotlyOutput("clusterplot3_plotly", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height),
                            "ggplot2" = plotOutput("clusterplot3", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height)
                        ), align="center"
                    ),
                    column(4,
                        switch(input$plot_type,
                            "plotly" = plotlyOutput("featureplot3_plotly", height = ui_plot_height, width = ui_plot_width),
                            "ggplot2" = plotOutput("featureplot3", height = ui_plot_height, width = ui_plot_width)
                        ), align="center"
                    ),
                    column(4,
                        shinydashboard::box(width = 12,
                            switch(input$plot_type,
                                "plotly" = plotlyOutput("dotplot3_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("dotplot3", height = ui_plot_height, width = ui_plot_width)
                            ), style = sprintf("height:%spx; overflow-y:scroll;",
                                            if (!is.null(session$clientData$output_dotplot3_width)) {
                                                session$clientData$output_dotplot3_width
                                            } else {
                                                session$clientData$output_dotplot3_plotly_width
                                            }
                                       ) #
                        ), align="center"
                    )
                )
            )
        } else if (input$layout_type == 'vertical') {
            list (
                fluidRow(
                    column(4,
                        tabsetPanel(id = "tabset_r1c1", selected = "Clusters",
                            tabPanel("Clusters",
                                switch(input$plot_type,
                                    "plotly" = plotlyOutput("clusterplot1_plotly", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height),
                                    "ggplot2" = plotOutput("clusterplot1", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height)
                                )
                            ),
                            tabPanel("Gene Expression",
                                switch(input$plot_type,
                                "plotly" = plotlyOutput("featureplot1_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("featureplot1", height = ui_plot_height, width = ui_plot_width)
                                )
                            )
                        ), align="center"
                    ),
                    column(4,
                        tabsetPanel(id = "tabset_r1c2", selected = "Clusters",
                            tabPanel("Clusters",
                                switch(input$plot_type,
                                    "plotly" = plotlyOutput("clusterplot2_plotly", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height),
                                    "ggplot2" = plotOutput("clusterplot2", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height)
                                )
                            ),
                            tabPanel("Gene Expression",
                                switch(input$plot_type,
                                    "plotly" = plotlyOutput("featureplot2_plotly", height = ui_plot_height, width = ui_plot_width),
                                    "ggplot2" = plotOutput("featureplot2", height = ui_plot_height, width = ui_plot_width)
                                )
                            )
                        ), align="center"
                    ),
                    column(4,
                        tabsetPanel(id = "tabset_r1c3", selected = "Clusters",
                            tabPanel("Clusters",
                                switch(input$plot_type,
                                    "plotly" = plotlyOutput("clusterplot3_plotly", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height),
                                    "ggplot2" = plotOutput("clusterplot3", height = ui_plot_height, width = ui_plot_width), # %>% withSpinner(proxy.height = spinner_height)
                                )
                            ),
                            tabPanel("Gene Expression",
                                switch(input$plot_type,
                                    "plotly" = plotlyOutput("featureplot3_plotly", height = ui_plot_height, width = ui_plot_width),
                                    "ggplot2" = plotOutput("featureplot3", height = ui_plot_height, width = ui_plot_width)
                                )
                            )
                        ), align="center"
                    )
                ),
                fluidRow(
                    column(4,
                        shinydashboard::box(width = 12,
                            switch(input$plot_type,
                                "plotly" = plotlyOutput("dotplot1_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("dotplot1", height = ui_plot_height, width = ui_plot_width)
                            ), style=sprintf("height:%spx; overflow-y:scroll;", dotplot_height)#sprintf("height:100%; overflow-y:scroll;", dotplot_height) #
                        ), align="center"
                    ),
                    column(4,
                        shinydashboard::box(width = 12,
                            switch(input$plot_type,
                                "plotly" = plotlyOutput("dotplot2_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("dotplot2", height = ui_plot_height, width = ui_plot_width)
                            ), style=sprintf("height:%spx; overflow-y:scroll;", dotplot_height) #
                        ), align="center"
                    ),
                    column(4,
                        shinydashboard::box(width = 12,
                            switch(input$plot_type,
                                "plotly" = plotlyOutput("dotplot3_plotly", height = ui_plot_height, width = ui_plot_width),
                                "ggplot2" = plotOutput("dotplot3", height = ui_plot_height, width = ui_plot_width)
                            ), style=sprintf("height:%spx; overflow-y:scroll;", dotplot_height) #
                        ), align="center"
                    )
                )
            )
        }
    })



    # -----------------------------
    # Load datasets.

    # gene_data <- reactiveValues()
    # No need to use reactiveValues here.
    # Note that values taken from the reactiveValues object are reactive, but the reactiveValues object itself is not

    observeEvent({
        input$dataset_1
    }, {
        if (input$dataset_1 %in% datasets) {
            current_index <- as.numeric(input$dataset_1)
            if (is.null(data_list[[current_index]])) {
                # Use <<- to modify global variable (shared across sessions)
                data_list[[current_index]] <<- read_data(json_data[[current_index]])
                logging::loginfo("loaded dataset #%s.", current_index)
            }
            data1 <<- data_list[[current_index]]
        } else {
            data1 <<- NULL
        }
        output$featureplot1 <- NULL
        output$featureplot1_plotly <- NULL
        output$dotplot1 <- NULL
        output$dotplot1_plotly <- NULL
    })

    observeEvent({
        input$dataset_2
    }, {
        if (input$dataset_2 %in% datasets) {
            current_index <- as.numeric(input$dataset_2)
            if (is.null(data_list[[current_index]])) {
                # Use <<- to modify global variable (shared across sessions)
                data_list[[current_index]] <<- read_data(json_data[[current_index]])
                logging::loginfo("loaded dataset #%s.", current_index)
            }
            data2 <<- data_list[[current_index]]
        } else {
            data2 <<- NULL
        }
        output$featureplot2 <- NULL
        output$featureplot2_plotly <- NULL
        output$dotplot2 <- NULL
        output$dotplot2_plotly <- NULL
    })

    observeEvent({
        input$dataset_3
    }, {
        if (input$dataset_3 %in% datasets) {
            current_index <- as.numeric(input$dataset_3)
            if (is.null(data_list[[current_index]])) {
                # Use <<- to modify global variable (shared across sessions)
                data_list[[current_index]] <<- read_data(json_data[[current_index]])
                logging::loginfo("loaded dataset #%s.", current_index)
            }
            data3 <<- data_list[[current_index]]
        } else {
            data3 <<- NULL
        }
        output$featureplot3 <- NULL
        output$featureplot3_plotly <- NULL
        output$dotplot3 <- NULL
        output$dotplot3_plotly <- NULL
    })

    observeEvent({
        input$dataset_1
        input$dataset_2
        input$dataset_3
    }, {
        gene_names <<- ""
        all_genes <<- as.list(c("", sort(unique(c(data1$genes, data2$genes, data3$genes)))))
        if (input$gene_symbol == "") {
            # first time loading
            if (default_single_gene == "first" && length(all_genes) > 1) gene_names <<- all_genes[[2]]
            if (default_single_gene %in% all_genes) gene_names <<- default_single_gene
        }

        #print(all_genes[1:5])
        print(gene_names)
        print(length(all_genes))
        updateSelectizeInput(session, 'gene_symbol', choices = all_genes, selected = gene_names, server = F)

        if (input$layout_type == 'vertical') {
            updateTabsetPanel(session, "tabset_r1c1", selected = "Clusters")
            updateTabsetPanel(session, "tabset_r1c2", selected = "Clusters")
            updateTabsetPanel(session, "tabset_r1c3", selected = "Clusters")
        }

    })


    # ---------------------------
    # update gene selections.

    observeEvent(input$gene_symbol, {
        gene_names <<- input$gene_symbol
        #print(gene_names)
    })

    gene_list <- ""
    observeEvent(input$gene_list_submit, {
        gene_list <<- trimws(strsplit(input$gene_list, '\n')[[1]])
        if (length(gene_list) == 1 && gene_list != "") {
            # simply trigger the gene_symbol box
            if (gene_list %in% all_genes) {
                updateSelectizeInput(session, 'gene_symbol', selected = gene_list, server = F)
            }
        }
        #print(gene_list)
    })

    # ----------------------------
    # Manual-scaling figure size

    observeEvent({  # this may be redundant
        input$auto_scaling_check
        input$plot_width
    }, {
        if(!input$auto_scaling_check) {
            plot_width <<- function() {input$plot_width}
        }
    })


    # ----------------------------
    # Auto-scaling figure size (check all figures)

    all_ggplots <- c("clusterplot1", "clusterplot2", "clusterplot3",
                    "featureplot1", "featureplot2", "featureplot3",
                    "dotplot1", "dotplot2", "dotplot3")

    for (plot_name in c(all_ggplots)) {
        #print(plot_name)
        observe({
            if (input$auto_scaling_check) {
                plot_width <<- function() {
                    gg_width <- sprintf("output_%s_width", plot_name)
                    pl_width <- sprintf("output_%s_plotly_width", plot_name)
                    switch(
                        input$plot_type,
                        "ggplot2" = if (is.null(session$clientData[[gg_width]]) && !is.null(session$clientData[[pl_width]]))
                                        {session$clientData[[pl_width]]}
                                    else {session$clientData[[gg_width]]},
                        "plotly" = if (is.null(session$clientData[[pl_width]]) && !is.null(session$clientData[[gg_width]]))
                                        {session$clientData[[gg_width]]}
                                    else {session$clientData[[pl_width]]}
                    )
                }
            } else {
                plot_width <<- function() {input$plot_width}
            }
            #print(plot_width())
            #print(scale_factor())
        })
    }

    # older implementation only checking the size of the first plot
    # observeEvent({
    #     input$auto_scaling_check
    #     input$plot_width
    # }, {
    #     if(input$auto_scaling_check) {
    #         plot_width <<- function() {
    #                 session$clientData$output_clusterplot1_width
    #         }
    #     } else {
    #         plot_width <<- function() {input$plot_width}
    #     }
    #     #print(plot_width())
    # })

    # ----------------------
    # calculate scale factor
    scale_factor <- function() {plot_width() / dpi / plot_inch}

    # ----------------------
    # Functions for plotting

    ClusterPlotScaled <- function(data, embedding, title_use, colors_use, n_cells, raw_ptsize, scale_factor) {
        scaled_ptsize <- function() {scale_factor ^ 2 * raw_ptsize}
        p <- DimPlot(data, reduction.use = embedding, do.label = T, pt.size = scaled_ptsize(),
                label.size = 4 * scale_factor, cols.use = colors_use, no.legend = T, no.axes = TRUE, do.return = TRUE,
                vector.friendly = FALSE)
        p <- p +
            labs(title = title_use) +
            theme(plot.title = element_text(hjust = 0.5, face="bold", size = 14 * scale_factor))
        p
    }

    SingleFeaturePlotScaled <- function(data, gene, embedding, raw_ptsize, scale_factor) {
        # only one gene
        scaled_ptsize <- function() {scale_factor ^ 2 * raw_ptsize}
        gene_exp <- FetchData(data, gene)[,1]
        cell_order <- names(gene_exp)[gene_exp %>% order()]
        plist <- FeaturePlot(object = data, features.plot = gene, cells.use = cell_order,
                    cols.use = c("grey90", "red"), pt.size = scaled_ptsize(), pch.use = 19, no.axes = TRUE,
                    reduction.use = embedding, nCol=1, do.return = TRUE)
        p <- plist[[1]] + theme(plot.title = element_text(size = 14 * scale_factor))
        p
    }

    DotPlotScaled <- function(data, genes, scale_factor) {
            print(scale_factor)
            #print(scale_factor * dpi)
            scaled_width <- function() {scale_factor * dpi * 0.633 * (1.8 + n_distinct(data@ident)/3.5)}
            scaled_height <- function() {scale_factor * dpi * 0.633 * (1.2 + length(unique(genes))/4) + 120}
            scaled_ptsize <- function() {scale_factor ^ 2 * 3.8}
            pdf(file=NULL)
            p <- DotPlot2(data, color_scaling = "zero-one", size_scaling = "area",
                        genes.plot = genes, legend.position = "bottom", cols.use = c("white", "magenta"), do.return = T,
                        dot.min=0.01, dot.scale = scaled_ptsize(), axis.label.size = 7 * scale_factor,
                        horizontal = FALSE) +
                        theme(legend.direction='vertical',
                              legend.text=element_text(size=7 * scale_factor),
                              legend.title=element_text(size=10 * scale_factor),
                              legend.key.size = unit(0.5, 'lines'))
            a <- dev.off()
            list(p = p, width = scaled_width, height = scaled_height)
        }


    # ---------------------------
    #  Cluster Plots

    observeEvent({
        input$dataset_1
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
        #session$clientData$output_clusterplot1_width
        #session$clientData$output_clusterplot1_plotly_width
    }, {

        if (!is.null(data1$seurat_data)) {
            #output$title1 <- renderUI({h4(data1$name)})
            # output$description1 <- renderUI({ sprintf("%s cells", data1$ncells) })
            clusterplot1 <<- list(
                p = ClusterPlotScaled(data1$seurat_data, data1$embedding, data1$name,
                            data1$colors, data1$ncells, data1$pt_size, scale_factor()),
                width = plot_width,
                height = plot_width
            )

            if (input$plot_type == 'ggplot2') {
                output$clusterplot1 <- renderPlot({
                    clusterplot1$p
                }, width = clusterplot1$width, height = clusterplot1$height, res = dpi)
                output$clusterplot1_plotly <- NULL
            } else if (input$plot_type == 'plotly') {
                output$clusterplot1_plotly <- renderPlotly({
                    ggplotly(clusterplot1$p, width = clusterplot1$width(), height = clusterplot1$height(),
                            tooltip = c("factor(x = ident)"),
                            source = "clusterplot1_plotly")
                })
                output$clusterplot1 <- NULL
            }
        } else {
            output$title1 <- NULL
            # output$description1 <- NULL
            output$clusterplot1 <- NULL
            output$clusterplot1_plotly <- NULL
        }
    })

    observeEvent({
        input$dataset_2
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
        #session$clientData$output_clusterplot2_width
        #session$clientData$output_clusterplot2_plotly_width
    }, {
        if (!is.null(data2$seurat_data)) {
            #output$title2 <- renderUI({h4(data2$name)})
            # output$description2 <- renderUI({ sprintf("%s cells", data2$ncells) })
            clusterplot2 <<- list(
                p = ClusterPlotScaled(data2$seurat_data, data2$embedding, data2$name,
                            data2$colors, data2$ncells, data2$pt_size, scale_factor()),
                width = plot_width,
                height = plot_width
            )

            if (input$plot_type == 'ggplot2') {
                output$clusterplot2 <- renderPlot({
                    clusterplot2$p
                }, width = clusterplot2$width, height = clusterplot2$height, res = dpi)
                output$clusterplot2_plotly <- NULL
            } else if (input$plot_type == 'plotly') {
                output$clusterplot2_plotly <- renderPlotly({
                    ggplotly(clusterplot2$p, width = clusterplot2$width(), height = clusterplot2$height(),
                            tooltip = c("factor(x = ident)"),
                            source = "clusterplot2_plotly")
                })
                output$clusterplot2 <- NULL
            }
        } else {
            output$title2 <- NULL
            # output$description2 <- NULL
            output$clusterplot2 <- NULL
            output$clusterplot2_plotly <- NULL
        }
    })

    observeEvent({
        input$dataset_3
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
        #input$dimension[1]
    }, {
        if (!is.null(data3$seurat_data)) {
            #output$title3 <- renderUI({h4(data3$name)})
            # output$description3 <- renderUI({ sprintf("%s cells", data3$ncells) })
            clusterplot3 <<- list(
                p = ClusterPlotScaled(data3$seurat_data, data3$embedding, data3$name,
                            data3$colors, data3$ncells, data3$pt_size, scale_factor()),
                width = plot_width,
                height = plot_width
            )

            if (input$plot_type == 'ggplot2') {
                output$clusterplot3 <- renderPlot({
                    clusterplot3$p
                }, width = clusterplot3$width, height = clusterplot3$height, res = dpi)
                output$clusterplot3_plotly <- NULL
            } else if (input$plot_type == 'plotly') {
                output$clusterplot3_plotly <- renderPlotly({
                    ggplotly(clusterplot3$p, width = clusterplot3$width(), height = clusterplot3$height(),
                            tooltip = c("factor(x = ident)"),
                            source = "clusterplot3_plotly")
                })
                output$clusterplot3 <- NULL
            }
        } else {
            output$title3 <- NULL
            # output$description3 <- NULL
            output$clusterplot3 <- NULL
            output$clusterplot3_plotly <- NULL
        }
    })

    # ---------------------------
    #  Gene Plots

    observeEvent({
        input$dataset_1
        input$gene_symbol
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
    }, {
        if (length(gene_names) == 1 && gene_names != "") {
            if (!is.null(data1$seurat_data) && (gene_names %in% data1$genes)) {
                featureplot1 <<- list(
                    p=SingleFeaturePlotScaled(data1$seurat_data, gene_names, data1$embedding, data1$pt_size, scale_factor()),
                    width=plot_width,
                    height=plot_width
                )

                if (input$plot_type == 'ggplot2') {
                    output$featureplot1 <- renderPlot({
                        featureplot1$p
                    }, width = featureplot1$width, height = featureplot1$height, res = dpi)
                    output$featureplot1_plotly <- NULL
                } else if (input$plot_type == 'plotly') {
                    output$featureplot1_plotly <- renderPlotly({
                        ggplotly(featureplot1$p, width = featureplot1$width(), height = featureplot1$height(),
                                tooltip = c("gene"),
                                source = "featureplot1_plotly")
                    })
                    output$featureplot1 <- NULL
                }

            } else {
                output$featureplot1 <- NULL
                output$featureplot1_plotly <- NULL
            }
        }
    })


    observeEvent({
        input$dataset_2
        input$gene_symbol
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
    }, {
        if (length(gene_names) == 1 && gene_names != "") {
            if (!is.null(data2$seurat_data) && (gene_names %in% data2$genes)) {
                featureplot2 <<- list(
                    p=SingleFeaturePlotScaled(data2$seurat_data, gene_names, data2$embedding, data2$pt_size, scale_factor()),
                    width=plot_width,
                    height=plot_width
                )

                if (input$plot_type == 'ggplot2') {
                    output$featureplot2 <- renderPlot({
                        featureplot2$p
                    }, width = featureplot2$width, height = featureplot2$height, res = dpi)
                    output$featureplot2_plotly <- NULL
                } else if (input$plot_type == 'plotly') {
                    output$featureplot2_plotly <- renderPlotly({
                        ggplotly(featureplot2$p, width = featureplot2$width(), height = featureplot2$height(),
                                tooltip = c("gene"),
                                source = "featureplot2_plotly")
                    })
                    output$featureplot2 <- NULL
                }
            } else {
                output$featureplot2 <- NULL
                output$featureplot2_plotly <- NULL
            }
        }
    })

    observeEvent({
        input$dataset_3
        input$gene_symbol
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
    }, {
        if (length(gene_names) == 1 && gene_names != "") {
            if (!is.null(data3$seurat_data) && (gene_names %in% data3$genes)) {
                featureplot3 <<- list(
                    p=SingleFeaturePlotScaled(data3$seurat_data, gene_names, data3$embedding, data3$pt_size, scale_factor()),
                    width=plot_width,
                    height=plot_width
                )

                if (input$plot_type == 'ggplot2') {
                    output$featureplot3 <- renderPlot({
                        featureplot3$p
                    }, width = featureplot3$width, height = featureplot3$height, res = dpi)
                    output$featureplot3_plotly <- NULL
                } else if (input$plot_type == 'plotly') {
                    output$featureplot3_plotly <- renderPlotly({
                        ggplotly(featureplot3$p, width = featureplot3$width(), height = featureplot3$height(),
                                tooltip = c("gene"),
                                source = "featureplot3_plotly")
                    })
                    output$featureplot3 <- NULL
                }
            } else {
                output$featureplot3 <- NULL
                output$featureplot3_plotly <- NULL
            }
        }
    })


    # ------------------------
    # Switch from cluster plots to gene plots when a gene is selected

    observeEvent({
        input$gene_symbol
    }, {
        if (length(gene_names) == 1 && gene_names != "") {
            if (input$layout_type == 'vertical') {
                updateTabsetPanel(session, "tabset_r1c1", selected = "Gene Expression")
                updateTabsetPanel(session, "tabset_r1c2", selected = "Gene Expression")
                updateTabsetPanel(session, "tabset_r1c3", selected = "Gene Expression")
            }
        }
    })

    # ------------------------
    # Dot Plots

    observeEvent({
        input$dataset_1
        input$gene_list_submit
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
    }, {

        if (length(gene_list) > 1 || (length(gene_list) == 1 && gene_list != "")) {
            if (!is.null(data1$seurat_data)) {
                dotplot1 <<- DotPlotScaled(data1$seurat_data, gene_list, scale_factor())

                if (input$plot_type == 'ggplot2') {
                    output$dotplot1 <- renderPlot({
                        dotplot1$p
                    }, width = dotplot1$width, height = dotplot1$height, res = dpi)
                    output$dotplot1_plotly <- NULL
                } else if (input$plot_type == 'plotly') {
                    output$dotplot1_plotly <- renderPlotly({
                        ggplotly(dotplot1$p, width = dotplot1$width() + scale_factor() * dpi * 0.55,
                                height = dotplot1$height() - 90, source = "dotplot1_plotly")
                    })
                    output$dotplot1 <- NULL
                }
            }
        }
        if (is.null(data1$seurat_data)) {
            output$dotplot1 <- NULL
            output$dotplot1_plotly <- NULL
        }
    })

    observeEvent({
        input$dataset_2
        input$gene_list_submit
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
    }, {

        if (length(gene_list) > 1 || (length(gene_list) == 1 && gene_list != "")) {
            if (!is.null(data2$seurat_data)) {
                dotplot2 <<- DotPlotScaled(data2$seurat_data, gene_list, scale_factor())

                if (input$plot_type == 'ggplot2') {
                    output$dotplot2 <- renderPlot({
                        dotplot2$p
                    }, width = dotplot2$width, height = dotplot2$height, res = dpi)
                    output$dotplot2_plotly <- NULL
                } else if (input$plot_type == 'plotly') {
                    output$dotplot2_plotly <- renderPlotly({
                        ggplotly(dotplot2$p, width = dotplot2$width() + scale_factor() * dpi * 0.55,
                                height = dotplot2$height() - 90, source = "dotplot2_plotly")
                    })
                    output$dotplot2 <- NULL
                }
            }
        }
        if (is.null(data2$seurat_data)) {
            output$dotplot2 <- NULL
            output$dotplot2_plotly <- NULL
        }
    })

    observeEvent({
        input$dataset_3
        input$gene_list_submit
        input$auto_scaling_check
        input$plot_width
        input$plot_type
        scale_factor()
    }, {

        if (length(gene_list) > 1 || (length(gene_list) == 1 && gene_list != "")) {
            if (!is.null(data3$seurat_data)) {
                dotplot3 <<- DotPlotScaled(data3$seurat_data, gene_list, scale_factor())

                if (input$plot_type == 'ggplot2') {
                    output$dotplot3 <- renderPlot({
                        dotplot3$p
                    }, width = dotplot3$width, height = dotplot3$height, res = dpi)
                    output$dotplot3_plotly <- NULL
                } else if (input$plot_type == 'plotly') {
                    output$dotplot3_plotly <- renderPlotly({
                        ggplotly(dotplot3$p, width = dotplot3$width() + scale_factor() * dpi * 0.55,
                                height = dotplot3$height() - 90, source = "dotplot3_plotly")
                    })
                    output$dotplot3 <- NULL
                }
            }
        }
        if (is.null(data3$seurat_data)) {
            output$dotplot3 <- NULL
            output$dotplot3_plotly <- NULL
        }
    })

    # ----------------------------
    # Click on dot plot and show the gene on t-SNE/UMAP
    observe({
        # Get subset based on selection
        click_dotplot_1 <- event_data("plotly_click", source = "dotplot1_plotly")
        click_dotplot_2 <- event_data("plotly_click", source = "dotplot2_plotly")
        click_dotplot_3 <- event_data("plotly_click", source = "dotplot3_plotly")

        # If NULL dont do anything
        click_data <- NULL
        if (!is.null(click_dotplot_1)) click_data <- click_dotplot_1
        if (!is.null(click_dotplot_2)) click_data <- click_dotplot_2
        if (!is.null(click_dotplot_3)) click_data <- click_dotplot_3

        if (is.null(click_data)) return(NULL)
        print(click_data)
        #print(click_dotplot_2)
        gene_selected <- rev(unique(gene_list))[click_data$y]  # need to be unique
        if (gene_selected %in% all_genes) {
            updateSelectizeInput(session, 'gene_symbol', selected = gene_selected, server = F)
        }
    })

    # ----------------------------
    # Download figures

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
