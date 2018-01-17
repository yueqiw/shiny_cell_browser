library(tidyverse)
library(Seurat)
library(rjson)

json_data <- fromJSON(file = './config.json')
umi_seurat_1 <- readRDS(json_data$file_1)
umi_seurat_2 <- readRDS(json_data$file_2)

plot_size <- 500

function(input, output, session) {

    gene_name <- eventReactive(input$submit, {
        trimws(input$gene_symbol)
    })

    output$value_1 <- renderPrint({
        if (gene_name() %in% umi_seurat_1@data %>% rownames()) {
            ""
        } else {
            sprintf("%s is not detected ( < 2 umi's in  < 2 cells) in scSeq dataset.", gene_name())
        }
    })

    output$exp_level_1 <- renderTable({
        umi_seurat_1@hvg.info[gene_name(),] %>%
            rownames_to_column(var = 'gene') %>%
            mutate(gene.sum = sum(umi_seurat_1@data[gene_name(),]),
                   n_cells = sum(umi_seurat_1@data[gene_name(),]>0))

    })

    output$plot1 <- renderPlot({
        p <- FeaturePlot(object = umi_seurat_1, features.plot = c(gene_name()),
                    cols.use = c("lightgrey", "blue"), pt.size = 0.8,
                    reduction.use = "tsne_p30_s4", nCol=1, do.return = TRUE)
        p[[1]]
    }, width = plot_size, height = plot_size)

    output$plot2 <- renderPlot({
        a <- gene_name()
        #p <- DimPlot(umi_seurat, reduction.use = "tsne", do.label = T, pt.size = 0.8, label.size = 5, colors.use = rainbow(16), no.legend = T)
        p <- DimPlot(umi_seurat_1, reduction.use = "tsne_p30_s4", do.label = T, pt.size = 0.8, label.size = 6, colors.use = rainbow(n_distinct(umi_seurat_1@ident)), no.legend = T)
        p
    }, width = plot_size, height = plot_size)

    output$value_2 <- renderPrint({
        if (gene_name() %in% umi_seurat_2@data %>% rownames()) {
            ""
        } else {
            sprintf("%s is not detected ( < 2 umi's in  < 2 cells) in scSeq dataset.", gene_name())
        }
    })

    output$exp_level_2 <- renderTable({
        umi_seurat_2@hvg.info[gene_name(),] %>%
            rownames_to_column(var = 'gene') %>%
            mutate(gene.sum = sum(umi_seurat_2@data[gene_name(),]),
                   n_cells = sum(umi_seurat_2@data[gene_name(),]>0))

    })

    output$plot3 <- renderPlot({
        p <- FeaturePlot(object = umi_seurat_2, features.plot = c(gene_name()),
                    cols.use = c("lightgrey", "blue"), pt.size = 0.2,
                    reduction.use = "tsne_p50_s2", nCol=1, do.return = TRUE)
        p[[1]]
    }, width = plot_size, height = plot_size)

    output$plot4 <- renderPlot({
        a <- gene_name()
        #p <- DimPlot(umi_seurat, reduction.use = "tsne", do.label = T, pt.size = 0.8, label.size = 5, colors.use = rainbow(16), no.legend = T)
        p <- DimPlot(umi_seurat_2, reduction.use = "tsne_p50_s2", do.label = T, pt.size = 0.2, label.size = 6, colors.use = rainbow(n_distinct(umi_seurat_2@ident)), no.legend = T)
        p
    }, width = plot_size, height = plot_size)

}
