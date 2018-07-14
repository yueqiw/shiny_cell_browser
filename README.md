
## Shiny Single Cell Browser

Interactive visualization of single cell RNAseq datasets

07-14-2018

  - Add Plotly interactive mode. Click on a gene in the dot plot to show it's distribution on t-SNE/UMAP.
  - The interactive mode require the new ggplot2 3.0.0 and dev branch of plotly. `devtools::install_github("ropensci/plotly")`
  - Non-interactive mode works for both ggplot2 2.2.1 and 3.0.0 

07-13-2018

  - Visualize two datasets simultaneously. Can easily switch beteen more datasets from dropdown menu.
  - Visualize cluster distribution on UMAP/t-SNE plots.
  - Plot the expression pattern of individual marker genes on UMAP/t-SNE embeddings.
  - Plot cluster-averaged expression of gene lists using dot plots.
  - Automatic resizing/scaling of figures to fit different browser windows and screen resolutions.
  - Export publication-quality figures in PDF and PNG formats. (use manual scaling for consistency)
  - Specify pre-analyzed datasets in the JSON config file (see `example_config.json`).
  - Currently support Seurat format.
  
### Dependencies

  - R==3.4.3 or 3.5.0
  - Seurat==2.3.2
  - rjson==0.2.15
  - shiny==1.1.0
  - shinyjs==1.0
  - tidyverse==1.2.1
  - ggplot2==2.2.1

For plotly interactive mode, install ggplot2 3.0.0 and developmental branch of plotly
  - ggplot2==3.0.0
  - plotly==4.7.1
