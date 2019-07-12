## Updates
08-15-2018
  - Full application rewrite for speed and simplicity in viewing only 1 dataset at a time.  Table added for differential expression navigation.  Download features still need to be re-implemented.

07-17-2018
  - Moved data loading to the global environment rather than per user session. Added more config options. 
  
07-15-2018
  - Visualize up to 3 datasets simultaneously and interactively in the same window. 
  - Both horizontal and vertical layouts are available. 

07-14-2018

  - Add Plotly interactive mode. Click on a gene in the dot plot to show it's distribution on t-SNE/UMAP.
  - The interactive mode requires the new ggplot2 `v3.0.0` and the dev branch of plotly. 
  - Non-interactive mode works for both ggplot2 `v2.2.1` and `v3.0.0`

07-13-2018

  - Visualize two datasets simultaneously. Can easily switch beteen more datasets from the dropdown menu.
  - Plot cluster distribution on UMAP/t-SNE plots.
  - Plot the expression pattern of individual marker genes on UMAP/t-SNE embeddings.
  - Plot cluster-averaged expression of marker gene lists using dot plots.
  - Automatic resizing/scaling of figures to fit different browser windows and screen resolutions.
  - Export publication-quality figures in PDF and PNG formats. (recommend using manual resizing for consistency)
  - Specify pre-analyzed datasets in the JSON config file (see `example_config.json`).
  - Currently support [Seurat](https://github.com/satijalab/seurat) data format.
