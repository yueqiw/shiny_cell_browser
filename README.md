
## Shiny Single Cell Browser

Open-source interactive visualization of single cell RNAseq datasets. Features include:

  - Visualize cluster distribution, marker gene expression and cluster-averaged expression of gene lists. 
  - Interactive visualization: clicking on individual genes in the gene list plot shows their expression on t-SNE/UMAP plot. 
  - Specify pre-analyzed datasets ([Seurat](https://github.com/satijalab/seurat) format) in the JSON config file for automatic data loading.

## Setting up and launch the App
  
  - Download the App, `git clone https://github.com/yueqiw/shiny_cell_browser.git`
  - Install dependencies as listed [below](#dependencies). 
  - To update, `cd shiny_cell_browser` then `git pull`
  - Store Seurat data objects as `.rds` files 
  - Optionally, store cluster colors as a vector in `seurat_data@misc[[sprintf("%s_colors", cluster_name)]]`
  - Specify the file paths and clustering results by creating `data/config.json` file. Follow the example in [`data/example_config.json`](data/example_config.json). The App will automatically open the files specified in `startup_data1/2/3`. 
  - To launch Single Cell Browser locally, run the following code.  
  ```
  cd shiny_cell_browser
  R -e "shiny::runApp('./', port=1234)
  ## or store the lunch script in run_app.sh and run the following
  ./run_app.sh 
  ```
  - This should launch the web browser at `http://127.0.0.1:1234/`
  - If you want other computers in the local network to access the web app, run `R -e "shiny::runApp('./', host='0.0.0.0' port=1234)`. Then visit `your-ip-address:1234`
  
Example `config.json` file: 

```
{
    "data": [
        {
            "name": "My favorite sample",
            "file": "path/to/seurat/data.rds",
            "clusters": "res.1",
            "embedding": "umap_1"
        },
        {
            "name": "My 2nd favorite sample",
            "file": "path/to/seurat/data.rds",
            "clusters": "res.1",
            "embedding": "tsne_1"
        }
    ],
    "config": {
        "startup_data1": 1,
        "startup_data2": 2,
        "startup_data3": "none",
        "default_layout": "horizontal",
        "default_viz_mode": "interactive",
        "default_single_gene": "none"
    }
}
```

### Dependencies

Check the Dockerfile.
  
## Updates
08-15-2018
  - Full application rewrite for speed and simplicity.  Table added for differential expression navigation.  Download features still need to be re-implemented.

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



