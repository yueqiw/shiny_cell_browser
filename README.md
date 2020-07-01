
## Shiny Single Cell Browser

Interactive visualization of single cell RNAseq datasets. 

  - Visualize cluster distribution, marker gene expression, and cluster-averaged expression of lists of genes. 
  - Select or click on a gene to show its expression on t-SNE/UMAP plots, select a cluster to show its marker genes.
  - Specify pre-analyzed datasets ([Seurat 2 or 3](https://github.com/satijalab/seurat) format) in the JSON config file as data source. Easily switch between differen datasets.

## Setting up and launch the App
  
  - Download the App, `git clone https://github.com/yueqiw/shiny_cell_browser.git`.
  - The master branch currently supports Seurat2. For Seurat3 datasets, check out the [seurat3](https://github.com/yueqiw/shiny_cell_browser/tree/seurat3) branch.
  - Install dependencies as listed [below](#dependencies).
  - Prepare data
    - Store Seurat data objects as `.rds` files.
    - Store marker gene differential expression table in `.csv` files (column names must contain `gene` and `cluster`).
    - Optionally, store cluster colors as a vector in `seurat_data@misc[[sprintf("%s_colors", cluster_name)]]`.
  - Specify the file paths and parameters by creating a `data/config.json` file. Follow the example in [`data/example_config.json`](data/example_config.json). The App will load the files on startup. 
  - To launch Single Cell Browser locally, run the following code.  
  ```
  cd shiny_cell_browser
  R -e "shiny::runApp('./', port=1234)
  ## or store the lunch script in run_app.sh and run the following
  ./run_app.sh 
  ```
  - This should launch the web browser at `http://127.0.0.1:1234/`. For other computers in the local network to access the web app, run `R -e "shiny::runApp('./', host='0.0.0.0' port=1234)`. Then visit `your-ip-address:1234`.
  - The App can be deployed on a web server using Docker or [shinyapps.io](https://www.shinyapps.io).
  
Example `config.json` file: 

```
{
    "data": [
        {
            "name": "My 1st sample",
            "file": "path/to/seurat/data.rds",
            "cluster": "res.1",
            "embedding": "umap",
            "diff_ex_cluster": "res.1", 
	    "diff_ex_file": "path/to/differential_expression/markers.csv"
        },
        {
            "name": "My 2nd sample",
            "file": "path/to/seurat/data2.rds",
            "cluster": "res.1_rename1",
            "embedding": "tsne",
            "diff_ex_cluster": "res.1", 
            "diff_ex_file": "path/to/differential_expression/markers.csv",

            "cluster_name_mapping": {
                "C1": "Neurons",
                "C2": "Astrocytes",
                "C3": "Neural Progenitors",
                "Note": "cluster_name_mapping is optional"
            },
            "pt_size": 2,
            "font_scale": 0.75
        }
    ],
    "config": {
        "ui_title": "Single Cell Browser",
        "title_link_text": "Optional subtitle (e.g. your lab)",
        "title_link_url": "http://optional-link-to-your-lab.com"
    }

}

```

### Dependencies

Check the Dockerfile.
  
## Updates

see [updates.md](UPDATES.md)



