# Other input files for DRMN

In order to run DRMN, we need: 

* A "lineage tree" file, showing the relationships between cell types.
Example: example_input/tree.txt
The format is: first column is child cell type, second column is parent cell type.
For a linear tree cell -> cell2 -> cell3, the tree is:
```
cell2   cell
cell3   cell2
```

   * To run DRMN on only one cell type, make your cell type the root, with no children.
```
NULL    cell
```

* A "celltype_order" file, with each cell type listed on its own line. This order needs to match the order of genes in the ortho_map file. Example: example_input/order.txt

* An OGIDs file, which matches genes between cell types. Although DRMN doesn't have different species, it uses different versions of the gene names for each cell type. An easy way to set this up is to use $name for the gene in the first cell type, and ${name}\_${cell} for the other cell types.
Example: example_input/ogids.txt

* A config file that lists out three input filenames for each cell type. Example: atac_qmotif_chromatin_9marks_k7/config.txt
The first line should be the location of the files specified in the rest of the config file.
The format is (tab delimited), each cell type-specific file needs to use the cell type-specific names for genes, matching the names in the OGIDs file:
```
cell    initial_clusterassign.txt     expression_file.txt    feature_data.txt
```

   * Cluster assignment format (tab-delim, 2 cols, no header):
```
g_cell    cID
```
   You should generate the initial clusters independently per cell type using k-means, then re-order them from lowest to highest. This is necessary so that the interpretation of each state matches between cell types. You can use the enclosed script doInitClust.m to do it in matlab.

   * Expression format (tab-delim, 2 cols, has one header line, where the header for col2 is the cell type name):
```
Gene          cell
g_cell    log_expr
```
