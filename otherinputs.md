# DRMN input files

In order to run DRMN, we need the following input files: 

1. A "lineage tree" file, showing the relationships between cell types.
Example: example_input/tree.txt
The format is: first column is child cell type, second column is parent cell type.
For a linear tree cell -> cell2 -> cell3, the tree is. [See example tree](../blob/master/example_input/tree.txt):
```
cell2   cell
cell3   cell2
```

**_Note_** To run DRMN on only one cell type, make your cell type the root, with no children.
```
NULL    cell
```

2. An OGIDs file, which matches genes between cell types. DRMN uses different versions of the gene names for each cell type. An easy way to set this up is to use $name for the gene in the first cell type, and ${name}\_${cell} for the other cell types.
[See example OGIDs file](example_input/ogids.txt)

3. A celltype_order file, with each cell type listed on its own line. This order needs to match the order of genes in the OGIDs file. [See example](../blob/master/example_input/order.txt)

4. A config file that lists out three input filenames for each cell type. [See example](../blob/master/example_input/atac_qmotif_chromatin_9marks_k7/config.txt)
The first line should be the location of the files specified in the rest of the config file. All subsequent lines have four tab-delimited columns: cell-type name, initial clusterassignment file, the expression file and the feature file. These files are described below. 
```
cell    initial_clusterassign.txt     expression_file.txt    feature_data.txt
```

5. **initial clusterassign file**. This is a two-column tab-delimited file, first column is the cell-type specific gene name, second column is the cluster assignment. The gene name needs to match the corresponding entries in the OGIDs file::
```
g_cell    cID
```
The initial cluster assignment files are generated independently per cell type using k-means with the same `k`, then re-ordered such that the lowest cluster ID has the lowest expression and the highest cluster ID has the highest expression. This is necessary so that the interpretation of each state matches between cell types. We provide a matlab script [doInitClust.m](doInitClust.m) to do this in matlab.

6. **expression file**. This is a two-column tab-delimited file. The first line is a header, the first column is the gene name, and the second column is the expression value. The second column of the header line is the cell type name:
```
Gene          cell
g_cell    log_expr
```

7. **feature file**. 
The feature file is a tab-delmited file. The first line has two entries: first is the number of features and the second is the number of genes. The remaining lines are three column, tab-delimited lines. Each line has entries as follows [See example in tar file](example_input/atac_qmotif_chromatin_9marks_k7.tar.gz):

```
feature1 g_cell 1.2
feature1 h_cell 0.3
.
.
feature40 g_cell 5
```
