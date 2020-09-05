# Preparing input data for DRMN

The regulatory features could be context-specific or context independent. For context-independent features, we use motif networks, by finding instances of motifs in promoter regions of genes.

![alt text](example_input/motif_small.png "Motif instances in gene promoter.")

For context-specific features, we can aggregate regulatory signals (different histone modification, or chromatin accessibility).

![alt text](example_input/signal_small.png "Aggregated signals in gene promoter.")

Furthermore, we can aggregate chromatin accessibility signals (e.g. DNase/ATAC-seq) signal in motif instances.

![alt text](example_input/qmotif_small.png "Q-Motif, aggregated signal in motif instances in gene promoter.")

We used [aggregateSignal](https://github.com/Roy-lab/aggregateSignal) to aggregate signal in promoter regions or motif instances. Briefly, we use bedtools to convert bam files to count files:
```
bedtools genomecov -ibam input.bam -bg -pc > output.counts
```
and aggregateSignal programs calculate coverage in input regions. See the [repository](https://github.com/Roy-lab/aggregateSignal) for more details on how to use the program.

Furthermore, for each feature, we log transform and quantile normalize the values across cell lines/time points. Additionally, we add cell line/time point specific suffixes to gene names in order (so gene names will be specific to a cell line while regulators will be the same for all cell lines). 

In the feature file, the first line is number of regulators and number of genes in the file (tab-delim).
The rest of the file will be in 3 columns format (tab-delim):
```
feature g_cell  value
```
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
