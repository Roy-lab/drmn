# DRMN: Dynamic Regulatory Module Networks

The enclosed example script and data show you how to run DRMN and format your input files.
Script: run_example.sh
Input data: example_input/

You will need:

* A "lineage tree" file, showing the relationships between cell types.
Example: example_input/lineage_tree.txt
The format is: first column is child cell type, second column is parent cell type.
For a linear tree cell -> cell2 -> cell3, the tree is:
cell2   cell
cell3   cell2

To run DRMN on only one cell type, make your cell type the root, with no children.
NULL    cell
Example: example_input/mef_tree.txt

* A "celltype_order" file, with each cell type listed on its own line. 
This order needs to match the order of genes in the ortho_map file.
Example: example_input/celltype_order.txt

* An OGIDs file, which matches genes between cell types. 
Although DRMN doesn't have different species, it uses different versions of the gene names for each cell type. An easy way to set this up is to use $name for the gene in the first cell type, and ${name}_${cell} for the other cell types.
Example: example_input/ogids_expr_100.txt

* A config file that lists out three input filenames for each cell type.
Example: example_input/histonly_k3_config.txt
The format is (tab delimited):
cell    initial_clusterassign.txt     expression_file.txt    feature_data.txt

Each cell type-specific file needs to use the cell type-specific names for genes, matching the names in the OGIDs file:

** Cluster assignment format (tab-delim, 2 cols, no header):
g_cell    cID
You should generate the initial clusters independently per cell type using k-means, then re-order them from lowest to highest. This is necessary so that the interpretation of each state matches between cell types.
You can use the enclosed script doInitClust.m to do it in matlab.

** Expression format (tab-delim, 2 cols, has one header line, where the header for col2 is the cell type name):
Gene          cell
g_cell    log_expr

** Feature data format (tab-delim, 3 cols, no header)
feature     g_cell    value

# Usage and other parameters

Usage: ./learnDRMN celltype_order ogids_file null k lineage_tree config rand[none|yes|<int>] outputDir mode[learn|learnCV|learnCV:<int>:<int>:<int>|generate|visualize] srcnode inittype[uniform|branchlength] p_diagonal_nonleaf [const_cov(double), selfInit]

0. celltype_order as described above
1. ogids as described above
2. not used (legacy argument that we should get rid of, but haven't yet)
3. k: number of states (I think this needs to match the initial cluster assignments)
4. lineage tree as described above
5. config file as described above
6. Whether to randomize the membership of the initial clusters. We have been using "none" for DRMN.
7. output directory name, which you MUST create in advance.
8. mode:
learn: train DRMN on all data
learnCV: do 3-fold CV in serial
learnCV:fold:seed[:nfolds]: do the fold of CV corresponding to 'fold' value. Without specifying nfolds, it will do 3-fold CV and accept fold=0, 1, or 2. If you want to do more folds, give a higher value of nfolds. If you want to do CV with disjoint test sets, you need to use the same seed for each fold.
generate: not implemented
visualize: not implemented
9. srcnode: Cell type to use as the reference cell type. The gene names for this cell type will be used in some output files.
10. p_diagonal_nonleaf : The prior for genes maintaining their state assignment between two adjacent cell types. We use 0.8 by default. IT doesn't seem to affect DRMN very much so far.

Optional arguments (it will figure out based on the value type which one(s) you are using):
11. float(const_cov) : If you are getting module switches (as we did with Escarole on Randy and Morten's datasets), you may need to fix the covariances to a constant value. But, we didn't need to do this for Chronis DRMN. Try 0.2 to start if you need to.

12. selfInit: If you include this exact string, DRMN will initialize the module parameters based on the data for each cell type separately. If you don't use this, then it will initialize all module params from the srcnode cell type. If the distributions of your data differ between cell types (say you used log zero mean expression), then you will definitely need to use this. If the distributions match, then it probably doesn't matter. I use it anyway.


for Least_Dirty:
./learnDRMN order.txt ogids.txt null 3 tree.txt config.txt none out/ learnCV:0:12345 esc uniform 0.8 selfInit LEASTDIRTY 25 50
for Least_L21:
./learnDRMN order.txt ogids.txt null 3 tree.txt config.txt none out/ learnCV:0:12345 esc uniform 0.8 selfInit LEASTL21 25 50
and for original method:
./learnDRMN order.txt ogids.txt null 3 tree.txt config.txt none out/ learnCV:0:12345 esc uniform 0.8 selfInit

