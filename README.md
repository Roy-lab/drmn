# DRMN: Dynamic Regulatory Module Networks

[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

Dynamic Regulatory Module Networks (DRMN) is a computational framework to infer context-specific regulatory network for cell lineage or a time course. It can incorporate context-specific features (such as histone modification) and context-independent features (such as motif networks). In order to handle small number of expression samples (1 per time point/cell line) we first cluster genes into groups of co-expressed genes, and then infer regulatory program for each module, with additional constraint that time points and cell lines that are close to each other should have similar regulatory program and modules.

![alt text](example_input/drmn_overview.png "Overview of DRMN. Given a cell lineage (or a time course), and context-specific and context-independent features for each cell line, DRMN infers modules of coexpressed genes in each cell line and infers a regulatory program for each module. DRMN allows for change in module assignment of genes across cell lines based on similarity of cell lines and changes in expression of genes.")


For data pre-processing steps needed for feature generation, see [feature generation](dataprocessing.md). For explanation of different input files, see [other input files for DRMN](otherinputs.md).

The enclosed example script[run_example.sh](run_example.sh) shows how to run DRMN on an example input dataset[example_input].


# Usage and other parameters

```
Usage: ./learnDRMN celltype_order ogids_file null k lineage_tree config rand[none|yes|<int>] outputDir mode[learn|learnCV|learnCV:<int>:<int>:<int>|generate|visualize] srcnode inittype[uniform|branchlength] p_diagonal_nonleaf [selfInit leasttype[LEASTFUSED|GREEDY] p1 p2 p3]
```

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

11. selfInit (optional): If you include this exact string, DRMN will initialize the module parameters based on the data for each cell type separately. If you don't use this, then it will initialize all module params from the srcnode cell type. If the distributions of your data differ between cell types (say you used log zero mean expression), then you will definitely need to use this. If the distributions match, then it probably doesn't matter. I use it anyway.

12. leasttype(optional): The multitask regression algorithm, LEASTFUSED for fused lasso, and if it is GREEDY or not specified, it ran the greedy hill climbing algorithm.
   * fused lasso has 3 hyper parameters (sparsity, fused penalty, and group penalty).

for LEASTFUSED:
```
./learnDRMN order.txt ogids.txt null 3 tree.txt config.txt none out/ learnCV:0:12345 esc uniform 0.8 selfInit LEASTFUSED 25 50 50
```
and for GREEDY method:
```
./learnDRMN order.txt ogids.txt null 3 tree.txt config.txt none out/ learnCV:0:12345 esc uniform 0.8 selfInit
```
