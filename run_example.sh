#!/bin/bash
# Shows some running examples for DRMN.

out=example_out
mkdir -p $out

# Number of clusters
k=3

### Run DRMN with all 3 cell types. Learn one model from all data. ###
./learnDRMN example_input/celltype_order.txt example_input/ogids_expr_100.txt null $k example_input/lineage_tree.txt example_input/histonly_k${k}_config.txt none $out learn ips uniform 0.8 selfInit

exit

##  Run fold 0 out of {0,1,2} of cross-validation
seed=12345
./learnDRMN example_input/celltype_order.txt example_input/ogids_expr_100.txt null $k example_input/lineage_tree.txt example_input/histonly_k${k}_config.txt none $out learnCV:0:${seed} ips uniform 0.8 selfInit

exit

### Run DRMN on one cell type only.
mycell=mef
./learnDRMN example_input/celltype_order.txt example_input/ogids_expr_100.txt null $k example_input/${cell}_tree.txt example_input/histonly_k${k}_config.txt none $out learn $mycell uniform 0.8 selfInit


