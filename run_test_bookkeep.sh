#!/bin/bash
# Compare bookkeeping version to original

orig=test_out_orig
bk=test_out_bookkeep

mkdir -p $orig $bk

#echo gdb ./learnDRMN_orig 
#echo test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none $orig learnCV:0:123 ips uniform 0.8 #> ${orig}/run_orig.log

#echo gdb ./learnDRMN 
#echo test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none $bk learnCV:0:123 ips uniform 0.8 #> ${bk}/run_bk.log
date
time ./learnDRMN_orig test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none $bk learnCV:0:123 ips uniform 0.8 &> ${orig}/run_orig.log
date
time ./learnDRMN test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none $bk learnCV:0:123 ips uniform 0.8 &> ${bk}/run_bk.log

orig=test_out_orig_chromatin_motif
bk=test_out_bookkeep_chromatin_motif
mkdir -p $orig $bk
date
time ./learnDRMN_orig test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/motif_chromatin_config.txt none $orig learnCV:0:123 ips uniform 0.8 &> ${orig}/run_orig.log
date
time ./learnDRMN test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/motif_chromatin_config.txt none $bk learnCV:0:123 ips uniform 0.8 &> ${bk}/run_bk.log


exit
echo gdb ./learnDRMN 
echo run test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none test_out learnCV:0:123 ips uniform 0.8
