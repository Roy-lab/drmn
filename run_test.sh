#!/bin/bash
# run sanity check with histone mark features only, 

if [[ $# != 1 ]]; then
	echo "usage: run_test.sh outdir"
	exit
fi
out=$1
mkdir -p $out

echo gdb ./learnDRMN 
echo run test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none $out learnCV:0:123 ips uniform 0.8



exit
########################################3

echo gdb ./learnDRMN 
echo run test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/mef_tree.txt test_data/histonly_config.txt none test_short learnCV:0:123 ips uniform 0.8



exit

# compare individually run folds to togetherly run (fold=-1)
for fold in 0 1 2
do
	./learnDRMN test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none test_out_parallel learnCV:${fold}:123 ips uniform 0.8 &> run_test_fold${fold}.log &
done
./learnDRMN test_data/celltype_order.txt test_data/ogids_expr_100.txt test_data/all_regulator_ids.txt 3 test_data/lineage_tree.txt test_data/histonly_config.txt none test_out_serial learnCV:-1:123 ips uniform 0.8 &> run_test_allfolds.log
echo "done"

# compare test genes
for fold in 0 1 2
do
	for foldSer in 0 1 2
	do
		printf "p${fold} to s${foldSer}, pips test genes shared "
		comm <(cut -f 1 test_out_parallel/fold${fold}/drmn/pips_pred_test.txt | sort) <(cut -f 1 test_out_serial/fold${foldSer}/drmn/pips_pred_test.txt | sort) -1 -2 | wc -l
		printf "p${fold} to s${foldSer}, pips test genes not shared "
		comm <(cut -f 1 test_out_parallel/fold${fold}/drmn/pips_pred_test.txt | sort) <(cut -f 1 test_out_serial/fold${foldSer}/drmn/pips_pred_test.txt | sort) -3 | wc -l
	done
done
