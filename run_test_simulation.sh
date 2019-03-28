#!/bin/bash
# test run for simulation
set -u 

cluster=/mnt/ws/sysbio/roygroup/shared/projects/drmn/ancillary_code/learnMoE/learnMoE_tab/learnMoE

fold=0
k=3 # 1000 genes cannot support 5+ clusters, so keep it small if you are using the small input
out=test_simulation_k${k}
inconfig=test_data/histonly_config.txt

mkdir -p $out

# run it once on all data
./learnDRMN test_data/celltype_order.txt test_data/ogids_expr_1000.txt null ${k} test_data/lineage_tree.txt $inconfig none $out learn ips uniform 0.8

# assemble data for subsequent application to predicted expr
newconfig=${out}/sim_config.txt
printf "" > $newconfig
while read cell clust expr feat; do
	# locate expression
	cellexp=${out}/drmn/${cell}_sample_maxclust_exprtab.txt
	if [[ ! -e $cellexp ]]; then echo "ERROR: Cannot find ${cellexp}"; exit 1; fi
	
	# let's just use the old cluster assigmnent instead of making a new one

	# add line to new config file
	printf "%s\t%s\t%s\t%s\n" $cell $clust $cellexp $feat >> $newconfig
done < $inconfig

cat $newconfig

# Now run on the simulation in CV
simout=${out}/simulation
mkdir -p $simout
./learnDRMN test_data/celltype_order.txt test_data/ogids_expr_1000.txt null ${k} test_data/lineage_tree.txt $newconfig none $simout learnCV ips uniform 0.8 

