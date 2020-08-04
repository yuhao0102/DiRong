#/bin/bash

source ~/module_load.sh
rm record_time_atm
rm record_time_ocn

atm_nml="atm_demo.nml"
ocn_nml="ocn_demo.nml"

num_process=(100 200 400 800 1200 1600)

for((ii=0;ii<6;ii++))
do
	for((k=0;k<${#num_process[@]};k++))
	do
		num_proc=${num_process[$k]}
		sed -i 's/num_point.*/num_point=4000000/' atm_demo.nml
		sed -i 's/grid_file_name.*/grid_file_name=CUBE_grid_0.3.nc/' ocn_demo.nml
		sed -i 's/weight_file_index.*/weight_file_index=32/' atm_demo.nml
		echo 4000000 CUBE_grid_0.3.nc $num_proc
	
		node_count=0
		num_total_proc=`expr $num_proc + $num_proc`
		NUM_CORES_PER_NODE=64
		let node_count=num_total_proc%NUM_CORES_PER_NODE
		if (( $node_count == 0 )); then
		    let node_num=num_total_proc/NUM_CORES_PER_NODE
		else
		    let node_num=num_total_proc/NUM_CORES_PER_NODE+1
		fi
		conf_file="job.conf"
		rm $conf_file
		cat >> $conf_file << EOF    
0-`expr $num_proc - 1` `pwd`/atm_demo
$num_proc-`expr $num_proc + $num_proc - 1` `pwd`/ocn_demo
EOF

		srun -p amd_256 -N $node_num -n $num_total_proc -o log/output${k}.log -e log/error${k}.log -J mct_experiment -l --multi-prog $conf_file
		echo 4000000 CUBE_grid_0.3.nc >> record_time_atm
		echo 4000000 CUBE_grid_0.3.nc >> record_time_ocn
	done
done
