#/bin/bash

rm record_time_atm
rm record_time_ocn

atm_nml="atm_demo.nml"
ocn_nml="ocn_demo.nml"

num_process=(60 120 240 480)
atm_grid_size=(500000 1000000 2000000 4000000)
ocn_grid=("CUBE_grid_2.5.nc" "CUBE_grid_1.nc" "CUBE_grid_0.3.nc" "CUBE_grid_0.1.nc")

for((k=0;k<${#num_process[@]};k++))
do
	num_proc=${num_process[$k]}
	for((i=0;i<${#atm_grid_size[@]};i++))
	do
	    if [ $i = 0 -o $i = 1 ];then
                continue
            fi
	    sed -i 's/num_point.*/num_point='${atm_grid_size[$i]}'/' atm_demo.nml
	    for((j=0;j<${#ocn_grid[@]};j++));
	    do
		sed -i 's/grid_file_name.*/grid_file_name='${ocn_grid[$j]}'/' ocn_demo.nml
		sed -i 's/weight_file_index.*/weight_file_index='$i$j'/' atm_demo.nml
		echo ${atm_grid_size[$i]} ${ocn_grid[$j]} $num_proc
	
		node_count=0
		num_total_proc=`expr $num_proc + $num_proc`
		NUM_CORES_PER_NODE=24
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

		srun -p Big_Job -N $node_num -n $num_total_proc -o log/output${k}${i}${j}.log -e log/error${k}${i}${j}.log -J mct_experiment -l --multi-prog $conf_file
		echo ${atm_grid_size[$i]} ${ocn_grid[$j]} >> record_time_atm
		echo ${atm_grid_size[$i]} ${ocn_grid[$j]} >> record_time_ocn
	    done
	done
done
