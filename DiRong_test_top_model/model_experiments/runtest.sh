#!/bin/bash

num_process=(100 200 400 800 1200 1600)

rm run/atm/atm_demo/data/record_time
rm run/ocn/ocn_demo/data/record_time
source source_env.sh

for((ii=0;ii<6;ii++))
do
	for((k=0;k<${#num_process[@]};k++))
	do
		sed -i '28s/.*num_total_proc=.*/num_total_proc='${num_process[$k]}'/' config/common/case.conf
		sed -i '35s/.*num_total_proc=.*/num_total_proc='${num_process[$k]}'/' config/common/case.conf
		./configure
		./runcase
		grep -rin "DEFAULT_WGT_of.* exists" CCPL_dir/run/CCPL_logs/by_executables/atm_demo/atm_demo.CCPL.log.* > /dev/null
		if [ $? = 0 ]; then
		        echo "used atm weight file" >> run/atm/atm_demo/data/record_time
		fi
		grep -rin "DEFAULT_WGT_of.* exists"  CCPL_dir/run/CCPL_logs/by_executables/ocn_demo/ocn_demo.CCPL.log.* > /dev/null
		if [ $? = 0 ]; then
		        echo "used ocn weight file" >> run/atm/atm_demo/data/record_time
		fi
		echo "finish the $ii ${num_process[$k]} process"
	done
done 
