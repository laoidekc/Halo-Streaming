#!/bin/bash
make
for i in 2 22 44
do
	qsub -N "procs_$i"-l nodes=1:ppn=$i -v num_nodes=1,ppn=$i run_injection_testing.sh
done

qsub -hold_jid procs_2,procs_22,procs_44 ./results_script
