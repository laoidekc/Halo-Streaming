#!/bin/bash
make clean
make swan
qsub -l nodes=$1:ppn=$2 -v num_nodes=$1,ppn=$2 swan_sub.sh
watch qstat -u $USER
