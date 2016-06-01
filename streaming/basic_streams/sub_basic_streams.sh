#!/bin/bash
make
qsub -l nodes=1:ppn=3 -v num_nodes=1,ppn=3 run_basic_streams.sh
watch qstat -u $USER
