#!/bin/bash
make
qsub -l nodes=1:ppn=2 -v num_nodes=1,ppn=2 extraction_testing.sh
watch qstat -u $USER
