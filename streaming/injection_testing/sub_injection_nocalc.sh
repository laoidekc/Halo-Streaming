#!/bin/bash
make
qsub -l nodes=1:ppn=2 -v num_nodes=1,ppn=2 injection_testing_nocalc.sh
qstat -u $USER
