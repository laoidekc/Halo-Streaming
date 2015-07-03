#!/bin/bash
make clean
make
qsub -pe mpi $1 halo.sge $1
