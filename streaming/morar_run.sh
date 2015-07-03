#!/bin/bash
rm out.bin out.dat halo_out.bin halo_out.dat
make
qsub -pe mpi $1 halo.sge $1
