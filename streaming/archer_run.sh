#!/bin/bash
rm out.bin out.dat halo_out.bin halo_out.dat
make
qsub -q short -l select=$1 halo_streaming.pbs