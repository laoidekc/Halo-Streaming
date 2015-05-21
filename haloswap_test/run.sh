#!/bin/bash
rm halo_out.bin halo_out.dat
make
mpiexec -n $1 ./halo_swapping
hexdump -v -e '10/8 "%f "' -e '"\n"' < halo_out.bin > halo_out.dat
