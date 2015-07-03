#!/bin/bash
rm out.bin out.dat halo_out.bin halo_out.dat
make
mpiexec -n $1 ./halo_streaming $1
hexdump -v -e '10/8 "%f "' -e '"\n"' < out.bin > out.dat
hexdump -v -e '10/8 "%f "' -e '"\n"' < halo_out.bin > halo_out.dat
diff out.dat halo_out.dat
