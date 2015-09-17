#!/bin/bash
make clean
make
mpiexec -n $1 ./halo_streaming $1
hexdump -v -e '10/8 "%f "' -e '"\n"' < streaming_out.bin > streaming_out.dat
hexdump -v -e '10/8 "%f "' -e '"\n"' < exchange_out.bin > exchange_out.dat
diff streaming_out.dat exchange_out.dat
