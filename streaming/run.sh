#!/bin/bash
make
mpiexec -n $1 ./halo_streaming
