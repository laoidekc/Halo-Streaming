#!/bin/bash
make clean
bolt -n $1 -N $2 -t 5:0:0 -o halo_streaming.bolt -j halo_streaming -A d68 ./halo_streaming $3
make archer
qsub halo_streaming.bolt
watch qstat -u $USER
