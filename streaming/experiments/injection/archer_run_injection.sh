#!/bin/bash
make clean
for i in 1 10 100 1000 10000 100000 1000000
do
	bolt -n $1 -N $2 -t 0:2:0 -o injection.bolt -j injection -A d35-clk ./injection_testing 1000000 100 $i $1
	make
	qsub -a 1745 injection.bolt
done
watch qstat -u $USER
