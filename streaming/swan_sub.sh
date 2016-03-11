#PBS -l walltime=00:01:00
#PBS -j oe
cd $PBS_O_WORKDIR

MPIPROG=`basename $PBS_JOBNAME .pbs`
num_procs=$(($num_nodes*$ppn))
echo '--------------------------------------------------------------------------------'

echo 'Running MPI program' $MPIPROG 'on ' $num_procs ' processes'

echo 'Started at' `date`
echo '--------------------------------------------------------------------------------'

(time aprun -n $num_procs ./halo_streaming $num_procs) 2>&1

echo '--------------------------------------------------------------------------------'
echo 'Finished at' `date`
