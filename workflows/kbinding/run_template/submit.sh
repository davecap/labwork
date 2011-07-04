#!/bin/bash
#PBS -l nodes=${nodes}:ppn=${ppn},walltime=${walltime}${pbs_extra}
#PBS -N ${job_name}

# $PBS_O_WORKDIR
# $PBS_JOBID

cd $job_dir

module load gcc/4.4.0
module load python       
source ~/ENV/bin/activate

python ${pydr_path} -j $PBS_JOBID --pbs-nodefile $PBS_NODEFILE
