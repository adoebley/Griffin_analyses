#!/bin/bash

#SBATCH --job-name=run_readCounter_1mb.sh
#SBATCH --partition=campus-new
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=1:0:0
#SBATCH -o slurm_logs//run_readCounter_1mb.sh.o.%A.%a

export PATH=/app/bin:$PATH
source /app/Lmod/lmod/lmod/init/bash
export LMOD_PACKAGE_PATH=/app/Lmod
module use /app/easybuild/modules/all
module use /app/Modules/modulefiles
ml R/3.6.0-foss-2016b-fh1

eval `sed -n "${SLURM_ARRAY_TASK_ID}{p;q;}" /fh/fast/ha_g/user/adoebley/ichorCNA/createPoN/1mb_delfi_HD/scripts/run_readCounter_1mb.sh`
