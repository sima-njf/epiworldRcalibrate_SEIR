#!/bin/sh
#SBATCH --job-name=part4_array
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/part4_scratch/part4_array/02-output-%A-%a.out
#SBATCH --array=1-100
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=part4_array
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  /uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/part4_scratch/part4_array/00-rscript.r
