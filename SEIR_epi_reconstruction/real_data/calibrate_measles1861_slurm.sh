#!/bin/bash
#SBATCH --job-name=measles1861_5method
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/real_data/measles1861_5method_%j.log

cd ~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction
export LD_LIBRARY_PATH="/uufs/chpc.utah.edu/common/home/u1418987/.local/share/r-miniconda/lib:$LD_LIBRARY_PATH"
export SKIP_WAVE1=1
export SKIP_MEASLES=1

Rscript -e '.libPaths(c("/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4", .libPaths())); source("real_data/calibrate_measles1861_5method.R")'
