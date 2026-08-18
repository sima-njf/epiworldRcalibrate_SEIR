#!/bin/bash
#SBATCH --job-name=current61_5method
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/real_data/current61_5method_%j.log

cd ~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction

export LD_LIBRARY_PATH="/uufs/chpc.utah.edu/common/home/u1418987/.local/share/r-miniconda/lib:$LD_LIBRARY_PATH"

Rscript -e '.libPaths(c("/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4", .libPaths())); source("real_data/calibrate_current61_5method.R")'
