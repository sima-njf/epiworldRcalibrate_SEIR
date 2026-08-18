#!/bin/bash
#SBATCH --job-name=seir_utah_calib
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/real_data/seir_utah_calib_%j.log

cd ~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction

export LD_LIBRARY_PATH="/uufs/chpc.utah.edu/common/home/u1418987/.local/share/r-miniconda/lib:$LD_LIBRARY_PATH"

Rscript -e '.libPaths(c("/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4", .libPaths())); source("real_data/seir_utah_calibration.R")'
