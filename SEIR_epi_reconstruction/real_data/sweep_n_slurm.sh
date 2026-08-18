#!/bin/bash
#SBATCH --job-name=sweep_n
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/real_data/sweep_n_%j.log
export LD_LIBRARY_PATH="/uufs/chpc.utah.edu/common/home/u1418987/.local/share/r-miniconda/lib:$LD_LIBRARY_PATH"
cd /uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction
Rscript -e '.libPaths(c("/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4", .libPaths())); source("real_data/sweep_n_bilstm.R")'
