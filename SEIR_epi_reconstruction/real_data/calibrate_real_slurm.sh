#!/bin/bash
#SBATCH --job-name=real_5method
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction/real_data/calibrate_real_%j.log

cd ~/sima/epiworldRcalibrate_SEIR/SEIR_epi_reconstruction

# epiworldR.so needs a newer libstdc++ (GLIBCXX_3.4.29) than the system one;
# r-miniconda bundles one.
export LD_LIBRARY_PATH="/uufs/chpc.utah.edu/common/home/u1418987/.local/share/r-miniconda/lib:$LD_LIBRARY_PATH"
export SKIP_MEASLES=1

Rscript -e '.libPaths(c("/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4", .libPaths())); source("real_data/calibrate_real_5method.R")'
