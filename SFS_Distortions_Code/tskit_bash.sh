#!/bin/bash

#SBATCH --output=tskitOutputGrowth2000_compPlots_scaled.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --time=40:00:00
#SBATCH --partition=standard
#SBATCH --account=masel

conda init
conda activate dadi_env

#dadi-cli InferDM --fs dadi_SFS.fs --output-prefix dadi_output --model bottlegrowth_1d --nomisid --lbounds 0.01 0.01 0.001 --ubounds 100 100 10 --optimizations 20
python tskitAnalysis_growth.py
