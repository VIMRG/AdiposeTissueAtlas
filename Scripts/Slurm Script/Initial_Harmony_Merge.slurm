#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem=200G
#SBATCH --output=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/Initial_Merge_Harmony_7.11_%A.out
#SBATCH --error=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/Initial_Merge_Harmony_7.11_%A.err
#SBATCH --mail-user=samuel.bailin@vumc.org
#SBATCH --mail-type=FAIL

# Load Module
module load Anaconda2/2019.10

#  Activate Conda
eval "$(conda shell.bash hook)"
conda activate R_env

Rscript --no-save /data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/Initial_Merge.R
