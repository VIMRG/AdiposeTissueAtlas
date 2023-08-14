#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=100G
#SBATCH --output=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/Logs/DGE_Integrated_4.17_%A_%a.out
#SBATCH --error=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/Logs/DGE_Integrated_4.17_%A_%a.err
#SBATCH --mail-user=samuel.bailin@vumc.org
#SBATCH --mail-type=FAIL
##SBATCH --constraint="sandybridge|haswell|skylake"
#SBATCH --array=1-50%7

# Load Module
module load Anaconda2/2019.10

#  Activate Conda
eval "$(conda shell.bash hook)"
conda activate R_env

Rscript --no-save /data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/Code/Supplemental/ClusterDGE.R
