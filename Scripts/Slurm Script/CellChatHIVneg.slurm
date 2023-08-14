#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=18:00:00
#SBATCH --mem=240G
#SBATCH --output=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/CellChat_4.14_%A.out
#SBATCH --error=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/CellChat_4.14_%A.err
#SBATCH --mail-user=samuel.bailin@vumc.org
#SBATCH --mail-type=FAIL
##SBATCH --constraint="sandybridge|haswell|skylake"

# Load Module
module load Anaconda2/2019.10

#  Activate Conda
eval "$(conda shell.bash hook)"
conda activate R_env

Rscript --no-save /data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/CodeRevised/CellChat_HIVneg.R
