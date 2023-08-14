#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/Dialogue_4.20_%A.out
#SBATCH --error=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/Dialogue_4.20_%A.err
#SBATCH --mail-user=samuel.bailin@vumc.org
#SBATCH --mail-type=FAIL

# Load Module
module load Anaconda2/2019.10

#  Activate Conda
eval "$(conda shell.bash hook)"
conda activate R_env

Rscript --no-save /data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/CodeRevised/Dialogue.R
#Rscript --no-save /data/p_koethe_lab/Atlas_AT/HATIM_Analysis/1.7.23/CodeRevised/Dialogue_HIVneg.R