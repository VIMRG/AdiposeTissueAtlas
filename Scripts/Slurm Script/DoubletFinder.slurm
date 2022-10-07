#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:30:00
#SBATCH --mem=55G
#SBATCH --output=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/DoubletFinder_%A_a%.out
#SBATCH --error=/data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Logs/DoubletFinder_%A_a%.err
#SBATCH --mail-user=samuel.bailin@vumc.org
#SBATCH --mail-type=FAIL
#SBATCH --array=1-23%23

# Print slurm Array ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

#  Activate Conda
eval "$(conda shell.bash hook)"
conda activate R_env

Rscript --no-save /data/p_koethe_lab/Atlas_AT/HATIM_Analysis/6.27/Code/DoubletFinder.R