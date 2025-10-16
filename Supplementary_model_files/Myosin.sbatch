#!/bin/bash
#SBATCH --job-name=myosin_app
#SBATCH --output=myosin_app.out
#SBATCH --error=myosin_app.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --account=bidone
#SBATCH --partition=granite
#SBATCH --qos=granite
#SBATCH --array=1-8%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u1201504@umail.utah.edu

# Load MATLAB module
module load matlab

# Go to folder containing MATLAB functions
cd /uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V1_DifferentClutch

# Run MATLAB and ensure the folder is on the path
srun --exclusive -n 1 matlab -nodisplay -r "addpath(genpath(pwd)); MyosinApplication(${SLURM_ARRAY_TASK_ID}); exit"

wait
