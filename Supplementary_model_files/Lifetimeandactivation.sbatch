#!/bin/bash
#SBATCH --job-name=Lifetimeandactivation       # Job name matches MATLAB file
#SBATCH --output=Lifetimeandactivation.out    # Single output file
#SBATCH --error=Lifetimeandactivation.err     # Single error file
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                             # One MATLAB instance per array task
#SBATCH --cpus-per-task=16                     # 16 cores available to MATLAB
#SBATCH --mem=128G
#SBATCH --account=bidone
#SBATCH --partition=granite
#SBATCH --qos=granite
#SBATCH --array=1-54%16                        # 54 total simulations, max 16 running simultaneously
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u1201504@umail.utah.edu

# Load MATLAB module
module load matlab

# Go to folder containing MATLAB scripts
cd /uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V1_DifferentClutch

# Run MATLAB in batch mode, add folder to MATLAB path, and pass the array task ID
srun --exclusive -n 1 matlab -nodisplay -nosplash -r "addpath(genpath(pwd)); Lifetimeandactivation(${SLURM_ARRAY_TASK_ID}); exit"

wait
