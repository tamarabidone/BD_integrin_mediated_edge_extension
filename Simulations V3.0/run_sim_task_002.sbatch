#!/bin/bash
#SBATCH --job-name=sim_task_002
#SBATCH --output=sim_task_002_%A_%a.out  # Unique output per task
#SBATCH --error=sim_task_002_%A_%a.err   # Unique error per task
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6                      # 8 independent tasks
#SBATCH --cpus-per-task=1                # 1 core per task
#SBATCH --mem=128G
#SBATCH --account=bidone
#SBATCH --partition=granite
#SBATCH --qos=granite
#SBATCH --array=1-120%16               # Run 100 cases, max 8 at a time
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u1201504@umail.utah.edu

module load matlab
cd /uufs/chpc.utah.edu/common/home/bidone-group3/Remi/Simulations_V3.0

# Run MATLAB, passing the array task ID as an argument
srun --exclusive -n 1 matlab -nodisplay -r "SimulationTaskList_002($SLURM_ARRAY_TASK_ID); exit" &
wait  # Ensure all tasks complete before exiting
