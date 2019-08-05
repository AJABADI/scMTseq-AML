#!/bin/bash
# CMD+F the name of the file and change it if needed

## ssh ajabadi@spartan.hpc.unimelb.edu.au
## cd /data/cephfs/punim0613/AL/FULL_PATH_TO_template_name
## sbatch template_name.slurm
## squeue --job 7883870

# To give your job a name
#SBATCH --job-name=template_name

# For R need to run on 4 CPUs
#SBATCH --ntasks=2
##SBATCH --mem 8000

# set your walltime=hours:minutes:seconds
#SBATCH -t 00:30:00

# Specify your email address to be notified of progress.
#SBATCH --mail-user=aljaabadi@gmail.com
#SBATCH --mail-type=ALL

# Load the environment variables for R
module load pandoc/2.0.3
module load R/3.5.2-spartan_gcc-6.2.0

# The command to actually run the job
R --vanilla < template_name_render.R
