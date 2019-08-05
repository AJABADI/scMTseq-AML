#!/bin/bash
# have you rsynced the local changes?
# CMD+F the name of the file and change it if needed

## ssh ajabadi@spartan.hpc.unimelb.edu.au
## cd /data/cephfs/punim0613/AL/sc_sPLS/sc_sPLS_Rproj/02_Methylation/01_Preprocessing_Methylation/03_call-contexts-and-sample-stats
## sbatch 03_call-contexts-and-sample-stats.slurm
## squeue --job 7889276
## ls /data/cephfs/punim0613/AL/sc_sPLS/data/scBSseq/met/parsed/cov_context

# To give your job a name
#SBATCH --job-name=03_call-contexts-and-sample-stats

# For R need to run on 4 CPUs
#SBATCH --ntasks=2
##SBATCH --mem 10000

# set your walltime=hours:minutes:seconds
#SBATCH -t 02:30:00

# Specify your email address to be notified of progress.
#SBATCH --mail-user=aljaabadi@gmail.com
#SBATCH --mail-type=ALL

# Load the environment variables for R
module load pandoc/2.0.3
module load R/3.5.2-spartan_gcc-6.2.0

# The command to actually run the job
R --vanilla < 03_call-contexts-and-sample-stats_render.R
