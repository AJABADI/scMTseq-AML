#!/bin/bash
# CMD+F the name of the file and change it if needed

## ssh ajabadi@spartan.hpc.unimelb.edu.au
## cd /data/cephfs/punim0613/AL/sc_sPLS/sc_sPLS_Rproj/02_Methylation/

## sbatch 02_QC_Methylation.slurm
## squeue --job 7890452
## cat slurm-xxxx
## ls ../../../../data/scBSseq/met/parsed/tmp/
## cd ~/Documents/_Projects/sc_sPLS
## rsync -auvz --exclude-from='**/.*' ~/Documents/_Projects/sc_sPLS/sc_sPLS_Rproj/02_Methylation/02_Methylation/02_QC_Methylation/ ajabadi@spartan.hpc.unimelb.edu.au:/data/cephfs/punim0613/AL/sc_sPLS/sc_sPLS_Rproj/02_Methylation/02_Methylation/02_QC_Methylation/

# To give your job a name
#SBATCH --job-name=02_QC_Methylation

# For R need to run on 4 CPUs
#SBATCH --ntasks=4
#SBATCH --mem 10000

# set your walltime=hours:minutes:seconds
#SBATCH -t 02:30:00

# Specify your email address to be notified of progress.
#SBATCH --mail-user=aljaabadi@gmail.com
#SBATCH --mail-type=ALL

# Load the environment variables for R
module load pandoc/2.0.3
module load R/3.5.2-spartan_gcc-6.2.0

# The command to actually run the job
R --vanilla < 02_QC_Methylation.R
