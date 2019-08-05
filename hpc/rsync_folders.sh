#!/bin/bash
# 
src=/Users/alabadi/Documents/_Projects/sc_sPLS/sc_sPLS_Rproj/02_Methylation/01_Preprocessing_Methylation/03_call-contexts-and-sample-stats/
dest=ajabadi@spartan.hpc.unimelb.edu.au:/data/cephfs/punim0613/AL/sc_sPLS/sc_sPLS_Rproj/02_Methylation/01_Preprocessing_Methylation/03_call-contexts-and-sample-stats/
echo "Source folder with trailing slash to rsync: \n"
read src
echo "Destination folder with trailing slash to rsync: \n"
read dest
rsync -avz --progress --exclude='.*' $src $dest
