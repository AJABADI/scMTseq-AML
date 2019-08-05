#!/bin/bash
read -p "Confirm that you want to sync all the sc_sPLS files from HPC to mac with Y or y " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rsync -avuz --progress --exclude='.*' --exclude="**/tmp" --exclude="**/*.cov" ajabadi@spartan.hpc.unimelb.edu.au:/data/cephfs/punim0613/AL/sc_sPLS/ ~/Documents/_Projects/sc_sPLS/
fi