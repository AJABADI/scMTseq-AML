#!/bin/bash
# bash /Users/alabadi/Documents/_Projects/sc_sPLS/pre_rsync_backup.sh
# bash /Users/alabadi/Documents/_Projects/sc_sPLS/_sc_sPLS_mac2hpc.sh
cd ~/Documents/_Projects/sc_sPLS
echo "Enter ssh password for Spartan at promot"
rsync -avz --exclude-from='exclude-mac2hpc.txt' --delete ~/Documents/_Projects/sc_sPLS/ ajabadi@spartan.hpc.unimelb.edu.au:/data/cephfs/punim0613/AL/sc_sPLS/
# no delete and keep the most recent in two
rsync -auvz --exclude-from='exclude-mac2hpc.txt' ~/Documents/_Projects/sc_sPLS/ ajabadi@spartan.hpc.unimelb.edu.au:/data/cephfs/punim0613/AL/sc_sPLS/
## -b get a backup of the ones existing in both adding ~ to names
## --progress
## go the directory and list files
ssh ajabadi@spartan.hpc.unimelb.edu.au