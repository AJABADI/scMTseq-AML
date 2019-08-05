#!/bin/bash
## UPDATES TO compress_and_backup.sh in: /Users/alabadi/Documents/_Projects/Cheatsheet/GitHub/01-computer/bash/rsync/compress_and_backup.sh
###################################################
## backs up all the specified filetypes from FOLDER in DEST as compressed files and logs the changes
## and user's message in DEST/backup_log.log - removed the oldest ones as specified in 'awk'.
## cd to FOLDER's parent dir first
###################################################
##INPUT FOLDER parent dir
cd ~/Documents/_Projects/sc_sPLS
## name of source folder, destination, and the file
##INPUT basename
FOLDER=sc_sPLS_Rproj
##INPUT full path
DEST=__local
backup_log=$DEST/backup.log
## extensions being backedup: log, R, Rmd
## customise the number of backups to keep in last line

############### for example purpose only ^^^^.
## make sure a log file exists in $DEST
if [ ! -f $backup_log ]; then
    printf "initialising the log file on... $(date +%A-%d-%B-%Y-%H-%M) \n" >> $backup_log
fi

##INPUT backup file name
BACKUP="${DEST}/${FOLDER}_$(date +%Y%m%d).tar.gz"
printf "\nfull backup name: $BACKUP \n \n"
## user input for the run
read -p "Enter a message for the backup file of the designated files: "  msg
printf " $(date +%A-%d-%B-%Y-%H-%M) $msg  backedupfie name: $BACKUP \n\n" | cat - $backup_log > temp && mv temp $backup_log

##INPUT all file extensions to compress
find ./$FOLDER -name "*.log" -o -name "*.R" -o -name "*.Rmd" | tar -cf $BACKUP -T -

## INPUT keep the last ? backups only
rm $(ls -td $DEST/*.tar.gz | awk 'NR>50')
## copy as pathname and go bash <full-path.sh>