#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT=/data/ATCA/ATCA_datareduction/
EPOCH=epoch1
BAND=L
TARGET=J015445

export PROJECT
export EPOCH
export BAND
export TARGET

cd $PROJECT/processing/

if [[ ! -e "${PROJECT}${TARGET}" ]]
then 
    echo "No folder for source found, making them"
    mkdir ${PROJECT}${TARGET}
    mkdir ${PROJECT}${TARGET}/casa_files
    mkdir ${PROJECT}${TARGET}/cal_tables
    mkdir ${PROJECT}${TARGET}/images 
fi

/data/bin/casa-6.2.0-124/bin/python3 ${PROJECT}bin/run_process.py 


