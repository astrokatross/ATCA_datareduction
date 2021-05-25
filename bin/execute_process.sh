#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT=/data/var_analysis/ATCA/Code/
EPOCH=epoch1
BAND=L
PRIMARY_CALIBRATOR=1934_cal_l
SECONDARY_CALIBRATOR=j2327-4510
TARGET=J001513
TARGET_NAME=J001513

export PROJECT
export SRC
export EPOCH
export BAND
export PRIMARY_CALIBRATOR
export SECONDARY_CALIBRATOR
export TARGET
export TARGET_NAME

cd $PROJECT/processing/

/data/bin/casa-6.1.2-7-pipeline-2020.1.0.36/bin/python3 $PROJECT/bin/run_process.py 
