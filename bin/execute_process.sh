#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT=/data/var_analysis/ATCA/Code/
EPOCH=epoch3
BAND=X
TARGET=J024838

export PROJECT
export EPOCH
export BAND
export TARGET

cd $PROJECT/processing/

/data/bin/casa-6.1.2-7-pipeline-2020.1.0.36/bin/python3 $PROJECT/bin/run_process.py 
