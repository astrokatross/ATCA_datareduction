#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT=/data/ATCA/ATCA_datareduction/
EPOCH=epoch5
BAND=L
TARGET=J215436

export PROJECT
export EPOCH
export BAND
export TARGET

cd $PROJECT/processing/

/data/bin/casa-6.2.0-124/bin/python3 ${PROJECT}bin/run_process.py 


