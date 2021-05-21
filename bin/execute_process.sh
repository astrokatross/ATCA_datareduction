#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT="/data/var_analysis/ATCA/Code/"
SRC=""
EPOCH=""
BAND=""
PRIMARY_CALIBRATOR=""
SECONDARY_CALIBRATOR=""
TARGET=""
TARGET_NAME=""

export PROJECT,SRC,EPOCH,BAND,PRIMARY_CALIBRATOR,SECONDARY_CALIBRATOR,TARGET,TARGET_NAME

cd $PROJECT/processing/

casapy -c $PROJECT/bin/run_process.py 

# sourcepar = [
#     0.1, 12.6, -0.9
# ] 