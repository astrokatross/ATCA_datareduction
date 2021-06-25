#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT=/data/ATCA/Code/
EPOCH=epoch2
BAND=L
TARGET=J042502

export PROJECT
export EPOCH
export BAND
export TARGET

cd $PROJECT/processing/

python3 $PROJECT/bin/run_process.py 
