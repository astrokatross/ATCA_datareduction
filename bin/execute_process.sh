#!/bin/bash
# This is a script to execute the python/casa script that will do the processing
# By K.Ross 

PROJECT=/data/var_analysis/ATCA/Code/
EPOCH=epoch1
BAND=C
PRIMARY_CALIBRATOR=1934_cal_CX
SECONDARY_CALIBRATOR=j2327-4510
TARGET=001513
TARGET_NAME=J001513
LOG_FILE=$PROJECT/casa_logs/${EPOCH}_${BAND}_${TARGET_NAME}.log

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

# casaplotms(vis=visname,field=pri,xaxis="channel",yaxis="amp",correlation="xy,yx",ydatacolumn="data",plotfile=f"{img_dir}{epoch}_{ATCA_band}_{pri}_ampvschan_pre_RFIflag_.png",showgui=False,overwrite=True)