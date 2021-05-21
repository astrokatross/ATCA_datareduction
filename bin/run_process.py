#!/usr/bin/python
# This script calls process.py with all functions to analyse ATCA data and executes them in order, only changes needed to script should just be what stepyou need/commenting out whatever you don't need 
# TODO: Correct the sourcepar to be less hardcoded
# By K.Ross 19/5/21

# Importing relevant python packages
import process

data_dir = os.environ['PROJECT']
src = os.environ['SRC']
epoch = os.environ['EPOCH']
ATCA_band = os.environ["BAND"]
pri = os.environ['PRIMARY_CALIBRATOR']
sec = os.environ['SECONDARY_CALIBRATOR']
tar = os.environ['TARGET']
tar_nm = os.environ['TARGET_NAME']
sourcepar = [
    0.1, 12.6, -0.9
]

# Defining constant variables for all sources
src_dir = data_dir + src
process_dir = data_dir + "processing/"
img_dir = data_dir + src + "images/"
visname = data_dir + "{0}_{1}.ms".format(epoch, ATCA_band)
msname = data_dir + "{0}_{1}_{2}.ms".format(epoch, ATCA_band, tar_nm)
targetms = data_dir + "{0}_{1}_{2}_img.ms".format(epoch, ATCA_band, tar_nm)
tar_ms = data_dir + "{0}_{1}_{2}_selfcal.ms".format(epoch, ATCA_band, tar_nm)
ref = "CA04"
if ATCA_band == "L":
    n_spw = 8
    if_centre = 0
elif ATCA_band == "C":
    n_spw = 5
    if_centre = 0
elif ATCA_band == "X":
    n_spw = 4
    if_centre = 1


# Uncomment whichever step you don't need to do 
process.flag_ms(visname,epoch,ATCA_band,pri,sec,tar,tar_nm)
split_ms(visname, msname,  if_centre, epoch, ATCA_band, pri,sec,tar,tar_nm)
calibrate_ms(msname, epoch, ATCA_band, ref, pri,sec,tar,tar_nm)
applycal_ms(msname, epoch, ATCA_band, pri,sec,tar)
inspectpostcal_ms(msname, epoch, ATCA_band, pri,sec,tar)
flagcal_ms(msname, epoch, ATCA_band, pri,sec)
flagcaltar_ms(msname, epoch, ATCA_band, pri,sec,tar)
imgmfs_ms(msname, targetms, epoch, ATCA_band, tar,tar_nm)
img_ms(targetms, epoch, ATCA_band, n_spw,tar)
slefcal_ms(targetms, epoch, ATCA_band, n_spw,tar)
pbcor_ms(targetms, epoch, ATCA_band, n_spw,tar)
measureflux_ms(targetms, tar_ms, epoch, ATCA_band, sourcepar,n_spw,tar, tar_nm)
general_cleanup(process_dir, img_dir, src_dir)