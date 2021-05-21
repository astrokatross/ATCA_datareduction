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

# Setting sourcepar dictionary to measrue flux
sourcepar_dict = {
    "J001513" : [0.1,12.6,-0.9],
    "J015445" : [0.15,14.1,-2.7],
    "J020507" : [0.3,2.6,0.7],
    "J021246" : [0.15,12.9,-5.2],
    "J022744" : [0.25,6,0.3],
    "J024838" : [0.1,1.5,0.5],
    "J032213" : [0.35,2.6,-3.4],
    "J032836" : [0.35,4,-2.5],
    "J033023" : [0.15,0,0],
    "J042502" : [0.1,4.8,0],
    "J044033" : [0.2,8.2,-0.3],
    "J044737" : [0.25,3.4,-2.2],
    "J052824" : [0.2,8,-3.9],
    "J223933" : [0.2,10,-3.1],
    "J224408" : [0.25,8.3,-5.2]
}

# Defining constant variables for all sources
src_dir = data_dir + src
process_dir = data_dir + "processing/"
img_dir = data_dir + src + "images/"
visname = data_dir + "{0}_{1}.ms".format(epoch, ATCA_band)
msname = data_dir + "{0}_{1}_{2}.ms".format(epoch, ATCA_band, tar_nm)
targetms = data_dir + "{0}_{1}_{2}_img.ms".format(epoch, ATCA_band, tar_nm)
tar_ms = data_dir + "{0}_{1}_{2}_selfcal.ms".format(epoch, ATCA_band, tar_nm)
ref = "CA04"
sourcepar = sourcepar_dict[tar_nm]
if ATCA_band == "L":
    n_spw = 8
    if_centre = 0
elif ATCA_band == "C":
    n_spw = 5
    if_centre = 0
elif ATCA_band == "X":
    n_spw = 4
    if_centre = 1

print("Here we go! Time to analyse some ATCA data!")
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