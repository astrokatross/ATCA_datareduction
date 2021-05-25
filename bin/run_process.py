#!/data/bin/casa-6.1.2-7-pipeline-2020.1.0.36/bin/python3
# This script calls process.py with all functions to analyse ATCA data and executes them in order, only changes needed to script should just be what stepyou need/commenting out whatever you don't need
# TODO: Change all the plotms commands in process.py to be separate commands here using casaplotms.plotms

# By K.Ross 19/5/21

# Importing relevant python packages
import os
import process

data_dir = str(os.environ["PROJECT"])
epoch = str(os.environ["EPOCH"])
ATCA_band = str(os.environ["BAND"])
pri = str(os.environ["PRIMARY_CALIBRATOR"])
sec = str(os.environ["SECONDARY_CALIBRATOR"])
tar = str(os.environ["TARGET"])
tar_nm = str(os.environ["TARGET_NAME"])
print(tar_nm)
# Setting sourcepar dictionary to measrue flux
sourcepar_dict = {
    "J001513": [0.1, 12.6, -0.9],
    "J015445": [0.15, 14.1, -2.7],
    "J020507": [0.3, 2.6, 0.7],
    "J021246": [0.15, 12.9, -5.2],
    "J022744": [0.25, 6, 0.3],
    "J024838": [0.1, 1.5, 0.5],
    "J032213": [0.35, 2.6, -3.4],
    "J032836": [0.35, 4, -2.5],
    "J033023": [0.15, 0, 0],
    "J042502": [0.1, 4.8, 0],
    "J044033": [0.2, 8.2, -0.3],
    "J044737": [0.25, 3.4, -2.2],
    "J052824": [0.2, 8, -3.9],
    "J223933": [0.2, 10, -3.1],
    "J224408": [0.25, 8.3, -5.2],
}

# Defining constant variables for all sources
src_dir = f"{data_dir}{tar_nm}"
process_dir = f"{data_dir}processing/"
img_dir = f"{data_dir}{tar_nm}/images/"
visname = f"{data_dir}data/{epoch}_{ATCA_band}.ms"
msname = f"{data_dir}data/{epoch}_{ATCA_band}_{tar_nm}.ms"
targetms = f"{data_dir}data/{epoch}_{ATCA_band}_{tar_nm}_img.ms"
tar_ms = f"{data_dir}data/{epoch}_{ATCA_band}_{tar_nm}_selfcal.ms"
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
process.flag_ms(img_dir, visname, epoch, ATCA_band, pri, sec, tar, tar_nm)
# process.split_ms(visname, msname, if_centre, epoch, ATCA_band, pri, sec, tar, tar_nm)
# process.calibrate_ms(msname, epoch, ATCA_band, ref, pri, sec, tar, tar_nm)
# process.applycal_ms(msname, epoch, ATCA_band, pri, sec, tar)
# process.inspectpostcal_ms(msname, epoch, ATCA_band, pri, sec, tar)
# process.flagcal_ms(msname, epoch, ATCA_band, pri, sec)
# process.flagcaltar_ms(msname, epoch, ATCA_band, pri, sec, tar, tar_nm)
# process.imgmfs_ms(msname, targetms, epoch, ATCA_band, n_spw, tar, tar_nm)
# process.img_ms(targetms, epoch, ATCA_band, n_spw, tar, tar_nm)
# process.slefcal_ms(targetms, epoch, ATCA_band, n_spw, tar, tar_nm)
# process.pbcor_ms(targetms, epoch, ATCA_band, n_spw, tar, tar_nm)
# process.measureflux_ms(targetms, tar_ms, epoch, ATCA_band, sourcepar, n_spw, tar,
#                tar_nm)
# process.general_cleanup(process_dir, img_dir, src_dir)

