#!/usr/bin/python
# This script calls process.py with all functions to analyse ATCA data and executes them in order, only changes needed to script should just be what stepyou need/commenting out whatever you don't need
# TODO: Compelte dictionaries so it has pri cal,sec cal for each epoch and maybe target name? That would be nice
# By K.Ross 19/5/21

# Importing relevant python packages
import os
import process

data_dir = str(os.environ["PROJECT"])
epoch = str(os.environ["EPOCH"])
ATCA_band = str(os.environ["BAND"])
tar = str(os.environ["TARGET"])
print(f"Target: {tar}\nEpoch: {epoch}\nATCA band: {ATCA_band}")

# Setting sourcepar dictionary to measrue flux
# Order of entries: secondary cal, (e1L, e1CX/e2L, e3L/e3CX/e4L/e4CX)
source_dict = {
    "J001513": ["j2327-4510", "j2327-4510", "2327-459", (0.1, 12.6, -0.9)],
    "J015445": ["1bac0237-233", "back0237", "0237-233", (0.15, 14.1, -2.7)],
    "J020507": ["1hed0238-084", "head0238", "0238-084", (0.3, 2.6, 0.7)],
    "J021246": ["1bac0237-233", "back0237", "0237-233", (0.15, 12.9, -5.2)],
    "J022744": ["1hed0238-084", "head0238", "0238-084", (0.25, 6, 0.3)],
    "J024838": ["1bac0237-233", "back0237", "0237-233", (0.1, 1.5, 0.5)],
    "J032213": ["1bc0355-483", "1breast0355", "0355-483", (0.35, 2.6, -3.4)],
    "J032836": ["1eye0310-150", "eye0310", "0310-150", (0.35, 4, -2.5)],
    "J033023": ["1eye0310-150", "eye0310", "0310-150", (0.15, 0, 0)],
    "J042502": ["1n0445-221", "Xneck0445", "0445-221", (0.1, 4.8, 0)],
    "J044033": ["1bc0355-483", "1breast0355", "0355-483", (0.2, 8.2, -0.3)],
    "J044737": ["1n0445-221", "Xneck0445", "0445-221", (0.25, 3.4, -2.2)],
    "J052824": ["2bc0528-250", "2breast0528", "0528-250", (0.2, 8, -3.9)],
    "J223933": ["j2327-4510", "j2327-4510", "2327-459", (0.2, 10, -3.1)],
    "J224408": ["2243-123", "2243-123", "2240-260", (0.25, 8.3, -5.2)],
}


# Defining constant variables for all sources
src_dir = f"{data_dir}{tar}"
process_dir = f"{data_dir}processing/"
img_dir = f"{data_dir}{tar}/images/"
visname = f"{data_dir}data/{epoch}_{ATCA_band}.ms"
msname = f"{data_dir}data/{epoch}_{ATCA_band}_{tar}.ms"
targetms = f"{data_dir}data/{epoch}_{ATCA_band}_{tar}_img.ms"
tar_ms = f"{data_dir}data/{epoch}_{ATCA_band}_{tar}_selfcal.ms"
ref = "CA04"
sourcepar = source_dict[tar][3]

if epoch == "epoch1":
    if ATCA_band == "L":
        sec = source_dict[tar][0]
    elif ATCA_band == "C" or ATCA_band == "X":
        sec = source_dict[tar][1]
elif epoch == "epoch2":
    sec = source_dict[tar][1]
elif epoch == "epoch3" or epoch == "epoch4":
    sec = source_dict[tar][2]


if ATCA_band == "L":
    n_spw = 8
    pri = "1934_cal_l"
elif ATCA_band == "C":
    n_spw = 5
    pri = "1934_cal_CX"
elif ATCA_band == "X":
    n_spw = 4
    pri = "1934_cal_CX"

print("Here we go! Time to analyse some ATCA data!")
# Uncomment whichever step you don't need to do
# Initial flagging and creating targetms
# process.flag_ms(img_dir, visname, epoch, ATCA_band, pri, sec, tar, tar)
process.split_ms(
    src_dir,
    img_dir,
    visname,
    msname,
    epoch,
    ATCA_band,
    pri,
    sec,
    tar,
    tar,
)

# Calibrate, and apply cal ms using primary and secondary
process.calibrate_ms(src_dir, msname, epoch, ATCA_band, ref, pri, sec, tar, tar)
process.applycal_ms(src_dir, msname, epoch, ATCA_band, pri, sec, tar)

# Post cal inspection and flagging
process.inspectpostcal_ms(img_dir, msname, epoch, ATCA_band, pri, sec, tar)
process.flagcal_ms(img_dir, msname, epoch, ATCA_band, pri, sec)
process.flagcaltar_ms(src_dir, img_dir, msname, epoch, ATCA_band, pri, sec, tar, tar)

# Imaging of target
process.imgmfs_ms(src_dir, msname, targetms, epoch, ATCA_band, n_spw, tar, tar)
process.img_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar, tar)
process.slefcal_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar, tar)

# Post image analysis: pbcor, measure flux
process.pbcor_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar, tar)
process.measureflux_ms(
    src_dir, targetms, tar_ms, epoch, ATCA_band, sourcepar, n_spw, tar, tar
)
