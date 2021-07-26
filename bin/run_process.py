#!/usr/bin/python3
# This script calls process.py with all functions to analyse ATCA data and executes them in order, only changes needed to script should just be what stepyou need/commenting out whatever you don't need
# By K.Ross 19/5/21

# Importing relevant python packages
import process
import os


# data_dir = str(os.environ["PROJECT"])
# epoch = str(os.environ["EPOCH"])
# ATCA_band = str(os.environ["BAND"])
# tar = str(os.environ["TARGET"])
data_dir = "/data/ATCA/ATCA_datareduction/"
epoch = "epoch5"
ATCA_band = "L"
tar = "J215436"

print(f"Target: {tar}\nEpoch: {epoch}\nATCA band: {ATCA_band}")

# Setting sourcepar dictionary to measrue flux
source_dict = {
    "J001513": ["2327-459", (0.1, 12.7, -0.5)],
    "J015445": ["0237-233", (0.15, 14.1, -2.7)],
    "J020507": ["0238-084", (0.3, 2.6, 0.7)],
    "J021246": ["0237-233", (0.15, 12.9, -5.2)],
    "J022744": ["0238-084", (0.25, 6, 0.3)],
    "J024838": ["0237-233", (0.1, 1.5, 0.5)],
    "J032213": ["0355-483", (0.35, 2.6, -3.4)],
    "J032836": ["0310-150", (0.35, 4, -2.5)],
    "J033023": ["0310-150", (0.15, 0, 0)],
    "J042502": ["0445-221", (0.1, 4.6, 1.9)],
    "J044033": ["0355-483", (0.2, 8.2, -0.3)],
    "J044737": ["0445-221", (0.25, 3.4, -2.2)],
    "J052824": ["0528-250", (0.2, 8, -3.9)],
    "J223933": ["2327-459", (0.2, 10, -3.1)],
    "J224408": ["2240-260", (0.25, 8.3, -5.2)],
    "J215436": ["2211-388", (0.3, 0, 0)],
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
export_pngs = True
sec = source_dict[tar][0]
sourcepar = source_dict[tar][1]

if ATCA_band == "L":
    n_spw = 8
    pri = "1934_cal_l"
    visname = f"{data_dir}data/{epoch}_{ATCA_band}_cut.ms"
elif ATCA_band == "C":
    n_spw = 5
    pri = "1934_cal_cx"
elif ATCA_band == "X":
    n_spw = 4
    pri = "1934_cal_cx"

print("Here we go! Time to analyse some ATCA data!")
# Uncomment whichever step you don't need to do
# Initial flagging
# process.flag_ms(visname)

# # Split to make its own ms
# process.split_ms(
#     src_dir,
#     img_dir,
#     visname,
#     msname,
#     epoch,
#     ATCA_band,
#     pri,
#     sec,
#     tar,
# )

# # Calibrate, and apply cal ms using primary and secondary
# process.calibrate_ms(src_dir, msname, epoch, ATCA_band, ref, pri, sec, tar)
# process.applycal_ms(src_dir, msname, epoch, ATCA_band, pri, sec, tar)

# # Post cal inspection and flagging
# process.flagcal_ms(img_dir, msname, epoch, ATCA_band, pri, sec)
# process.flagcaltar_ms(src_dir, img_dir, msname, epoch, ATCA_band, pri, sec, tar)

# Imaging of target
process.imgmfs_ms(src_dir, msname, targetms, epoch, ATCA_band, n_spw, tar)
# process.img_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar)
# process.slefcal_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar)

# # Post image analysis: pbcor, measure flux
# process.pbcor_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar)
# process.measureflux_ms(
#     src_dir, targetms, tar_ms, epoch, ATCA_band, sourcepar, n_spw, tar
# )

# # Post analysis
# if export_pngs is True:
#     process.inspection_plots(
#         src_dir, img_dir, visname, msname, epoch, ATCA_band, pri, sec, tar
#     )
# process.export_fitspng(src_dir, n_spw, epoch, ATCA_band, tar)
