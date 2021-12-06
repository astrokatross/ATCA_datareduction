#!/usr/bin/python3
# This script calls process.py with all functions to analyse ATCA data and executes them in order, only changes needed to script should just be what stepyou need/commenting out whatever you don't need
# By K.Ross 19/5/21

# TODO: introduce epoch processing 

# Importing relevant python packages
import process
import os


data_dir = str(os.environ["PROJECT"])
run_epoch = str(os.environ["EPOCH"])
ATCA_band = str(os.environ["BAND"])
tar = str(os.environ["TARGET"])
calibrator = str(os.environ["CALIBRATOR"])

print(f"Target: {tar}\nATCA band: {ATCA_band}")

# Setting sourcepar dictionary to measrue flux
source_dict = {
    "J001513": ["2327-459", (0.1, 12, -0.22), (0.5, 0, 0)],
    "J015445": ["0237-233", (0.15, 0, 0), (0.5, 0, 0)],
    "J020507": ["0159-117", (0.3, 2.6, 0.7), (0.5, 0, 0)],
    "J021246": ["0237-233", (0.15, 12.9, -5.2), (0.5, 0, 0)],
    "J022744": ["0238-084", (0.25, 6, 0.3), (0.5, 0, 0)],
    "J024838": ["0237-233", (0.1, 1.5, 0.5), (0.5, 0, 0)],
    "J032213": ["0355-483", (0.35, 2.6, -3.4), (0.5, 0, 0)],
    "J032836": ["0310-150", (0.35, 4, -2.5), (0.5, 0, 0)],
    "J033023": ["0336-019", (0.15, 1.1, 0.86), (0.5, 0, 0)],
    "J042502": ["0445-221", (0.1, 4.6, 1.9), (0.5, 0, 0)],
    "J044033": ["0355-483", (0.2, 8.2, -0.3), (0.5, 0, 0)],
    "J044737": ["0445-221", (0.25, 3.4, -2.2), (0.5, 0, 0)],
    "J052824": ["0528-250", (0.2, 8, -3.9), (0.5, 0, 0)],
    "J223933": ["2327-459", (0.2, 10, -3.1), (0.5, 0, 0)],
    "J224408": ["2240-260", (0.25, 8.3, -5.2), (0.5, 0, 0)],
    "J215436": ["2211-388", (0.3, 1.1, -0.34), (0.5, 0, 0)],
    "1934_cal_l": ["1934_cal_l", (1, 0, 0), (0.5, 0, 0)],
    "j032237": ["j0307-4857", (0.1, 0.0, 0.0), (0.5, 0.0, 0.0)],
}

# Defining constant variables for all sources
ref = "CA04"
export_pngs = True
sec = source_dict[tar][0]
sourcepar = source_dict[tar][1]
src_dir = f"{data_dir}{tar}"
visname = f"{data_dir}data/atca_2020_{ATCA_band}.ms"
msname = f"{data_dir}data/{tar}_{ATCA_band}.ms"

if ATCA_band == "L":
    n_spw = 8
    pri = "1934_cal_l"
elif ATCA_band == "C":
    n_spw = 5
    pri = "1934_cal_cx"
elif ATCA_band == "X":
    pri = "1934_cal_cx"
    n_spw = 4

print("Here we go! Time to analyse some ATCA data!")
# Uncomment whichever step you don't need to do
# Initial flagging
# process.flag_ms(visname)

# Split to make its own ms
# process.split_ms(
#     src_dir,
#     f"{src_dir}/images",
#     visname,
#     msname,
#     ATCA_band,
#     pri,
#     sec,
#     tar,
#     n_spw,
# )


# # Calibrate, and apply cal ms using primary and secondary
# process.calibrate_ms(src_dir, msname, ATCA_band, ref, pri, sec, tar)
# process.applycal_ms(src_dir, msname, ATCA_band, pri, sec, tar)

# # Post cal inspection and flagging
# process.flagcal_ms(f"{src_dir}/images", msname, ATCA_band, pri, sec)
# process.flagcaltar_ms(src_dir, msname, ATCA_band, pri, sec, tar)

# process.split_imgms(data_dir, tar, "", ATCA_band, n_spw)
imagems = f"2020-_{tar}_{ATCA_band}.ms"
imagename = f"2020-{tar}_{ATCA_band}"

# imagems = [f"{data_dir}data/previous_processing/epoch1_L_{tar}_img.ms", f"{data_dir}data/previous_processing/epoch2_L_{tar}_img.ms", f"{data_dir}data/previous_processing/epoch3_L_{tar}_img.ms", f"{data_dir}data/previous_processing/epoch4_L_{tar}_img.ms"]

# # process.imgmfs_ms(src_dir, imagems, imagename, ATCA_band, n_spw)
# process.img_ms(src_dir, imagems, imagename, ATCA_band, n_spw)
# process.slefcal_ms(src_dir, imagems, imagename, ATCA_band, n_spw)
process.measureflux_ms(
    src_dir, imagems, f"2020_{tar}_{ATCA_band}_postscal.ms", f"2020_{tar}_{ATCA_band}", ATCA_band, sourcepar, n_spw)

if run_epoch == "TRUE":
    for epoch in ["01", "03", "04", "05"]:
        process.split_imgms(data_dir, tar, epoch, ATCA_band, n_spw)
        imagems = f"2020-{epoch}_{tar}_{ATCA_band}.ms"
        imagename = f"2020-{epoch}_{tar}_{ATCA_band}"
        process.imgmfs_ms(src_dir, imagems, imagename, ATCA_band, n_spw)
        process.img_ms(src_dir, imagems, imagename, ATCA_band, n_spw)
        process.slefcal_ms(src_dir, imagems, imagename, ATCA_band, n_spw)
        process.measureflux_ms(
            src_dir, imagems, f"2020-{epoch}_{tar}_{ATCA_band}_postscal.ms", f"2020-{epoch}_{tar}_{ATCA_band}", ATCA_band, sourcepar, n_spw)



# Post image analysis: pbcor, measure flux
# process.pbcor_ms(src_dir, targetms, ATCA_band, n_spw, tar)


# # Post analysis
# if export_pngs is True:
#     process.inspection_plots(
#         src_dir, img_dir, visname, msname, targetms, epoch, ATCA_band, pri, sec, tar
#     )
# process.export_fitspng(src_dir, n_spw, epoch, ATCA_band, tar)

