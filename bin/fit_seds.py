#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import analysis_functs
import numpy as np
import datetime
import matplotlib.pyplot as plt 
import math 

# Source/run information
save_dir = "/data/ATCA/analysis/"
data_dir = "/data/ATCA/ATCA_datareduction/"
gleam_targets = [
    "GLEAM J001513-472706",
    # "GLEAM J015445-232950",
    "GLEAM J020507-110922",
    "GLEAM J021246-305454",
    "GLEAM J022744-062106",
    "GLEAM J024838-321336",
    "GLEAM J032213-462646",
    "GLEAM J032836-202138",
    "GLEAM J033023-074052",
    "GLEAM J042502-245129",
    "GLEAM J044033-422918",
    "GLEAM J044737-220335",
    "GLEAM J052824-331104",
    "GLEAM J223933-451414",
    "GLEAM J224408-202719",
]


for i in range(len(gleam_targets)):
    gleam_tar = gleam_targets[i]
    target = gleam_tar.strip("GLEAM ")[0:7]
    nu_p, err_nu_p = analysis_functs.run_everything(save_dir, data_dir, gleam_tar)

    if target == "J001513":
        plt.close()
        plt.clf()
        start_times = ["17:53:00", "19:31:30"]
        end_times = ["18:02:20", "19:40:50"]
        timeranges_001513 = analysis_functs.read_timeranges(start_times, end_times)
        outfile_dir = f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_C_2021-10-15"
        fluxes_j001513_c, err_fluxes_j001513_c, mod_j001513_c = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/2021-10-15_C_{target}", outfile_dir, timeranges_001513
        )
        outfile_dir = f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_X_2021-10-15"
        fluxes_j001513_x, err_fluxes_j001513_x, mod_j001513_x = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/2021-10-15_X_{target}", outfile_dir, timeranges_001513
        )
        fluxes_j001513 = np.stack((fluxes_j001513_c, fluxes_j001513_x))
        err_fluxes_j001513 = np.stack((err_fluxes_j001513_c, err_fluxes_j001513_x))
        mod_j001513 = np.stack((mod_j001513_c, mod_j001513_x))
    elif target == "J020507":
        start_times = ["17:01:50", "18:46:10"]
        end_times = ["17:11:20", "18:55:40"]
        timeranges_020507 = analysis_functs.read_timeranges(start_times, end_times)
        outfile_dir = f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_C_2021-10-15"
        fluxes_j020507_c, err_fluxes_j020507_c, mod_j020507_c = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/2021-10-15_C_{target}", outfile_dir, timeranges_020507
        )
        outfile_dir = f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_X_2021-10-15"
        fluxes_j020507_x, err_fluxes_j020507_x, mod_j020507_x = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/2021-10-15_X_{target}", outfile_dir, timeranges_020507
        )
        fluxes_j020507 = np.stack((fluxes_j020507_c, fluxes_j020507_x))
        err_fluxes_j020507 = np.stack((err_fluxes_j020507_c, err_fluxes_j020507_x))
        mod_j020507 = np.stack((mod_j020507_c, mod_j020507_x))
        scan_times_arr = np.arange(0, 10, 0.5)
        scan_times_0015_2 = np.arange(98.5, 108.5, 0.5)
        scan_times_0205_2 = np.arange(104.3333333, 114.3333333, 0.5)
        scan_times_src1 = np.stack((scan_times_arr, scan_times_0015_2))
        scan_times_src2 = np.stack((scan_times_arr, scan_times_0205_2))
        print("plotting!!!!")
        analysis_functs.plt_lightcurve(
            f"{save_dir}Plots/",
            scan_times_src1,
            scan_times_src2,
            fluxes_j001513,
            fluxes_j020507,
            err_fluxes_j001513,
            err_fluxes_j020507,
            mod_j001513,
            mod_j020507,
        )
