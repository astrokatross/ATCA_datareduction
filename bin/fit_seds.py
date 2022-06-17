#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import analysis_functs
import numpy as np
import matplotlib.pyplot as plt

# Source/run information
save_dir = "/data/ATCA/analysis/"
data_dir = "/data/ATCA/ATCA_datareduction/"
gleam_targets = [
    # "GLEAM J215436-410853",
    "GLEAM J001513-472706",
    "GLEAM J015445-232950",
    "GLEAM J020507-110922",
    # "GLEAM J021246-305454",
    "GLEAM J022744-062106",
    # "GLEAM J024838-321336",
    # "GLEAM J032213-462646",
    # "GLEAM J032836-202138",
    # "GLEAM J033023-074052",
    # "GLEAM J042502-245129",
    # "GLEAM J044033-422918",
    # "GLEAM J044737-220335",
    # "GLEAM J052824-331104",
    # "GLEAM J223933-451414",
    "GLEAM J224408-202719",
    # "GLEAM J215436-410853",
]

avg_logz = {}
for i in range(len(gleam_targets)):
    gleam_tar = gleam_targets[i]
    target = gleam_tar.strip("GLEAM ")[0:7]
    avg_logz_src = analysis_functs.run_everything(save_dir, data_dir, gleam_tar)
    avg_logz_src = np.array(avg_logz_src)
    avg_logz[gleam_tar] = np.around(avg_logz_src, decimals=1)
    if target == "J215436":
        start_times = [
            "09:37:44",
            "10:41:28",
            # "10:54:16",
            # "11:41:12",
            "11:54:00",
            # "12:40:56",
            "12:53:44",
            # "13:40:40",
            # "13:53:28",
            "14:40:48",
            "14:53:52",
        ]
        end_times = [
            "09:47:13",
            "10:50:47",
            # "11:03:45",
            # "11:50:41",
            "12:03:29",
            # "12:50:25",
            "13:03:13",
            # "13:50:09",
            # "14:02:57",
            "14:50:17",
            "15:03:21",
        ]
        timeranges_215436 = analysis_functs.read_timeranges(start_times, end_times)
        scan_times_arr = np.arange(0, 10, 0.5)
        scan_times_arr2 = np.arange(64, 74, 0.5)
        # scan_times_arr3 = np.arange(77, 87, 0.5)
        # scan_times_arr4 = np.arange(124, 134, 0.5)
        scan_times_arr5 = np.arange(137, 147, 0.5)
        # scan_times_arr6 = np.arange(183, 193, 0.5)
        scan_times_arr7 = np.arange(196, 206, 0.5)
        scan_times_arr8 = np.arange(243, 253, 0.5)
        scan_times_arr9 = np.arange(256, 266, 0.5)
        scan_times_arr10 = np.arange(303, 313, 0.5)
        scan_times_arr11 = np.arange(316, 326, 0.5)
        scan_times = np.concatenate(
            (
                scan_times_arr,
                scan_times_arr2,
                # scan_times_arr3,
                # scan_times_arr4,
                scan_times_arr5,
                # scan_times_arr6,
                scan_times_arr7,
                # scan_times_arr8,
                # scan_times_arr9,
                scan_times_arr10,
                scan_times_arr11,
            )
        )
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J215436/casa_files/{target}_X_epoch5"
        )
        (
            fluxes_j215436_x,
            err_fluxes_j215436_x,
            mod_j215436_x,
        ) = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/epoch5_X_{target}", outfile_dir, timeranges_215436
        )
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J215436/casa_files/{target}_C_epoch5"
        )
        (
            fluxes_j215436_c,
            err_fluxes_j215436_c,
            mod_j215436_c,
        ) = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/epoch5_C_{target}", outfile_dir, timeranges_215436
        )
        fluxes_j215436 = np.stack((fluxes_j215436_c, fluxes_j215436_x))
        err_fluxes_j215436 = np.stack((err_fluxes_j215436_c, err_fluxes_j215436_x))
        mod_j215436 = np.stack((mod_j215436_c, mod_j215436_x))
        # print(len(timeranges_215436))

        start_times_sec = [
            "09:35:20",
            "09:48:08",
            "10:39:04",
            "10:51:52",
            # "11:04:40",
            # "11:38:48",
            "11:51:36",
            "12:04:24",
            # "12:38:32",
            "12:51:20",
            "13:04:08",
            # "13:38:16",
            # "13:51:04",
            # "14:03:52",
            "14:38:16",
            "14:51:20",
            "15:04:32",
        ]
        end_times_sec = [
            "09:36:49",
            "09:49:37",
            "10:40:33",
            "10:53:21",
            # "11:06:09",
            # "11:40:17",
            "11:53:05",
            "12:05:53",
            # "12:40:01",
            "12:52:49",
            "13:05:37",
            # "13:39:45",
            # "13:52:33",
            # "14:05:21",
            "14:39:45",
            "14:52:49",
            "15:06:01",
        ]
        timeranges_215436_sec = analysis_functs.read_timeranges(
            start_times_sec, end_times_sec
        )
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J215436/casa_files/2211-388_X_epoch5"
        )
        (
            fluxes_j215436_sec_x,
            err_fluxes_j215436_sec_x,
            mod_j215436_sec_x,
        ) = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/epoch5_X_2211-388", outfile_dir, timeranges_215436_sec
        )
        scan_times_sec = np.concatenate(
            (
                np.arange(0, 2, 0.5),
                np.arange(13, 15, 0.5),
                np.arange(64, 66, 0.5),
                np.arange(76, 78, 0.5),
                # np.arange(89, 91, 0.5),
                # np.arange(123, 125, 0.5),
                np.arange(136, 138, 0.5),
                np.arange(149, 151, 0.5),
                # np.arange(183, 185, 0.5),
                np.arange(196, 198, 0.5),
                np.arange(209, 211, 0.5),
                # np.arange(243, 245, 0.5),
                # np.arange(256, 258, 0.5),
                # np.arange(268, 270, 0.5),
                np.arange(303, 305, 0.5),
                np.arange(316, 318, 0.5),
                np.arange(329, 331, 0.5),
            )
        )
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J215436/casa_files/2211-388_C_epoch5"
        )
        (
            fluxes_j215436_sec_c,
            err_fluxes_j215436_sec_c,
            mod_j215436_sec_c,
        ) = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/epoch5_C_2211-388", outfile_dir, timeranges_215436_sec
        )

        fluxes_j215436_sec = np.stack((fluxes_j215436_sec_c, fluxes_j215436_sec_x))
        err_fluxes_j215436_sec = np.stack(
            (err_fluxes_j215436_sec_c, err_fluxes_j215436_sec_x)
        )
        mod_j215436_sec = np.stack((mod_j215436_sec_c, mod_j215436_sec_x))
        # analysis_functs.plt_lightcurve_continual(
        #     f"{save_dir}Plots/",
        #     scan_times,
        #     fluxes_j215436,
        #     err_fluxes_j215436,
        #     mod_j215436,
        #     scan_times_sec,
        #     fluxes_j215436_sec,
        #     err_fluxes_j215436_sec,
        #     mod_j215436_sec,
        #     ext="png",
        # )
    if target == "J001513":
        plt.close()
        plt.clf()
        start_times = ["17:53:00", "19:31:30"]
        end_times = ["18:02:20", "19:40:50"]
        timeranges_001513 = analysis_functs.read_timeranges(start_times, end_times)
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_C_2021-10-15"
        )
        (
            fluxes_j001513_c,
            err_fluxes_j001513_c,
            mod_j001513_c,
        ) = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/2021-10-15_C_{target}", outfile_dir, timeranges_001513
        )
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_X_2021-10-15"
        )
        (
            fluxes_j001513_x,
            err_fluxes_j001513_x,
            mod_j001513_x,
        ) = analysis_functs.read_lightcurveflux(
            f"{data_dir}data/2021-10-15_X_{target}", outfile_dir, timeranges_001513
        )
        fluxes_j001513 = np.stack((fluxes_j001513_c, fluxes_j001513_x))
        err_fluxes_j001513 = np.stack((err_fluxes_j001513_c, err_fluxes_j001513_x))
        mod_j001513 = np.stack((mod_j001513_c, mod_j001513_x))
    elif target == "J020507":
        start_times = ["17:01:50", "18:46:10"]
        end_times = ["17:11:20", "18:55:40"]
        timeranges_020507 = analysis_functs.read_timeranges(start_times, end_times)
        outfile_dir = (
            f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_C_2021-10-15"
        )
        # (
        #     fluxes_j020507_c,
        #     err_fluxes_j020507_c,
        #     mod_j020507_c,
        # ) = analysis_functs.read_lightcurveflux(
        #     f"{data_dir}data/2021-10-15_C_{target}", outfile_dir, timeranges_020507
        # )
        # outfile_dir = (
        #     f"/data/ATCA/ATCA_datareduction/J020507/casa_files/{target}_X_2021-10-15"
        # )
        # (
        #     fluxes_j020507_x,
        #     err_fluxes_j020507_x,
        #     mod_j020507_x,
        # ) = analysis_functs.read_lightcurveflux(
        #     f"{data_dir}data/2021-10-15_X_{target}", outfile_dir, timeranges_020507
        # )
        # fluxes_j020507 = np.stack((fluxes_j020507_c, fluxes_j020507_x))
        # err_fluxes_j020507 = np.stack((err_fluxes_j020507_c, err_fluxes_j020507_x))
        # mod_j020507 = np.stack((mod_j020507_c, mod_j020507_x))
        # scan_times_arr = np.arange(0, 10, 0.5)
        # scan_times_0015_2 = np.arange(98.5, 108.5, 0.5)
        # scan_times_0205_2 = np.arange(104.3333333, 114.3333333, 0.5)
        # scan_times_src1 = np.stack((scan_times_arr, scan_times_0015_2))
        # scan_times_src2 = np.stack((scan_times_arr, scan_times_0205_2))
        # print("plotting!!!!")
        # analysis_functs.plt_lightcurve(
        #     f"{save_dir}Plots/",
        #     scan_times_src1,
        #     scan_times_src2,
        #     fluxes_j001513,
        #     fluxes_j020507,
        #     err_fluxes_j001513,
        #     err_fluxes_j020507,
        #     mod_j001513,
        #     mod_j020507,
        #     # ext="png",
        # )


print(avg_logz)
