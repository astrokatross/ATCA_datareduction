#!/usr/bin/python3
# This script is to plot and analyse the ATCA data reduction
# By K.Ross 20/1/21

import pandas as pd
import matplotlib.pyplot as plt
import CFigTools.CustomFigure as CF
import numpy as np
import gpscssmodels

plt.rcParams["axes.grid"] = False
plt.rcParams["font.family"] = "serif"

# Setting the frequencies (delete individual atca bands when processed)
freq_cont = np.linspace(0.01, 11, num=10000)
freq_mwa = [
    0.076,
    0.084,
    0.092,
    0.099,
    0.107,
    0.115,
    0.122,
    0.130,
    0.143,
    0.151,
    0.158,
    0.166,
    0.174,
    0.181,
    0.189,
    0.197,
    0.204,
    0.212,
    0.220,
    0.227,
]

freq_xtra = [
    0.076,
    0.084,
    0.092,
    0.099,
    0.150,
    0.408,
    0.843,
    1.400,
    0.074,
]  # [gleamx4, TGSS, MRC, SUMSS, NVSS, VLSSr]

freq_atca = [
    1.33,
    1.407,
    1.638,
    1.869,
    2.1,
    2.331,
    2.562,
    2.793,
    4.71,
    5.090,
    5.500,
    5.910,
    6.320,
    8.732,
    9.245,
    9.758,
    10.269,
]
xtra_fluxes = [
    "S_076",
    "S_084",
    "S_092",
    "S_099",
    "S_tgss",
    "S_mrc",
    "S_sumss",
    "S_nvss",
    "S_vlssr",
]
gleam_fluxes = ["S_076_err", "S_084_err", "S_092_err", "S_099_err"]
ext = ("yr1", "yr2")
channel = ("69", "93", "121", "145", "169")
for extension in ext:
    if extension == "yr1":
        mwa_yr1_fluxes = (
            f"S_107_{extension}",
            f"S_115_{extension}",
            f"S_123_{extension}",
            f"S_130_{extension}",
            f"S_143_{extension}",
            f"S_150_{extension}",
            f"S_158_{extension}",
            f"S_166_{extension}",
            f"S_174_{extension}",
            f"S_181_{extension}",
            f"S_189_{extension}",
            f"S_197_{extension}",
            f"S_204_{extension}",
            f"S_212_{extension}",
            f"S_220_{extension}",
            f"S_227_{extension}",
        )
        mwa_yr1_errors = (
            f"S_107_err_{extension}",
            f"S_115_err_{extension}",
            f"S_123_err_{extension}",
            f"S_130_err_{extension}",
            f"S_143_err_{extension}",
            f"S_150_err_{extension}",
            f"S_158_err_{extension}",
            f"S_166_err_{extension}",
            f"S_174_err_{extension}",
            f"S_181_err_{extension}",
            f"S_189_err_{extension}",
            f"S_197_err_{extension}",
            f"S_204_err_{extension}",
            f"S_212_err_{extension}",
            f"S_220_err_{extension}",
            f"S_227_err_{extension}",
        )
    if extension == "yr2":
        mwa_yr2_fluxes = (
            f"S_107_{extension}",
            f"S_115_{extension}",
            f"S_123_{extension}",
            f"S_130_{extension}",
            f"S_143_{extension}",
            f"S_150_{extension}",
            f"S_158_{extension}",
            f"S_166_{extension}",
            f"S_174_{extension}",
            f"S_181_{extension}",
            f"S_189_{extension}",
            f"S_197_{extension}",
            f"S_204_{extension}",
            f"S_212_{extension}",
            f"S_220_{extension}",
            f"S_227_{extension}",
        )
        mwa_yr2_errors = (
            f"S_107_err_{extension}",
            f"S_115_err_{extension}",
            f"S_123_err_{extension}",
            f"S_130_err_{extension}",
            f"S_143_err_{extension}",
            f"S_150_err_{extension}",
            f"S_158_err_{extension}",
            f"S_166_err_{extension}",
            f"S_174_err_{extension}",
            f"S_181_err_{extension}",
            f"S_189_err_{extension}",
            f"S_197_err_{extension}",
            f"S_204_err_{extension}",
            f"S_212_err_{extension}",
            f"S_220_err_{extension}",
            f"S_227_err_{extension}",
        )
quadplot_para_yr1 = ["norm_quad_yr1", "alpha_quad_yr1", "quad_curve_yr1"]
quadplot_para_yr2 = ["norm_quad_yr2", "alpha_quad_yr2", "quad_curve_yr2"]
# Reading in the MWA fluxes from master pop
source_dict = {
    "J001513": ["GLEAM J001513-472706"],
    "J015445": ["GLEAM J015445-232950"],
    "J020507": ["GLEAM J020507-110922"],
    "J021246": ["GLEAM J021246-305454"],
    "J022744": ["GLEAM J022744-062106"],
    "J024838": ["GLEAM J024838-321336"],
    "J032213": ["GLEAM J032213-462646"],
    "J032836": ["GLEAM J032836-202138"],
    "J033023": ["GLEAM J033023-074052"],
    "J042502": ["GLEAM J042502-245129"],
    "J044033": ["GLEAM J044033-422918"],
    "J044737": ["GLEAM J044737-220335"],
    "J052824": ["GLEAM J052824-331104"],
    "J223933": ["GLEAM J223933-451414"],
    "J224408": ["GLEAM J224408-202719"],
    "J215436": ["GLEAM J215436-410853"],
}
subchans_dict = {
    "69": ["072-080", "080-088", "088-095", "095-103"],
    "93": ["103-111", "111-118", "118-126", "126-134"],
    "121": ["139-147", "147-154", "154-162", "162-170"],
    "145": ["170-177", "177-185", "185-193", "193-200"],
    "169": ["200-208", "208-216", "216-223", "223-231"],
}


def read_gleam_fluxes(directory, tar):
    master_pop_pd = pd.read_csv(f"{directory}/master_pop_extended.csv")
    name = source_dict[tar][0]
    source_pd = master_pop_pd.query(f"Name=='{name}'")
    vip = np.array(source_pd["VIP"])
    moss = np.array(source_pd["MOSS"])
    mwa_flux_yr1 = np.array(source_pd.loc[:, source_pd.columns.isin(mwa_yr1_fluxes)])[0]
    mwa_flux_yr2 = np.array(source_pd.loc[:, source_pd.columns.isin(mwa_yr2_fluxes)])[0]
    mwa_err_yr1 = np.sqrt(
        (np.array(source_pd.loc[:, source_pd.columns.isin(mwa_yr1_errors)])[0]) ** 2
        + (0.01 * mwa_flux_yr1) ** 2
    )
    mwa_err_yr2 = np.sqrt(
        (np.array(source_pd.loc[:, source_pd.columns.isin(mwa_yr2_errors)])[0]) ** 2
        + (0.01 * mwa_flux_yr2) ** 2
    )

    freq_cont = np.linspace(10, 11000, num=10000)
    para_yr1 = np.array(source_pd.loc[:, source_pd.columns.isin(quadplot_para_yr1)])[0]
    para_yr1 = np.flip(para_yr1)
    yvals_yr1 = gpscssmodels.quad_plot(freq_cont, *para_yr1)
    para_yr2 = np.array(source_pd.loc[:, source_pd.columns.isin(quadplot_para_yr2)])[0]
    para_yr2 = np.flip(para_yr2)
    yvals_yr2 = gpscssmodels.quad_plot(freq_cont, *para_yr2)
    fluxes_extra = np.array(source_pd.loc[:, source_pd.columns.isin(xtra_fluxes)])[0]
    gleam_err = np.sqrt(
        (np.array(source_pd.loc[:, source_pd.columns.isin(gleam_fluxes)])[0]) ** 2
        + (0.01 * fluxes_extra[0:4]) ** 2
    )

    mwa_flux = [mwa_flux_yr1, mwa_flux_yr2, mwa_err_yr1, mwa_err_yr2]
    model_vals = [yvals_yr1, yvals_yr2]
    xtra_values = [fluxes_extra, gleam_err]
    values = [name, vip, moss]
    return mwa_flux, model_vals, xtra_values, values


def read_mwa_fluxes(directory, tar, epoch):
    name = source_dict[tar][0]
    mwa_flux = []
    mwa_errs = []
    for i in range(len(channel)):
        subchans = subchans_dict[channel[i]]
        chan = channel[i]
        for subchan in subchans:
            try:
                src_mwa_pd = pd.read_csv(
                    f"{directory}/{epoch}/{chan}/minimosaic/{tar}_{subchan}MHz_scaled_comp_xmatch.csv"
                )
                src_pd = src_mwa_pd.query(f"Name=='{name}'")
                mwa_flux_chan = np.array(src_pd["int_flux"])[0]
                mwa_errs_chan = np.sqrt(
                    (
                        np.array(src_pd["err_int_flux"]) ** 2
                        + (0.01 * mwa_flux_chan) ** 2
                    )
                )[0]
                mwa_flux.append(mwa_flux_chan)
                mwa_errs.append(mwa_errs_chan)
            except (FileNotFoundError, KeyError):
                print(
                    f"{directory}/{epoch}/{chan}/minimosaic/{tar}_{subchan}MHz_scaled_comp_xmatch.csv not found!"
                )
                mwa_flux.append(np.nan)
                mwa_errs.append(np.nan)
            except:
                print(f"{tar} {subchan} not found in catalogue! Maybe too faint")
                mwa_flux.append(np.nan)
                mwa_errs.append(np.nan)

    return mwa_flux, mwa_errs


def read_atca_fluxes(directory, tar):
    atca_fluxes = [
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
    ]
    epochs = ("epoch1", "epoch2", "epoch3", "epoch4", "epoch5")
    for i in range(0, 4):
        epoch = epochs[i]
        try:
            atca_Lband_pd = pd.read_csv(f"{directory}/{tar}/{tar}_{epoch}_L.csv")
            atca_Lband_epoch = np.array(atca_Lband_pd["# S_Lband"])
        except FileNotFoundError:
            print(f"No L-Band for {epoch}")
            atca_Lband_epoch = [[np.nan] * 8]

        try:
            atca_Cband_pd = pd.read_csv(f"{directory}/{tar}/{tar}_{epoch}_C.csv")
            atca_Cband_epoch = np.array(atca_Cband_pd["# S_Cband"])
        except FileNotFoundError:
            print(f"No C-Band for {epoch}")
            atca_Cband_epoch = [[np.nan] * 5]

        try:
            atca_Xband_pd = pd.read_csv(f"{directory}/{tar}/{tar}_{epoch}_X.csv")
            atca_Xband_epoch = np.array(atca_Xband_pd["# S_Xband"])
        except FileNotFoundError:
            print(f"No X-Band for {epoch}")
            atca_Xband_epoch = [[np.nan] * 4]

        atca_epoch = np.concatenate(
            (atca_Lband_epoch, atca_Cband_epoch, atca_Xband_epoch), axis=None
        )
        atca_epoch[np.where(atca_epoch < 0.0)] = np.nan
        atca_fluxes[i] = atca_epoch
    return atca_fluxes


def plt_sed(
    directory,
    tar,
    MWA=True,
    extra_surveys=True,
    save_fig=True,
):
    colours = ("red", "mediumseagreen", "C9", "C4", "C4")
    epoch_nms = ("Epoch1", "Epoch2", "Epoch3", "Epoch4", "Epoch5")
    mwa_flux, model_vals, xtra_values, values = read_gleam_fluxes("/data/MWA", tar)
    mwa_fluxes_e3, mwa_errors_e3 = read_mwa_fluxes("/data/MWA", tar, "epoch3")
    atca_fluxes = read_atca_fluxes(directory, tar)
    f = CF.sed_fig()
    if MWA is True:
        f.plot_spectrum(
            freq_mwa[4:20],
            mwa_flux[0],
            mwa_flux[2],
            marker="o",
            label="MWAYr1",
            marker_color="C6",
        )
        f.plot_spectrum(
            freq_mwa[4:20],
            mwa_flux[1],
            mwa_flux[3],
            marker="o",
            label="MWAYr2",
            marker_color="mediumblue",
        )
        f.plot_residuals(
            freq_mwa[4:20],
            mwa_flux[0],
            mwa_flux[1],
            mwa_flux[2],
            mwa_flux[3],
        )
        f.display_model(
            freq_cont,
            model_vals[0],
            color="C6",
            label=None,
            model_min=None,
            model_max=None,
            alpha_patch=0.2,
        )
        f.display_model(
            freq_cont,
            model_vals[1],
            color="mediumblue",
            label=None,
            model_min=None,
            model_max=None,
            alpha_patch=0.2,
        )
    print(mwa_fluxes_e3)
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e3,
        mwa_errors_e3,
        marker="o",
        marker_color=colours[2],
    )
    for i in range(0, 4):
        if i == 1:
            if tar in ["J001513"]:
                freq_atca = [
                    1.33,
                    1.407,
                    1.638,
                    1.869,
                    2.1,
                    2.331,
                    2.562,
                    2.793,
                    4.71,
                    5.090,
                    5.500,
                    5.910,
                    6.320,
                    8.732,
                    9.245,
                    9.758,
                ]
            else:
                freq_atca = [
                    1.33,
                    1.407,
                    1.638,
                    1.869,
                    2.1,
                    2.331,
                    2.562,
                    2.793,
                    4.71,
                    5.090,
                    5.500,
                    5.910,
                    6.320,
                    8.732,
                    9.245,
                    9.758,
                    10.269,
                ]
        else:
            freq_atca = [
                1.33,
                1.407,
                1.638,
                1.869,
                2.1,
                2.331,
                2.562,
                2.793,
                4.71,
                5.090,
                5.500,
                5.910,
                6.320,
                8.732,
                9.245,
                9.758,
                10.269,
            ]
        f.plot_spectrum(
            freq_atca,
            atca_fluxes[i],
            atca_fluxes[i] * 0.03,
            marker="o",
            marker_color=colours[i],
            label=epoch_nms[i],
        )

    if extra_surveys is True:
        fluxes_extra = xtra_values[0]
        gleam_err = xtra_values[1]
        f.plot_spectrum(
            freq_xtra[0:4],
            fluxes_extra[0:4],
            gleam_err,
            marker="o",
            marker_color="C6",
        )
        if fluxes_extra[8] != np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[8],
                fluxes_extra[8],
                marker="X",
                label="VLSSr",
                marker_color="red",
            )
        if fluxes_extra[4] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[4],
                fluxes_extra[4] * 0.001,
                marker="s",
                label="TGSS",
                marker_color="k",
            )
        if fluxes_extra[5] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[5],
                fluxes_extra[5],
                marker="*",
                label="MRC",
                marker_color="green",
            )
        if fluxes_extra[6] == np.nan:
            print("No SUMSS")
            pass
        else:
            f.plot_point(
                freq_xtra[6],
                fluxes_extra[6] * 0.001,
                marker="p",
                label="SUMSS",
                marker_color="orange",
            )
        if fluxes_extra[7] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[7],
                fluxes_extra[7] * 0.001,
                marker="D",
                label="NVSS",
                marker_color="navy",
            )

    f.legend(values[1], values[2], loc="lower center")
    f.title(values[0])
    f.format(xunit="GHz")

    if save_fig is True:
        f.save(f"{directory}/SEDs/{tar}_sed", ext="png")
    return
