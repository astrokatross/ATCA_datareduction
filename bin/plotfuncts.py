#!/usr/bin/python3
# This script is to plot and analyse the ATCA data reduction
# By K.Ross 20/1/21

import pandas as pd
import matplotlib.pyplot as plt
import CFigTools.CustomFigure as CF
import numpy as np
import gpscssmodels
import os 

plt.rcParams["axes.grid"] = False
plt.rcParams["font.family"] = "serif"

# Setting the frequencies (delete individual atca bands when processed)
freq_atca = [
    1.200,
    1.454,
    1.711,
    1.968,
    2.225,
    2.482,
    2.739,
    2.996,
    4.680,
    5.090,
    5.500,
    5.910,
    6.320,
    8.732,
    9.245,
    9.758,
    10.269,
]
freq_mwa = [
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
freq_mwa_e3 = [0.143,0.15,0.158,0.166,0.174,0.181,0.189,0.197,0.204,0.212,0.22,0.227]
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
# freq_lband = [1.310, 1.425, 1.725, 1.968, 2.225, 2.482, 2.739, 2.95]
# This is freequency for cut overall Lband 
freq_lband = [1.35, 1.483, 1.754, 1.985, 2.216, 2.447, 2.678, 2.909]
freq_cband = [4.680, 5.090, 5.500, 5.910, 6.320]
freq_xband = [8.732, 9.245, 9.758, 10.269]
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
mwa_yr1_fluxes = (
    "S_107_yr1",
    "S_115_yr1",
    "S_123_yr1",
    "S_130_yr1",
    "S_143_yr1",
    "S_150_yr1",
    "S_158_yr1",
    "S_166_yr1",
    "S_174_yr1",
    "S_181_yr1",
    "S_189_yr1",
    "S_197_yr1",
    "S_204_yr1",
    "S_212_yr1",
    "S_220_yr1",
    "S_227_yr1",
)

mwa_yr2_fluxes = (
    "S_107_yr2",
    "S_115_yr2",
    "S_123_yr2",
    "S_130_yr2",
    "S_143_yr2",
    "S_150_yr2",
    "S_158_yr2",
    "S_166_yr2",
    "S_174_yr2",
    "S_181_yr2",
    "S_189_yr2",
    "S_197_yr2",
    "S_204_yr2",
    "S_212_yr2",
    "S_220_yr2",
    "S_227_yr2",
)
mwa_yr1_errors = (
    "S_107_err_yr1",
    "S_115_err_yr1",
    "S_123_err_yr1",
    "S_130_err_yr1",
    "S_143_err_yr1",
    "S_150_err_yr1",
    "S_158_err_yr1",
    "S_166_err_yr1",
    "S_174_err_yr1",
    "S_181_err_yr1",
    "S_189_err_yr1",
    "S_197_err_yr1",
    "S_204_err_yr1",
    "S_212_err_yr1",
    "S_220_err_yr1",
    "S_227_err_yr1",
)
mwa_yr2_errors = (
    "S_107_err_yr2",
    "S_115_err_yr2",
    "S_123_err_yr2",
    "S_130_err_yr2",
    "S_143_err_yr2",
    "S_150_err_yr2",
    "S_158_err_yr2",
    "S_166_err_yr2",
    "S_174_err_yr2",
    "S_181_err_yr2",
    "S_189_err_yr2",
    "S_197_err_yr2",
    "S_204_err_yr2",
    "S_212_err_yr2",
    "S_220_err_yr2",
    "S_227_err_yr2",
)
mwa_epoch3_fluxes = ("S_143_e3",
    "S_150_e3",
    "S_158_e3",
    "S_166_e3",
    "S_174_e3",
    "S_181_e3",
    "S_189_e3",
    "S_197_e3",
    "S_204_e3",
    "S_212_e3",
    "S_220_e3",
    "S_227_e3")
mwa_epoch3_errors = ("err_S_143_e3",
    "err_S_150_e3",
    "err_S_158_e3",
    "err_S_166_e3",
    "err_S_174_e3",
    "err_S_181_e3",
    "err_S_189_e3",
    "err_S_197_e3",
    "err_S_204_e3",
    "err_S_212_e3",
    "err_S_220_e3",
    "err_S_227_e3")
quadplot_paras_yr1 = ["norm_quad_yr1", "alpha_quad_yr1", "quad_curve_yr1"]
quadplot_paras_yr2 = ["norm_quad_yr2", "alpha_quad_yr2", "quad_curve_yr2"]
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
    "J215436": ["GLEAM J215436-410853"]
}


def plt_sed(
    directory,
    tar,
    epoch1=True,
    epoch2=True,
    epoch3=True,
    epoch4=True,
    MWA=True,
    extra_surveys=False,
    models=False,
    save_fig=True,
):

    f = CF.sed_fig()
    name = source_dict[tar][0]
    master_pop_pd = pd.read_csv(f"{directory}/data/epoch3_cat.csv")
    if MWA is True:
        mwa_flux_yr1 = np.array(
            master_pop_pd.query(f"Name=='{name}'").loc[
                :, master_pop_pd.columns.isin(mwa_yr1_fluxes)
            ]
        )[0]
        mwa_flux_yr2 = np.array(
            master_pop_pd.query(f"Name=='{name}'").loc[
                :, master_pop_pd.columns.isin(mwa_yr2_fluxes)
            ]
        )[0]
        mwa_err_yr1 = np.sqrt(
            (
                np.array(
                    master_pop_pd.query(f"Name=='{name}'").loc[
                        :, master_pop_pd.columns.isin(mwa_yr1_errors)
                    ]
                )[0]
            )
            ** 2
            + (0.01 * mwa_flux_yr1) ** 2
        )
        mwa_err_yr2 = np.sqrt(
            (
                np.array(
                    master_pop_pd.query(f"Name=='{name}'").loc[
                        :, master_pop_pd.columns.isin(mwa_yr2_errors)
                    ]
                )[0]
            )
            ** 2
            + (0.01 * mwa_flux_yr2) ** 2
        )
        f.plot_spectrum(
            freq_mwa,
            mwa_flux_yr1,
            mwa_err_yr1,
            marker="o",
            label="MWAYr1",
            marker_color="C6",
        )
        f.plot_spectrum(
            freq_mwa,
            mwa_flux_yr2,
            mwa_err_yr2,
            marker="o",
            label="MWAYr2",
            marker_color="mediumblue",
        )
        f.plot_residuals(
            freq_mwa,
            mwa_flux_yr1,
            mwa_flux_yr2,
            mwa_err_yr1,
            mwa_err_yr2,
        )
        mwa_epoch3 = np.array(
            master_pop_pd.query(f"Name=='{name}'").loc[
                :, master_pop_pd.columns.isin(mwa_epoch3_fluxes)
            ]
        )[0]
        err_mwa_epoch3 = np.sqrt(
            (
                np.array(
                    master_pop_pd.query(f"Name=='{name}'").loc[
                        :, master_pop_pd.columns.isin(mwa_epoch3_errors)
                    ]
                )[0]
            )
            ** 2
            + (0.01 * mwa_epoch3) ** 2
        )
        f.plot_spectrum(
            freq_mwa_e3,
            mwa_epoch3,
            err_mwa_epoch3,
            marker="o",
            marker_color="C9",
        )

    # PLOTTING THE L BAND DATA!!!

    try:
        atca_Lband_e1_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch1_L.csv")
        atca_Lband_e1 = np.array(atca_Lband_e1_pd["# S_Lband"])
        err_atca_Lband_e1 = 0.03 * (atca_Lband_e1)
        for i in range(len(atca_Lband_e1)):
            if atca_Lband_e1[i] <=0:
                atca_Lband_e1[i]= np.nan
        f.plot_spectrum(
        freq_lband,
        atca_Lband_e1,
        err_atca_Lband_e1,
        marker="o",
        label="Epoch1",
        marker_color="red",
        )
    except:
        print("No L-Band for Epoch1")

    # try:
    #     atca_Lband_e2_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch2_L.csv")
    #     atca_Lband_e2 = np.array(atca_Lband_e2_pd["# S_Lband"])
    #     err_atca_Lband_e2 = 0.03 * (atca_Lband_e2)
    #     for i in range(len(atca_Lband_e2)):
    #         if atca_Lband_e2[i] <=0:
    #             atca_Lband_e2[i]= np.nan
    #     f.plot_spectrum(
    #     freq_lband,
    #     atca_Lband_e2,
    #     err_atca_Lband_e2,
    #     marker="o",
    #     label="Epoch2",
    #     marker_color="mediumseagreen",
    #     )
    # except:
    #     print("No L-Band for Epoch2")


    try:
        atca_Lband_e3_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch3_L.csv")
        atca_Lband_e3 = np.array(atca_Lband_e3_pd["# S_Lband"])
        err_atca_Lband_e3 = 0.03 * (atca_Lband_e3)
        for i in range(len(atca_Lband_e3)):
            if atca_Lband_e3[i] <=0:
                atca_Lband_e3[i]= np.nan
        f.plot_spectrum(
            freq_lband,
            atca_Lband_e3,
            err_atca_Lband_e3,
            marker="o",
            label="Epoch3",
            marker_color="C9",
        )
    except:
        print("No L-Band for Epoch3")


    try:
        atca_Lband_e4_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch4_L.csv")
        atca_Lband_e4 = np.array(atca_Lband_e4_pd["# S_Lband"])
        err_atca_Lband_e4 = 0.03 * (atca_Lband_e4)
        for i in range(len(atca_Lband_e4)):
            if atca_Lband_e4[i] <=0:
                atca_Lband_e4[i]= np.nan
        f.plot_spectrum(
            freq_lband,
            atca_Lband_e4,
            err_atca_Lband_e4,
            marker="o",
            label="Epoch4",
            marker_color="C4",
        )
    except:
        print("No L-Band for Epoch4")



    # PLOTTING THE C BAND DATA!!!
    try:
        atca_Cband_e1_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch1_C.csv")
        atca_Cband_e1 = np.array(atca_Cband_e1_pd["# S_Cband"])
        err_atca_Cband_e1 = 0.03 * (atca_Cband_e1)
        for i in range(len(atca_Cband_e1)):
            if atca_Cband_e1[i] <=0:
                atca_Cband_e1[i]= np.nan
        f.plot_spectrum(
        freq_cband,
        atca_Cband_e1,
        err_atca_Cband_e1,
        marker="o",
        marker_color="red",
        )
    except:
        print("No C-Band for Epoch1")

    try:
        atca_Cband_e2_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch2_C.csv")
        atca_Cband_e2 = np.array(atca_Cband_e2_pd["# S_Cband"])
        err_atca_Cband_e2 = 0.03 * (atca_Cband_e2)
        for i in range(len(atca_Cband_e2)):
            if atca_Cband_e2[i] <=0:
                atca_Cband_e2[i]= np.nan
        f.plot_spectrum(
        freq_cband,
        atca_Cband_e2,
        err_atca_Cband_e2,
        marker="o",
        marker_color="mediumseagreen",
        )
    except:
        print("No C-Band for Epoch2")


    try:
        atca_Cband_e3_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch3_C.csv")
        atca_Cband_e3 = np.array(atca_Cband_e3_pd["# S_Cband"])
        err_atca_Cband_e3 = 0.03 * (atca_Cband_e3)
        for i in range(len(atca_Cband_e3)):
            if atca_Cband_e3[i] <=0:
                atca_Cband_e3[i]= np.nan
        f.plot_spectrum(
            freq_cband,
            atca_Cband_e3,
            err_atca_Cband_e3,
            marker="o",
            marker_color="C9",
        )
    except:
        print("No C-Band for Epoch3")


    try:
        atca_Cband_e4_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch4_C.csv")
        atca_Cband_e4 = np.array(atca_Cband_e4_pd["# S_Cband"])
        err_atca_Cband_e4 = 0.03 * (atca_Cband_e4)
        for i in range(len(atca_Cband_e4)):
            if atca_Cband_e4[i] <=0:
                atca_Cband_e4[i]= np.nan
        f.plot_spectrum(
            freq_cband,
            atca_Cband_e4,
            err_atca_Cband_e4,
            marker="o",
            marker_color="C4",
        )
    except:
        print("No C-Band for Epoch4")



    # PLOTTING THE X BAND DATA!!!
    try:
        atca_Xband_e1_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch1_X.csv")
        atca_Xband_e1 = np.array(atca_Xband_e1_pd["# S_Xband"])
        err_atca_Xband_e1 = 0.03 * (atca_Xband_e1)
        for i in range(len(atca_Xband_e1)):
            if atca_Xband_e1[i] <=0:
                atca_Xband_e1[i]= np.nan
        f.plot_spectrum(
        freq_xband,
        atca_Xband_e1,
        err_atca_Xband_e1,
        marker="o",
        marker_color="red",
        )
    except:
        print("No X-Band for Epoch1")

    try:
        atca_Xband_e2_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch2_X.csv")
        atca_Xband_e2 = np.array(atca_Xband_e2_pd["# S_Xband"])
        err_atca_Xband_e2 = 0.03 * (atca_Xband_e2)
        for i in range(len(atca_Xband_e2)):
            if atca_Xband_e2[i] <=0:
                atca_Xband_e2[i]= np.nan
        f.plot_spectrum(
        freq_xband,
        atca_Xband_e2,
        err_atca_Xband_e2,
        marker="o",
        marker_color="mediumseagreen",
        )
    except:
        print("No X-Band for Epoch2")


    try:
        atca_Xband_e3_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch3_X.csv")
        atca_Xband_e3 = np.array(atca_Xband_e3_pd["# S_Xband"])
        err_atca_Xband_e3 = 0.03 * (atca_Xband_e3)
        for i in range(len(atca_Xband_e3)):
            if atca_Xband_e3[i] <=0:
                atca_Xband_e3[i]= np.nan
        f.plot_spectrum(
            freq_xband,
            atca_Xband_e3,
            err_atca_Xband_e3,
            marker="o",
            marker_color="C9",
        )
    except:
        print("No X-Band for Epoch3")


    try:
        atca_Xband_e4_pd = pd.read_csv(f"{directory}/{tar}/{tar}_epoch4_X.csv")
        atca_Xband_e4 = np.array(atca_Xband_e4_pd["# S_Xband"])
        err_atca_Xband_e4 = 0.03 * (atca_Xband_e4)
        for i in range(len(atca_Xband_e4)):
            if atca_Xband_e4[i] <=0:
                atca_Xband_e4[i]= np.nan
        f.plot_spectrum(
            freq_xband,
            atca_Xband_e4,
            err_atca_Xband_e4,
            marker="o",
            marker_color="C4",
        )
    except:
        print("No X-Band for Epoch4")


    if extra_surveys is True:
        fluxes_extra = np.array(
            master_pop_pd.query(f"Name=='{name}'").loc[
                :, master_pop_pd.columns.isin(xtra_fluxes)
            ]
        )[0]
        gleam_err = np.sqrt(
            (
                np.array(
                    master_pop_pd.query(f"Name=='{name}'").loc[
                        :, master_pop_pd.columns.isin(gleam_fluxes)
                    ]
                )[0]
            )
            ** 2
            + (0.01 * fluxes_extra[0:4]) ** 2
        )

        f.plot_spectrum(
            freq_xtra[0:4],
            fluxes_extra[0:4],
            gleam_err,
            marker="o",
            label="GLEAM",
            marker_color="grey",
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
                fluxes_extra[4]*0.001,
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
                fluxes_extra[6]*0.001,
                marker="p",
                label="SUMSS",
                marker_color="orange",
            )
        if fluxes_extra[7] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[7],
                fluxes_extra[7]*0.001,
                marker="D",
                label="NVSS",
                marker_color="navy",
            )

    if models is True:
        freq_cont = np.linspace(10, 11000, num=10000)
        paras_yr1 = np.array(
            master_pop_pd.query(f"Name=='{name}'").loc[
                :, master_pop_pd.columns.isin(quadplot_paras_yr1)
            ]
        )[0]
        paras_yr1 = np.flip(paras_yr1)
        yvals_yr1 = gpscssmodels.quad_plot(freq_cont, *paras_yr1)
        freq_cont = freq_cont * 0.001
        f.display_model(
            freq_cont,
            yvals_yr1,
            color="C6",
            label=None,
            model_min=None,
            model_max=None,
            alpha_patch=0.2,
        )
        paras_yr2 = np.array(
            master_pop_pd.query(f"Name=='{name}'").loc[
                :, master_pop_pd.columns.isin(quadplot_paras_yr2)
            ]
        )[0]
        paras_yr2 = np.flip(paras_yr2)
        freq_cont = np.linspace(10, 11000, num=10000)
        yvals_yr2 = gpscssmodels.quad_plot(freq_cont, *paras_yr2)
        freq_cont = freq_cont * 0.001
        f.display_model(
            freq_cont,
            yvals_yr2,
            color="mediumblue",
            label=None,
            model_min=None,
            model_max=None,
            alpha_patch=0.2,
        )

    vip = np.array(master_pop_pd.query(f"Name=='{name}'")["VIP"])
    moss = np.array(master_pop_pd.query(f"Name=='{name}'")["MOSS"])
    f.legend(vip, moss, loc="lower center")
    f.title(name)
    f.format(xunit="GHz")

    if save_fig is True:
        f.save(f"{directory}/SEDs/{tar}_sed_testing", ext="png")
    return
