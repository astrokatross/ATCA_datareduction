#!/usr/bin/python3
# This script is to plot and analyse the ATCA data reduction
# By K.Ross 20/1/21

import pandas as pd
import matplotlib.pyplot as plt
import CFigTools.CustomFigure as CF
import numpy as np
import gpscssmodels
import cmasher as cmr
import scipy.optimize as opt
import emcee
import corner

plt.rcParams["axes.grid"] = False
plt.rcParams["font.family"] = "serif"

freq_cont = np.linspace(0.01, 15, num=10000)
epochs = ["epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6"]
epoch_nms = ("2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20")
channel = ("69", "93", "121", "145", "169")
subchans_dict = {
    "69": ["072-080", "080-088", "088-095", "095-103"],
    "93": ["103-111", "111-118", "118-126", "126-134"],
    "121": ["139-147", "147-154", "154-162", "162-170"],
    "145": ["170-177", "177-185", "185-193", "193-200"],
    "169": ["200-208", "208-216", "216-223", "223-231"],
}


def read_mwa_fluxes(directory, tarcomp, name, epoch):
    mwa_flux = []
    mwa_errs = []
    for i in range(len(channel)):
        subchans = subchans_dict[channel[i]]
        chan = channel[i]
        for subchan in subchans:
            try:
                src_mwa_pd = pd.read_csv(
                    f"{directory}/{epoch}/{chan}/minimosaic/{tarcomp}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv"
                )
                mask = src_mwa_pd["Name"] == name
                src_pd = src_mwa_pd[mask]
                mwa_flux_chan = np.squeeze(src_pd["int_flux"].values)
                mwa_errs_chan = np.squeeze(
                    np.sqrt(src_pd["local_rms"]) ** 2 + (0.02 * mwa_flux_chan) ** 2
                )
                mwa_flux.append(mwa_flux_chan)
                mwa_errs.append(mwa_errs_chan)
            except (FileNotFoundError, KeyError):
                # print(
                # f"{directory}/{epoch}/{chan}/minimosaic/{tarcomp}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv not found!"
                # )
                mwa_flux.append(np.nan)
                mwa_errs.append(np.nan)
            except:
                # print(f"{name} {subchan} not found in catalogue! Maybe too faint")
                mwa_flux.append(np.nan)
                mwa_errs.append(np.nan)

    return np.array(mwa_flux), np.array(mwa_errs)


def read_atca_fluxes(directory, tar_dir, tar):
    atca_fluxes = [
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
    ]
    epochs = ("epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6")
    for i in range(0, 6):
        epoch = epochs[i]
        try:
            atca_Lband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_L.csv")
            atca_Lband_epoch = np.array(atca_Lband_pd["# S_Lband"])
        except FileNotFoundError:
            # print(f"No L-Band for {epoch}")
            atca_Lband_epoch = [[np.nan] * 8]

        try:
            atca_Cband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_C.csv")
            atca_Cband_epoch = np.array(atca_Cband_pd["# S_Cband"])
        except FileNotFoundError:
            # print(f"No C-Band for {epoch}")
            atca_Cband_epoch = [[np.nan] * 5]

        try:
            atca_Xband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_X.csv")
            atca_Xband_epoch = np.array(atca_Xband_pd["# S_Xband"])
        except FileNotFoundError:
            # print(f"No X-Band for {epoch}")
            atca_Xband_epoch = [[np.nan] * 4]

        atca_epoch = np.concatenate(
            (atca_Lband_epoch, atca_Cband_epoch, atca_Xband_epoch), axis=None
        )
        atca_epoch[np.where(atca_epoch < 0.0)] = np.nan
        atca_fluxes[i] = atca_epoch
    return atca_fluxes


def create_epochcat(directory, gleam_tar, tar, epoch):
    mwa_flux, err_mwa = read_mwa_fluxes("/data/MWA", tar, gleam_tar, epochs[epoch])
    atca_flux_all = read_atca_fluxes(directory, tar, tar)
    atca_flux = atca_flux_all[epoch]
    src_flux = np.hstack((mwa_flux, atca_flux))
    src_errs = np.hstack((err_mwa, atca_flux * 0.05))
    return src_flux, src_errs


def fit_model(freq, flux, err_flux, model):
    poptsingssa, pcovsingssa = opt.curve_fit(
        model, freq, flux, sigma=err_flux, maxfev=10000
    )
    return poptsingssa, pcovsingssa


def log_likelihood(params, freq, flux, err_flux, model):
    model_flux = model(freq, *params)
    sigma2 = err_flux ** 2
    return -0.5 * np.sum((flux - model_flux) ** 2 / sigma2 + np.log(sigma2))


def log_prior(theta):
    S_norm, beta, peak_freq = theta
    if 0.0 < S_norm < 20 and 0.0 < beta < 10.0 and 0.0 < peak_freq < 20:
        return 0.0
    return -np.inf


def log_probability(params, freq, flux, err_flux, model):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, freq, flux, err_flux, model)


def run_mcmc(nwalkers, ndim, p0, freq, flux, err_flux, chosen_model, niters=10000):
    # Setting up and running the MCMC
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        log_likelihood,
        args=(freq, flux, err_flux, chosen_model),
    )
    sampler.run_mcmc(p0, niters, progress=True)
    return sampler


def plot_mcmcqualitycheck(directory, figname, sampler, p0, labels):
    # Plotting the different walkers and how they sample the space/converge
    fig, axes = plt.subplots(len(labels), figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(len(labels)):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    plt.title(f"{figname} Parameter Chains")
    plt.savefig(f"{directory}/{figname}_paramchains.png", overwrite=True)
    plt.clf()
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    # Plotting the corner plot of the different variables
    corner.corner(flat_samples, labels=labels, truths=p0)
    plt.title()
    plt.savefig(f"{directory}/{figname}_corner.png", overwrite=True)
    plt.clf()
    return


def plot_sed(directory, freq_cont, freq, gleam_tar, tar, colors):
    # mwa_flux, model_vals, xtra_values, values = read_gleam_fluxes(
    #     "/data/MWA", gleam_tar
    # )
    src_epoch1, err_src_epoch1 = create_epochcat(directory, gleam_tar, tar, 0)
    src_epoch2, err_src_epoch2 = create_epochcat(directory, gleam_tar, tar, 1)
    src_epoch3, err_src_epoch3 = create_epochcat(directory, gleam_tar, tar, 2)
    src_epoch4, err_src_epoch4 = create_epochcat(directory, gleam_tar, tar, 3)
    src_epoch5, err_src_epoch5 = create_epochcat(directory, gleam_tar, tar, 4)
    src_epoch6, err_src_epoch6 = create_epochcat(directory, gleam_tar, tar, 5)
    # plotting SED
    f = CF.sed_fig()
    # f.plot_spectrum(
    #     freq[4:20],
    #     mwa_flux[0],
    #     mwa_flux[2],
    #     marker="o",
    #     label=epoch_nms[0],
    #     marker_color=colors[0],
    #     s=75,
    # )
    # f.plot_spectrum(
    #     freq[4:20],
    #     mwa_flux[1],
    #     mwa_flux[3],
    #     marker="o",
    #     label=epoch_nms[1],
    #     marker_color=colors[1],
    #     s=75,
    # )

    f.plot_spectrum(
        freq,
        src_epoch1,
        err_src_epoch1,
        marker="o",
        label=epoch_nms[2],
        marker_color=colors[2],
        s=75,
    )
    f.plot_spectrum(
        freq,
        src_epoch2,
        err_src_epoch2,
        marker="o",
        label=epoch_nms[3],
        marker_color=colors[3],
        s=75,
    )
    f.plot_spectrum(
        freq,
        src_epoch3,
        err_src_epoch3,
        marker="o",
        label=epoch_nms[4],
        marker_color=colors[4],
        s=75,
    )
    f.plot_spectrum(
        freq,
        src_epoch4,
        err_src_epoch4,
        marker="o",
        label=epoch_nms[5],
        marker_color=colors[5],
        s=75,
    )
    f.plot_spectrum(
        freq,
        src_epoch5,
        err_src_epoch5,
        marker="o",
        label=epoch_nms[6],
        marker_color=colors[6],
        s=75,
    )
    f.plot_spectrum(
        freq,
        src_epoch6,
        err_src_epoch6,
        marker="o",
        label=epoch_nms[7],
        marker_color=colors[7],
        s=75,
    )
    # f.plt_mcmcfits(sampler, chosen_model, freq_cont, color=colors[5])
    f.legend(0.0, 0.0, loc="lower center")
    f.title(gleam_tar)
    f.format(xunit="GHz")
    f.save(f"{directory}/SEDs/{tar}_sed", ext="png")
    return


def plot_epochsed(
    directory,
    epoch,
    freq_cont,
    freq,
    gleam_tar,
    tar,
    colors,
    sampler,
    chosen_model,
    model_nm,
):
    epoch_nms = ("2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20")
    epoch_nm = epoch_nms[epoch]
    src_flux, err_src_flux = create_epochcat(directory, gleam_tar, tar, 0)
    # plotting SED
    f = CF.sed_fig()
    f.plot_spectrum(
        freq,
        src_flux,
        err_src_flux,
        marker="o",
        label=epoch_nms[epoch],
        marker_color=colors[epoch],
        s=75,
    )
    f.plt_mcmcfits(sampler, chosen_model, freq_cont, color=colors[epoch])
    f.legend(loc="lower center")
    f.title(f"{gleam_tar} {epoch_nm} {model_nm}")
    f.format(xunit="GHz")
    f.save(f"{directory}/SEDs/{tar}_{epoch_nm}_{model_nm}_sed", ext="png")
    return


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
    20,
    8.6,
    4.8,
    0.8875,
]  # [gleamx4, TGSS, MRC, SUMSS, NVSS, VLSSr, AT20G, AT20G8, AT20G5]

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
    "S20",
    "S8",
    "S5",
    "total_flux_source",
]

NUM_COLORS = 8

# colors = []
# cm = pylab.get_cmap('gnuplot2')
# cm = cmr.bubblegum
colors = cmr.take_cmap_colors(
    "cmr.bubblegum", NUM_COLORS, cmap_range=(0.15, 0.85), return_fmt="hex"
)
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
    "J001513": ["GLEAM J001513-472706", "2327-459", "GLEAM J002130-465755", "J002130"],
    "J015445": ["GLEAM J015445-232950", "0237-233", "GLEAM J015010-225841", "J015010"],
    "J020507": ["GLEAM J020507-110922", "0238-084", "GLEAM J020510-113515", "J020510"],
    "J021246": ["GLEAM J021246-305454", "0237-233", "GLEAM J021622-300933", "J021622"],
    "J022744": ["GLEAM J022744-062106", "0238-084", "GLEAM J023253-064552", "J023253"],
    "J024838": ["GLEAM J024838-321336", "0237-233", "GLEAM J025040-313523", "J025040"],
    "J032213": ["GLEAM J032213-462646", "0355-483", "GLEAM J033336-443859", "J033336"],
    "J032836": ["GLEAM J032836-202138", "0310-150", "GLEAM J032843-205825", "J032843"],
    "J033023": ["GLEAM J033023-074052", "0310-150", "GLEAM J033156-071220", "J033156"],
    "J042502": ["GLEAM J042502-245129", "0445-221", "GLEAM J042256-250339", "J042256"],
    "J044033": ["GLEAM J044033-422918", "0355-483", "GLEAM J044149-430658", "J044149"],
    "J044737": ["GLEAM J044737-220335", "0445-221", "GLEAM J044313-220224", "J044313"],
    "J052824": ["GLEAM J052824-331104", "0528-250", "GLEAM J053114-334532", "J053114"],
    "J223933": ["GLEAM J223933-451414", "2327-459", "GLEAM J224320-454100", "J224320"],
    "J224408": ["GLEAM J224408-202719", "2240-260", "GLEAM J224252-200050", "J224252"],
    "J215436": ["GLEAM J215436-410853"],
    "j032237": ["GLEAM J032237-482010"],
}
subchans_dict = {
    "69": ["072-080", "080-088", "088-095", "095-103"],
    "93": ["103-111", "111-118", "118-126", "126-134"],
    "121": ["139-147", "147-154", "154-162", "162-170"],
    "145": ["170-177", "177-185", "185-193", "193-200"],
    "169": ["200-208", "208-216", "216-223", "223-231"],
}
epochs = ("epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6")
extensions = [
    "072_080",
    "080_088",
    "088_095",
    "095_103",
    "103_111",
    "111_118",
    "118_126",
    "126_134",
    "139_147",
    "147_154",
    "154_162",
    "162_170",
    "170_177",
    "177_185",
    "185_193",
    "193_200",
    "200_208",
    "208_216",
    "216_223",
    "223_231",
]
gleamx_fluxes = []
gleamx_errors = []
for subchans in extensions:
    gleamx_fluxes.append(f"int_flux_N_{subchans}MHz")
    gleamx_errors.append(f"local_rms_N_{subchans}MHz")


def read_gleam_fluxes(directory, name):
    master_pop_pd = pd.read_csv(f"{directory}/master_pop_extended.csv")
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


def read_gleamx_fluxes(directory, name):
    gleamx_pop_pd = pd.read_csv(
        f"{directory}/XG_D-27_20201001_joined_rescaled_cat_015445_xmatch.csv"
    )
    mask = gleamx_pop_pd["Name"] == name
    sub_csv = gleamx_pop_pd[mask]
    try:
        gleamx_flux = np.squeeze(sub_csv[gleamx_fluxes].values)
        gleamx_errs = np.squeeze(
            np.sqrt(sub_csv[gleamx_errors] ** 2 + (0.02 * gleamx_flux) ** 2)
        )
    #         print(gleamx_errs)
    except:
        gleamx_flux = np.array([np.nan] * 20)
        gleamx_errs = np.array([np.nan] * 20)

    return gleamx_flux, gleamx_errs


# def read_mwa_fluxes(directory, tarcomp, name, epoch):
#     mwa_flux = []
#     mwa_errs = []
#     for i in range(len(channel)):
#         subchans = subchans_dict[channel[i]]
#         chan = channel[i]
#         for subchan in subchans:
#             try:
#                 src_mwa_pd = pd.read_csv(
#                     f"{directory}/{epoch}/{chan}/minimosaic/{tarcomp}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv"
#                 )
#                 mask = src_mwa_pd["Name"] == name
#                 src_pd = src_mwa_pd[mask]
#                 mwa_flux_chan = np.squeeze(src_pd["int_flux"].values)
#                 mwa_errs_chan = np.squeeze(
#                     np.sqrt(src_pd["local_rms"]) ** 2 + (0.02 * mwa_flux_chan) ** 2
#                 )
#                 mwa_flux.append(mwa_flux_chan)
#                 mwa_errs.append(mwa_errs_chan)
#             except (FileNotFoundError, KeyError):
#                 # print(
#                 # f"{directory}/{epoch}/{chan}/minimosaic/{tarcomp}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv not found!"
#                 # )
#                 mwa_flux.append(np.nan)
#                 mwa_errs.append(np.nan)
#             except:
#                 # print(f"{name} {subchan} not found in catalogue! Maybe too faint")
#                 mwa_flux.append(np.nan)
#                 mwa_errs.append(np.nan)

#     return np.array(mwa_flux), np.array(mwa_errs)


# def read_atca_fluxes(directory, tar_dir, tar):
#     atca_fluxes = [
#         [np.nan] * 17,
#         [np.nan] * 17,
#         [np.nan] * 17,
#         [np.nan] * 17,
#         [np.nan] * 17,
#         [np.nan] * 17,
#     ]
#     epochs = ("epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6")
#     for i in range(0, 6):
#         epoch = epochs[i]
#         try:
#             atca_Lband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_L.csv")
#             atca_Lband_epoch = np.array(atca_Lband_pd["# S_Lband"])
#         except FileNotFoundError:
#             # print(f"No L-Band for {epoch}")
#             atca_Lband_epoch = [[np.nan] * 8]

#         try:
#             atca_Cband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_C.csv")
#             atca_Cband_epoch = np.array(atca_Cband_pd["# S_Cband"])
#         except FileNotFoundError:
#             # print(f"No C-Band for {epoch}")
#             atca_Cband_epoch = [[np.nan] * 5]

#         try:
#             atca_Xband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_X.csv")
#             atca_Xband_epoch = np.array(atca_Xband_pd["# S_Xband"])
#         except FileNotFoundError:
#             # print(f"No X-Band for {epoch}")
#             atca_Xband_epoch = [[np.nan] * 4]

#         atca_epoch = np.concatenate(
#             (atca_Lband_epoch, atca_Cband_epoch, atca_Xband_epoch), axis=None
#         )
#         atca_epoch[np.where(atca_epoch < 0.0)] = np.nan
#         atca_fluxes[i] = atca_epoch
#     return atca_fluxes


# Originally in Tingay, de Kool 2003
def singSSA(freq, S_norm, beta, peak_freq):  # Single SSA model
    return (
        S_norm
        * ((freq / peak_freq) ** (-(beta - 1) / 2))
        * (1 - np.exp(-((freq / peak_freq) ** (-(beta + 4) / 2))))
        / ((freq / peak_freq) ** (-(beta + 4) / 2))
    )


def fit_singSSA(freq, flux, err_flux):
    poptsingssa, pcovsingssa = opt.curve_fit(
        singSSA, freq, flux, sigma=err_flux, maxfev=10000
    )
    return poptsingssa, pcovsingssa


def fit_quadSSA(freq, flux, err_flux):
    poptquadssa, pcovquadssa = opt.curve_fit(
        gpscssmodels.singhomobremss, freq, flux, sigma=err_flux, maxfev=10000
    )
    return poptquadssa, pcovquadssa


# def create_epochcat(directory, gleam_tar, tar, epoch):
#     mwa_flux, err_mwa = read_mwa_fluxes("/data/MWA", tar, gleam_tar, epochs[epoch])
#     atca_flux_all = read_atca_fluxes(directory, tar, tar)
#     atca_flux = atca_flux_all[epoch]
#     src_flux = np.hstack((mwa_flux, atca_flux))
#     src_errs = np.hstack((err_mwa, atca_flux * 0.05))
#     return src_flux, src_errs


def plt_sed(
    directory,
    tar,
    MWA=True,
    extra_surveys=True,
    save_fig=True,
):
    colours = ("C3", "mediumseagreen", "C9", "C4", "C1", "C8")
    epoch_nms = ("Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20")
    atca_fluxes = read_atca_fluxes(directory, tar, tar)
    tarnm = source_dict[tar][0]
    if tar == "j032237":
        mwa_fluxes_e3, mwa_errors_e3 = read_mwa_fluxes(
            "/data/MWA", "J032213", tarnm, "epoch3"
        )
        mwa_fluxes_e4, mwa_errors_e4 = read_mwa_fluxes(
            "/data/MWA", "J032213", tarnm, "epoch4"
        )
        mwa_fluxes_e5, mwa_errors_e5 = read_mwa_fluxes(
            "/data/MWA", "J032213", tarnm, "epoch5"
        )
        mwa_fluxes_e6, mwa_errors_e6 = read_mwa_fluxes(
            "/data/MWA", "J032213", tarnm, "epoch6"
        )
    else:
        mwa_fluxes_e3, mwa_errors_e3 = read_mwa_fluxes(
            "/data/MWA", tar, tarnm, "epoch3"
        )
        mwa_fluxes_e4, mwa_errors_e4 = read_mwa_fluxes(
            "/data/MWA", tar, tarnm, "epoch4"
        )
        mwa_fluxes_e5, mwa_errors_e5 = read_mwa_fluxes(
            "/data/MWA", tar, tarnm, "epoch5"
        )
        mwa_fluxes_e6, mwa_errors_e6 = read_mwa_fluxes(
            "/data/MWA", tar, tarnm, "epoch6"
        )
        gleamxfluxes, gleamxerrors = read_gleamx_fluxes("/data/MWA", tarnm)
        print(gleamxfluxes)
    f = CF.sed_fig()
    if MWA is True:
        mwa_flux, model_vals, xtra_values, values = read_gleam_fluxes(
            "/data/MWA", tarnm
        )
        f.plot_spectrum(
            freq_mwa[4:20],
            mwa_flux[0],
            mwa_flux[2],
            marker="o",
            label="2013",
            marker_color="C6",
            alpha=0.5,
        )
        f.plot_spectrum(
            freq_mwa[4:20],
            mwa_flux[1],
            mwa_flux[3],
            marker="o",
            label="2014",
            marker_color="mediumblue",
            alpha=0.5,
        )
        # f.plot_residuals(
        #     freq_mwa[4:20],
        #     mwa_flux[0],
        #     mwa_flux[1],
        #     mwa_flux[2],
        #     mwa_flux[3],
        # )
        f.display_model(
            freq_cont,
            model_vals[0],
            color="C6",
            label=None,
            model_min=None,
            model_max=None,
            alpha_patch=0.2,
            alpha=0.5,
        )
        f.display_model(
            freq_cont,
            model_vals[1],
            color="mediumblue",
            label=None,
            model_min=None,
            model_max=None,
            alpha_patch=0.2,
            alpha=0.5,
        )
    for i in range(0, 6):
        f.plot_spectrum(
            freq_atca,
            atca_fluxes[i],
            atca_fluxes[i] * 0.03,
            marker="o",
            marker_color=colours[i],
            label=epoch_nms[i],
        )
    for i in [0, 1, 3, 4]:
        f.plot_residuals(
            freq_atca,
            atca_fluxes[2],
            atca_fluxes[i],
            atca_fluxes[2] * 0.03,
            atca_fluxes[i] * 0.03,
            color=colours[i],
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
            alpha=0.5,
        )
        if fluxes_extra[8] != np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[8],
                fluxes_extra[8],
                marker="X",
                # label="VLSSr",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[4] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[4],
                fluxes_extra[4] * 0.001,
                marker="s",
                # label="TGSS",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[5] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[5],
                fluxes_extra[5],
                marker="*",
                # label="MRC",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[6] == np.nan:
            print("No SUMSS")
            pass
        else:
            f.plot_point(
                freq_xtra[6],
                fluxes_extra[6] * 0.001,
                marker="p",
                # label="SUMSS",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[7] == np.nan:
            pass
        else:
            f.plot_point(
                freq_xtra[7],
                fluxes_extra[7] * 0.001,
                marker="D",
                # label="NVSS",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[9] == np.nan or fluxes_extra[9] == 0:
            pass
        else:
            f.plot_point(
                freq_xtra[9],
                fluxes_extra[9] * 0.001,
                marker="1",
                # label="NVSS",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[10] == np.nan or fluxes_extra[10] == 0:
            pass
        else:
            f.plot_point(
                freq_xtra[10],
                fluxes_extra[10] * 0.001,
                marker="2",
                # label="NVSS",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[11] == np.nan or fluxes_extra[11] == 0:
            pass
        else:
            f.plot_point(
                freq_xtra[11],
                fluxes_extra[11] * 0.001,
                marker="3",
                # label="NVSS",
                marker_color="k",
                alpha=1,
            )
        if fluxes_extra[12] == np.nan or fluxes_extra[12] == 0:
            pass
        else:
            f.plot_point(
                freq_xtra[12],
                fluxes_extra[12] * 0.001,
                marker="4",
                # label="NVSS",
                marker_color="k",
                alpha=1,
            )
    fluxes_extra = xtra_values[0]
    mwa_2013 = np.concatenate((fluxes_extra[0:4], mwa_flux[0]))
    gleam_err = xtra_values[1]
    mwa_2013_err = np.concatenate((gleam_err, mwa_flux[2]))
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e3,
        mwa_errors_e3,
        marker="o",
        s=40,
        marker_color=colours[2],
    )
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e4,
        mwa_errors_e4,
        marker="o",
        s=40,
        marker_color=colours[3],
    )
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e5,
        mwa_errors_e5,
        marker="o",
        s=40,
        # label=epoch_nms[4],
        marker_color=colours[4],
    )
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e6,
        mwa_errors_e6,
        marker="o",
        s=40,
        # label=epoch_nms[5],
        marker_color=colours[5],
    )
    f.plot_spectrum(
        freq_mwa,
        gleamxfluxes,
        gleamxerrors,
        marker="o",
        s=40,
        # label=epoch_nms[5],
        marker_color="k",
    )
    f.plot_residuals(
        freq_mwa[4:20],
        mwa_flux[0],
        mwa_flux[1],
        mwa_flux[2],
        mwa_flux[3],
        alpha=0.5,
        color="mediumblue",
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e3,
        mwa_2013_err,
        mwa_errors_e3,
        color=colours[2],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e4,
        mwa_2013_err,
        mwa_errors_e4,
        color=colours[3],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e5,
        mwa_2013_err,
        mwa_errors_e5,
        color=colours[4],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e6,
        mwa_2013_err,
        mwa_errors_e6,
        color=colours[5],
    )

    # if MWA == True:
    f.legend(values[1], values[2], loc="lower center")
    f.title(values[0])
    f.format(xunit="GHz")

    if save_fig is True:
        f.save(f"{directory}/SEDs/{tar}_sed", ext="png")
    return


def plt_secondary(directory, tar, save_fig=True):
    name = source_dict[tar][1]
    atca_fluxes = read_atca_fluxes(directory, tar, name)
    colours = ("C3", "mediumseagreen", "C9", "C4", "C1")
    epoch_nms = ("Jan20", "Mar20", "Apr20", "May20", "")
    f = CF.sed_fig()
    for i in range(0, 4):
        f.plot_spectrum(
            freq_atca,
            atca_fluxes[i],
            atca_fluxes[i] * 0.03,
            marker="o",
            marker_color=colours[i],
            label=epoch_nms[i],
        )
        f.plot_residuals(
            freq_atca,
            atca_fluxes[2],
            atca_fluxes[i],
            atca_fluxes[2] * 0.03,
            atca_fluxes[i] * 0.03,
            color=colours[i],
        )
    f.format(xunit="GHz")
    f.title(name)
    f.legend(0, 0, loc="lower left")
    if save_fig is True:
        f.save(f"{directory}/SEDs/secondary/{tar}_sec_sed", ext="png")
    return


def plt_nearby(directory, tar, MWA=True, extra_surveys=True, save_fig=True):
    colours = ("C3", "mediumseagreen", "C9", "C4", "C1", "C8")
    epoch_nms = ("Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20")
    name = source_dict[tar][3]
    secnm = source_dict[tar][2]
    mwa_fluxes_e3, mwa_errors_e3 = read_mwa_fluxes("/data/MWA", tar, secnm, "epoch3")
    mwa_fluxes_e4, mwa_errors_e4 = read_mwa_fluxes("/data/MWA", tar, secnm, "epoch4")
    mwa_fluxes_e5, mwa_errors_e5 = read_mwa_fluxes("/data/MWA", tar, secnm, "epoch5")
    mwa_fluxes_e6, mwa_errors_e6 = read_mwa_fluxes("/data/MWA", tar, secnm, "epoch6")
    mwa_flux, model_vals, xtra_values, values = read_gleam_fluxes("/data/MWA", secnm)
    fluxes_extra = xtra_values[0]
    mwa_2013 = np.concatenate((fluxes_extra[0:4], mwa_flux[0]))
    gleam_err = xtra_values[1]
    mwa_2013_err = np.concatenate((gleam_err, mwa_flux[2]))
    f = CF.sed_fig()
    f.plot_spectrum(
        freq_mwa[4:20],
        mwa_flux[0],
        mwa_flux[2],
        marker="o",
        label="2013",
        marker_color="C6",
    )
    f.plot_spectrum(
        freq_mwa[4:20],
        mwa_flux[1],
        mwa_flux[3],
        marker="o",
        label="2014",
        marker_color="mediumblue",
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
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e3,
        mwa_errors_e3,
        marker="o",
        label=epoch_nms[2],
        marker_color=colours[2],
    )
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e4,
        mwa_errors_e4,
        marker="o",
        label=epoch_nms[3],
        marker_color=colours[3],
    )
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e5,
        mwa_errors_e5,
        marker="o",
        label=epoch_nms[4],
        marker_color=colours[4],
    )
    f.plot_spectrum(
        freq_mwa,
        mwa_fluxes_e6,
        mwa_errors_e6,
        marker="o",
        label=epoch_nms[5],
        marker_color=colours[5],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e3,
        mwa_2013_err,
        mwa_errors_e3,
        color=colours[2],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e4,
        mwa_2013_err,
        mwa_errors_e4,
        color=colours[3],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e5,
        mwa_2013_err,
        mwa_errors_e5,
        color=colours[4],
    )
    f.plot_residuals(
        freq_mwa,
        mwa_2013,
        mwa_fluxes_e6,
        mwa_2013_err,
        mwa_errors_e6,
        color=colours[5],
    )
    f.format(xunit="GHz")
    f.title(secnm)
    f.legend(0, 0)
    if save_fig is True:
        f.save(f"{directory}/SEDs/nearby/{name}_sed", ext="png")
    return
