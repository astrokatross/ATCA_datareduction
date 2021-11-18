#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import pandas as pd
import numpy as np
import CFigTools.CustomFigure as CF
import gpscssmodels
import matplotlib.pyplot as plt
import ultranest
import json
import cmasher as cmr
from ultranest.plot import PredictionBand

freq_cont = np.linspace(0.01, 15, num=10000)
epochs = ["epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6", "2021-10-15"]
epoch_nms = (
    "2013",
    "2014",
    "Jan20",
    "Mar20",
    "Apr20",
    "May20",
    "July20",
    "Oct20",
)
channel = ("69", "93", "121", "145", "169")
subchans_dict = {
    "69": ["072-080", "080-088", "088-095", "095-103"],
    "93": ["103-111", "111-118", "118-126", "126-134"],
    "121": ["139-147", "147-154", "154-162", "162-170"],
    "145": ["170-177", "177-185", "185-193", "193-200"],
    "169": ["200-208", "208-216", "216-223", "223-231"],
}

xtra_fluxes = [
    "S_tgss",
    "S_mrc",
    "S_sumss",
    "S_nvss",
    "S_vlssr",
    "S20",
    "S8",
    "S5",
    "total_flux_source",  # RACS
]
exts = [
    "107",
    "115",
    "123",
    "130",
    "143",
    "150",
    "158",
    "166",
    "174",
    "181",
    "189",
    "197",
    "204",
    "212",
    "220",
    "227",
]
mwa_yr1_fluxes = ["S_076", "S_084", "S_092", "S_099"]
mwa_yr2_fluxes = []
mwa_yr1_errors = ["S_076_err", "S_084_err", "S_092_err", "S_099_err"]
mwa_yr2_errors = []
for ext in exts:
    mwa_yr1_fluxes.append(f"S_{ext}_yr1")
    mwa_yr2_fluxes.append(f"S_{ext}_yr2")
    mwa_yr1_errors.append(f"local_rms_{ext}_yr1")
    mwa_yr2_errors.append(f"local_rms_{ext}_yr2")


def read_extra_fluxes(directory, gleam_tar):
    master_pop_pd = pd.read_csv(f"{directory}/master_pop_extended.csv")
    mask = master_pop_pd["Name"] == gleam_tar
    src_pd = master_pop_pd[mask]
    fluxes_extra = np.squeeze(src_pd[xtra_fluxes].values)
    fluxes_extra[0] = fluxes_extra[0]*0.001
    fluxes_extra[2] = fluxes_extra[2]*0.001
    fluxes_extra[3] = fluxes_extra[3]*0.001
    fluxes_extra[5] = fluxes_extra[5]*0.001
    fluxes_extra[6] = fluxes_extra[6]*0.001
    fluxes_extra[7] = fluxes_extra[7]*0.001
    fluxes_extra[8] = fluxes_extra[8]*0.001
    return fluxes_extra


def read_gleam_fluxes(directory, name):
    master_pop_pd = pd.read_csv(f"{directory}/master_pop_extended.csv")
    mask = master_pop_pd["Name"] == name
    src_pd = master_pop_pd[mask]
    mwa_flux_yr1 = np.squeeze(src_pd[mwa_yr1_fluxes].values)
    err_mwa_flux_yr1 = np.squeeze(
        np.sqrt(src_pd[mwa_yr1_errors]) ** 2 + (0.02 * mwa_flux_yr1) ** 2
    )
    mwa_flux_yr2 = np.squeeze(src_pd[mwa_yr2_fluxes].values)
    err_mwa_flux_yr2 = np.squeeze(
        np.sqrt(src_pd[mwa_yr2_errors]) ** 2 + (0.02 * mwa_flux_yr2) ** 2
    )

    buffer = np.empty(4)
    buffer[:] = np.nan
    mwa_flux_yr2 = np.hstack((buffer, mwa_flux_yr2))
    err_mwa_flux_yr2 = np.hstack((buffer, err_mwa_flux_yr2))

    fluxes_extra = np.squeeze(src_pd[xtra_fluxes].values)
    fluxes_extra[0] = fluxes_extra[0]*0.001
    fluxes_extra[2] = fluxes_extra[2]*0.001
    fluxes_extra[3] = fluxes_extra[3]*0.001
    fluxes_extra[5] = fluxes_extra[5]*0.001
    fluxes_extra[6] = fluxes_extra[6]*0.001
    fluxes_extra[7] = fluxes_extra[7]*0.001
    fluxes_extra[8] = fluxes_extra[8]*0.001
    return mwa_flux_yr1, err_mwa_flux_yr1, mwa_flux_yr2, err_mwa_flux_yr2, fluxes_extra


def read_mwa_fluxes(directory, tarcomp, name, epoch):
    mwa_flux = []
    mwa_errs = []
    for i in range(len(channel)):
        subchans = subchans_dict[channel[i]]
        chan = channel[i]
        if epoch == "epoch3":
            percentage = 0.05
        else:
            percentage = 0.02
        for subchan in subchans:
            try:
                src_mwa_pd = pd.read_csv(
                    f"{directory}/{epoch}/{chan}/minimosaic/{tarcomp}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv"
                )
                mask = src_mwa_pd["Name"] == name
                src_pd = src_mwa_pd[mask]
                if src_pd.empty is False:
                    mwa_flux_chan = np.squeeze(src_pd["int_flux"].values)
                    mwa_errs_chan = np.squeeze(
                        np.sqrt(src_pd["local_rms"]) ** 2
                        + (percentage * mwa_flux_chan) ** 2
                    )
                    mwa_flux.append(mwa_flux_chan)
                    mwa_errs.append(mwa_errs_chan)
                else:
                    # print("empty pd")
                    mwa_flux.append(np.nan)
                    mwa_errs.append(np.nan)
                    pass
            except (FileNotFoundError, KeyError):
                # print(
                # f"{directory}/{epoch}/{chan}/minimosaic/{tarcomp}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv not found!"
                # )
                mwa_flux.append(np.nan)
                mwa_errs.append(np.nan)
                pass

    return np.array(mwa_flux), np.array(mwa_errs)


def read_atca_fluxes(directory, tar_dir, tar):
    atca_fluxes = [
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
    ]
    err_atca_fluxes = [
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
        [np.nan] * 17,
    ]
    epochs = ("epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6", "2021-10-15")
    for i in range(0, 7):
        epoch = epochs[i]
        try:
            atca_Lband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_L.csv")
            atca_Lband_epoch = np.array(atca_Lband_pd["# S_Lband"])
            err_atca_Lband_epoch = np.array(np.sqrt((0.05*atca_Lband_epoch)**2 + (0.0005**2)))
        except FileNotFoundError:
            # print(f"No L-Band for {epoch}")
            atca_Lband_epoch = [[np.nan] * 8]
            err_atca_Lband_epoch = [[np.nan] * 8]

        try:
            atca_Cband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_C.csv")
            atca_Cband_epoch = np.array(atca_Cband_pd["# S_Cband"])
            # print(atca_Cband_epoch)
            err_atca_Cband_epoch = np.array(np.sqrt((0.05*atca_Cband_epoch)**2 + (0.0004**2)))
            # print(err_atca_Cband_epoch)
        except FileNotFoundError:
            # print(f"No C-Band for {epoch}")
            atca_Cband_epoch = [[np.nan] * 5]
            err_atca_Cband_epoch = [[np.nan] * 5]

        try:
            atca_Xband_pd = pd.read_csv(f"{directory}/{tar_dir}/{tar}_{epoch}_X.csv")
            atca_Xband_epoch = np.array(atca_Xband_pd["# S_Xband"])
            err_atca_Xband_epoch = np.array(np.sqrt((0.05*atca_Xband_epoch)**2 + (0.0003**2)))
        except FileNotFoundError:
            # print(f"No X-Band for {epoch}")
            atca_Xband_epoch = [[np.nan] * 4]
            err_atca_Xband_epoch = [[np.nan] * 4]

        atca_epoch = np.concatenate(
            (atca_Lband_epoch, atca_Cband_epoch, atca_Xband_epoch), axis=None
        )
        err_atca_epoch = np.concatenate(
            (err_atca_Lband_epoch, err_atca_Cband_epoch, err_atca_Xband_epoch), axis=None
        )
        atca_epoch[np.where(atca_epoch < 0.0)] = np.nan
        atca_fluxes[i] = atca_epoch
        err_atca_fluxes[i] = err_atca_epoch
    return atca_fluxes, err_atca_fluxes


def create_epochcat(directory, tar, gleam_tar, epoch):
    mwa_flux, err_mwa = read_mwa_fluxes("/data/MWA", tar, gleam_tar, epochs[epoch])
    atca_flux_all, err_atca_flux_all = read_atca_fluxes(directory, tar, tar)
    atca_flux = np.array(atca_flux_all[epoch])
    err_atca_flux = err_atca_flux_all[epoch]
    src_flux = np.hstack((mwa_flux, atca_flux))
    src_errs = np.hstack((err_mwa, err_atca_flux))
    return src_flux, src_errs


def createfitflux(data_dir, gleam_target):
    target = gleam_target.strip("GLEAM ")[0:7]
    frequency = np.array(
        [
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
    )
    src_epoch1, err_src_epoch1 = create_epochcat(data_dir, target, gleam_target, 0)
    src_epoch2, err_src_epoch2 = create_epochcat(data_dir, target, gleam_target, 1)
    src_epoch3, err_src_epoch3 = create_epochcat(data_dir, target, gleam_target, 2)
    src_epoch4, err_src_epoch4 = create_epochcat(data_dir, target, gleam_target, 3)
    src_epoch5, err_src_epoch5 = create_epochcat(data_dir, target, gleam_target, 4)
    src_epoch6, err_src_epoch6 = create_epochcat(data_dir, target, gleam_target, 5)
    (
        mwa_flux_yr1,
        err_mwa_flux_yr1,
        mwa_flux_yr2,
        err_mwa_flux_yr2,
        extra_surveys,
    ) = read_gleam_fluxes("/data/MWA", gleam_target)
    if target == "J215436":
        fit_flux1 = np.hstack((mwa_flux_yr1, src_epoch4[20:37]))
        fit_flux2 = np.hstack((mwa_flux_yr2, src_epoch4[20:37]))
        fit_flux3 = np.hstack((mwa_flux_yr1, src_epoch4[20:37]))
        fit_flux4 = np.hstack((mwa_flux_yr1, src_epoch4[20:37]))
        fit_flux5 = src_epoch4
        fit_flux6 = src_epoch4
        fit_flux7 = np.hstack((src_epoch5[0:20], src_epoch4[20:37]))
        fit_flux8 = np.hstack((src_epoch6[0:20], src_epoch4[20:37]))
        err_fit_flux1 = np.hstack((err_mwa_flux_yr1, err_src_epoch4[20:37]))
        err_fit_flux2 = np.hstack((err_mwa_flux_yr2, err_src_epoch4[20:37]))
        err_fit_flux3 = np.hstack((err_mwa_flux_yr1, err_src_epoch4[20:37]))
        err_fit_flux4 = np.hstack((err_mwa_flux_yr1, err_src_epoch4[20:37]))
        err_fit_flux5 = err_src_epoch4
        err_fit_flux6 = err_src_epoch4
        err_fit_flux7 = np.hstack((err_src_epoch5[0:20], err_src_epoch4[20:37]))
        err_fit_flux8 = np.hstack((err_src_epoch6[0:20], err_src_epoch4[20:37]))
    else:
        fit_flux1 = np.hstack((mwa_flux_yr1, src_epoch4[20:37]))
        fit_flux2 = np.hstack((mwa_flux_yr2, src_epoch4[20:37]))
        fit_flux3 = np.hstack((mwa_flux_yr1, src_epoch1[20:37]))
        fit_flux4 = np.hstack((mwa_flux_yr1, src_epoch2[20:37]))
        fit_flux5 = src_epoch3
        fit_flux6 = src_epoch4
        fit_flux7 = np.hstack((src_epoch5[0:20], src_epoch4[20:37]))
        fit_flux8 = np.hstack((src_epoch6[0:20], src_epoch4[20:37]))
        err_fit_flux1 = np.hstack((err_mwa_flux_yr1, err_src_epoch4[20:37]))
        err_fit_flux2 = np.hstack((err_mwa_flux_yr2, err_src_epoch4[20:37]))
        err_fit_flux3 = np.hstack((err_mwa_flux_yr1, err_src_epoch1[20:37]))
        err_fit_flux4 = np.hstack((err_mwa_flux_yr1, err_src_epoch2[20:37]))
        err_fit_flux5 = err_src_epoch3
        err_fit_flux6 = err_src_epoch4
        err_fit_flux7 = np.hstack((err_src_epoch5[0:20], err_src_epoch4[20:37]))
        err_fit_flux8 = np.hstack((err_src_epoch6[0:20], err_src_epoch4[20:37]))
    fit_flux = np.stack(
        (
            fit_flux1,
            fit_flux2,
            fit_flux3,
            fit_flux4,
            fit_flux5,
            fit_flux6,
            fit_flux7,
            fit_flux8,
        )
    )
    err_fit_flux = np.stack(
        (
            err_fit_flux1,
            err_fit_flux2,
            err_fit_flux3,
            err_fit_flux4,
            err_fit_flux5,
            err_fit_flux6,
            err_fit_flux7,
            err_fit_flux8,
        )
    )
    fit_freq = []
    fit_flux_mask = []
    err_fit_flux_mask = []
    for i in range(8):
        mask = np.where(~np.isnan(fit_flux[i]))
        fit_flux_mask.append(fit_flux[i][mask])
        err_fit_flux_mask.append(err_fit_flux[i][mask])
        fit_freq.append(frequency[mask])
    return fit_flux_mask, err_fit_flux_mask, fit_freq


def createsrcflux(data_dir, gleam_target):
    target = gleam_target.strip("GLEAM ")[0:7]
    src_epoch1, err_src_epoch1 = create_epochcat(data_dir, target, gleam_target, 0)
    src_epoch2, err_src_epoch2 = create_epochcat(data_dir, target, gleam_target, 1)
    src_epoch3, err_src_epoch3 = create_epochcat(data_dir, target, gleam_target, 2)
    src_epoch4, err_src_epoch4 = create_epochcat(data_dir, target, gleam_target, 3)
    src_epoch5, err_src_epoch5 = create_epochcat(data_dir, target, gleam_target, 4)
    src_epoch6, err_src_epoch6 = create_epochcat(data_dir, target, gleam_target, 5)
    (
        mwa_flux_yr1,
        err_mwa_flux_yr1,
        mwa_flux_yr2,
        err_mwa_flux_yr2,
        extra_surveys,
    ) = read_gleam_fluxes("/data/MWA", gleam_target)
    atca_buffer = np.full(17, np.nan)
    mwa_flux_yr1 = np.hstack((mwa_flux_yr1, atca_buffer))
    err_mwa_flux_yr1 = np.hstack((err_mwa_flux_yr1, atca_buffer))
    mwa_flux_yr2 = np.hstack((mwa_flux_yr2, atca_buffer))
    err_mwa_flux_yr2 = np.hstack((err_mwa_flux_yr2, atca_buffer))
    # print(len(mwa_flux_yr1), len(mwa_flux_yr2), len(src_epoch1), len(src_epoch2), len(src_epoch3), len(src_epoch4), len(src_epoch5), len(src_epoch6))
    src_flux = np.stack(
        (
            mwa_flux_yr1,
            mwa_flux_yr2,
            src_epoch1,
            src_epoch2,
            src_epoch3,
            src_epoch4,
            src_epoch5,
            src_epoch6,
        )
    )
    err_src_flux = np.stack(
        (
            err_mwa_flux_yr1,
            err_mwa_flux_yr2,
            err_src_epoch1,
            err_src_epoch2,
            err_src_epoch3,
            err_src_epoch4,
            err_src_epoch5,
            err_src_epoch6,
        )
    )
    return src_flux, err_src_flux


def plot_epochsed(
    save_dir,
    freq,
    src_flux,
    err_src_flux,
    band,
    color,
    tar,
    model,
    epoch,
):
    # src_flux, err_src_flux = create_epochcat(data_dir, tar, gleam_tar, epoch)
    # plotting SED
    f = CF.sed_fig()
    f.plot_spectrum(
        freq,
        src_flux,
        err_src_flux,
        marker="o",
        marker_color=color,
        s=75,
    )
    band.line(color=color)
    band.shade(color=color, alpha=0.3)
    band.shade(color=color, q=0.01, alpha=0.1)
    f.title(f"{tar}")
    f.format(xunit="GHz")
    f.save(f"{save_dir}_sed", ext="png")
    return


def create_lnlike(freq, int_flux, err_flux, model):
    def lnlike(params):
        model_flux = model(freq, *params)
        sigma2 = err_flux ** 2
        likelihood = -0.5 * np.sum(
            (int_flux - model_flux) ** 2 / sigma2 + np.log(sigma2)
        )
        return likelihood

    return lnlike


def run_ultranest_mcmc(
    directory,
    parameters,
    freq,
    flux,
    err_flux,
    model,
    prior_transform,
    resume="resume-similar",
    run_num=1,
):
    log = ultranest.utils.make_run_dir(directory, run_num=run_num)
    # ultranest.utils.create_logger("ultranest", log_dir=directory)
    sampler = ultranest.ReactiveNestedSampler(
        param_names=parameters,
        loglike=create_lnlike(freq, flux, err_flux, model),
        transform=prior_transform,
        log_dir=log["run_dir"],
        resume=resume,
        storage_backend="hdf5",
    )
    return sampler
