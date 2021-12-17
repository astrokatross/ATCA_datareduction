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
epochs = ["2020-01", "2020-03", "2020-04", "2020-05", "2020-07", "2020-10"]
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
mwa_2013_fluxes = ["S_076", "S_084", "S_092", "S_099"]
mwa_2014_fluxes = []
mwa_2013_errors = ["S_076_err", "S_084_err", "S_092_err", "S_099_err"]
mwa_2014_errors = []
for ext in exts:
    mwa_2013_fluxes.append(f"S_{ext}_yr1")
    mwa_2014_fluxes.append(f"S_{ext}_yr2")
    mwa_2013_errors.append(f"local_rms_{ext}_yr1")
    mwa_2014_errors.append(f"local_rms_{ext}_yr2")


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


def read_mwa_fluxes(directory, tar, name, epochs):
    master_pop_pd = pd.read_csv(f"{directory}/master_pop_extended.csv")
    print(name)
    mask = master_pop_pd["Name"] == name
    src_pd = master_pop_pd[mask]

    mwa_flux = np.full([len(epochs), 20], np.nan)
    err_mwa = np.full([len(epochs), 20], np.nan)

    for j in range(len(epochs)):
        epoch = epochs[j]
        chan_flux = []
        err_chan_flux = []
        if epoch == "2013":
            chan_flux = np.squeeze(src_pd[mwa_2013_fluxes].values)
            err_chan_flux = np.squeeze(
                np.sqrt(src_pd[mwa_2013_errors]) ** 2 + (0.02 * chan_flux) ** 2
            )
        elif epoch == "2014":
            buffer = np.full(4, np.nan)
        
            chan_flux_temp = np.squeeze(src_pd[mwa_2014_fluxes].values)
            err_chan_flux_temp = np.squeeze(
                np.sqrt(src_pd[mwa_2014_errors]) ** 2 + (0.02 * chan_flux_temp) ** 2
            )
            chan_flux = np.hstack((buffer, chan_flux_temp))
            err_chan_flux = np.hstack((buffer, err_chan_flux_temp))
        else:
            for i in range(len(channel)):
                subchans = subchans_dict[channel[i]]
                chan = channel[i]
                if epoch == "2020-04":
                    percentage = 0.05
                else:
                    percentage = 0.02
                for subchan in subchans:
                    try:
                        src_mwa_pd = pd.read_csv(
                            f"{directory}/{epoch}/{chan}/minimosaic/{tar}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv"
                        )
                        mask = src_mwa_pd["Name"] == name
                        src_pd = src_mwa_pd[mask]
                        if src_pd.empty is False:
                            mwa_flux_chan = np.squeeze(src_pd["int_flux"].values)
                            mwa_errs_chan = np.squeeze(
                                np.sqrt(src_pd["local_rms"]) ** 2
                                + (percentage * mwa_flux_chan) ** 2
                            )
                            chan_flux.append(mwa_flux_chan)
                            err_chan_flux.append(mwa_errs_chan)
                        else:
                            chan_flux.append(np.nan)
                            err_chan_flux.append(np.nan)
                            pass
                    except (FileNotFoundError, KeyError):
                        chan_flux.append(np.nan)
                        err_chan_flux.append(np.nan)
                        pass

        mwa_flux[j] = chan_flux
        err_mwa[j] = err_chan_flux
    return mwa_flux, err_mwa


def read_atca_fluxes(directory, tar, epochs):
    atca_fluxes = np.full([len(epochs), 17], np.nan)
    # atca_fluxes[:] = np.nan

    err_atca_fluxes = np.full([len(epochs), 17], np.nan)
    for i in range(len(epochs)):
        epoch = epochs[i]
        try:
            atca_Lband_pd = pd.read_csv(f"{directory}/{tar}/{epoch}_{tar}_L.csv")
            atca_Lband_epoch = np.array(atca_Lband_pd["# S_Lband"])
            err_atca_Lband_epoch = np.array(np.sqrt((0.05*atca_Lband_epoch)**2 + (0.0005**2)))

        except FileNotFoundError:
            # print(f"No L-Band for {epoch}")
            atca_Lband_epoch = [[np.nan] * 8]
            err_atca_Lband_epoch = [[np.nan] * 8]

        try:
            atca_Cband_pd = pd.read_csv(f"{directory}/{tar}/{epoch}_{tar}_C.csv")
            atca_Cband_epoch = np.array(atca_Cband_pd["# S_Cband"])
            # print(atca_Cband_epoch)
            err_atca_Cband_epoch = np.array(np.sqrt((0.05*atca_Cband_epoch)**2 + (0.0004**2)))
            # print(err_atca_Cband_epoch)
        except FileNotFoundError:
            # print(f"No C-Band for {epoch}")
            atca_Cband_epoch = [[np.nan] * 5]
            err_atca_Cband_epoch = [[np.nan] * 5]

        try:
            atca_Xband_pd = pd.read_csv(f"{directory}/{tar}/{epoch}_{tar}_X.csv")
            atca_Xband_epoch = np.array(atca_Xband_pd["# S_Xband"])
            err_atca_Xband_epoch = np.array(np.sqrt((0.05*atca_Xband_epoch)**2 + (0.0003**2)))
        except FileNotFoundError:
            # print(f"No X-Band for {epoch}")
            atca_Xband_epoch = [[np.nan] * 4]
            err_atca_Xband_epoch = [[np.nan] * 4]

        atca_epoch = np.concatenate(
            (atca_Lband_epoch, atca_Cband_epoch, atca_Xband_epoch),axis=None)
        err_atca_epoch = np.concatenate(
            (err_atca_Lband_epoch, err_atca_Cband_epoch, err_atca_Xband_epoch),axis=None)
        # print(atca_epoch)
        # print(i)
        # print(atca_fluxes[i])
        # print(i)
        # atca_epoch[np.where(atca_epoch < 0.0)] = np.nan
        atca_fluxes[i] = atca_epoch
        err_atca_fluxes[i] = err_atca_epoch
        # print(i)
        # print(atca_epoch)
        # print(atca_fluxes[i])
    return atca_fluxes.squeeze(), err_atca_fluxes.squeeze()


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

    mwa_fluxes, err_mwa_fluxes = read_mwa_fluxes("/data/MWA", target, gleam_target, ["2013", "2014", "2020-04", "2020-05", "2020-07", "2020-10"])
    try:
        atca_fluxes, err_atca_fluxes = read_atca_fluxes(data_dir, target, ["2020"])
        atca_fluxes_sum = np.sum(atca_fluxes)
        # print(atca_fluxes_sum)
        if np.isnan(atca_fluxes_sum) == True:
            print("No 2020 csv, finding average")
            raise ValueError
    except ValueError:
        atca_fluxes_all, err_atca_fluxes_all = read_atca_fluxes(data_dir, target, ["2020-01", "2020-03", "2020-04", "2020-05"])
        atca_fluxes = np.nanmean(atca_fluxes_all, axis=0)
        err_atca_fluxes = np.nanmean(err_atca_fluxes_all, axis=0)
    # print(atca_fluxes)
    src_flux = []
    err_src_flux = []
    for i in range(len(mwa_fluxes)):
        src_flux.append(np.hstack((mwa_fluxes[i], atca_fluxes)))
        err_src_flux.append(np.hstack((err_mwa_fluxes[i], err_atca_fluxes)))
    fit_freq = []
    fit_flux_mask = []
    err_fit_flux_mask = []
    for i in range(len(src_flux)):
        mask = np.where(~np.isnan(src_flux[i]))
        fit_flux_mask.append(src_flux[i][mask])
        err_fit_flux_mask.append(err_src_flux[i][mask])
        fit_freq.append(frequency[mask])
    return fit_flux_mask, err_fit_flux_mask, fit_freq


def createsrcflux(data_dir, gleam_target, epochs=["2013", "2014", "2020-01", "2020-03", "2020-04", "2020-05", "2020-07", "2020-10"]):
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

    mwa_fluxes, err_mwa_fluxes = read_mwa_fluxes("/data/MWA", target, gleam_target, epochs)
    atca_fluxes, err_atca_fluxes = read_atca_fluxes(data_dir, target, epochs)
    try:
        atca_fluxes_2020, err_atca_fluxes_2020 = read_atca_fluxes(data_dir, target, ["2020"])
        atca_fluxes_sum = np.sum(atca_fluxes_2020)
        # print(atca_fluxes_sum)
        if np.isnan(atca_fluxes_sum) == True:
            print("No 2020 csv, finding average")
            raise ValueError
    except ValueError:
        atca_fluxes_all, err_atca_fluxes_all = read_atca_fluxes(data_dir, target, ["2020-01", "2020-03", "2020-04", "2020-05"])
        atca_fluxes_2020 = np.nanmean(atca_fluxes_all, axis=0)
        err_atca_fluxes_2020 = np.nanmean(err_atca_fluxes_all, axis=0)
    # atca_fluxes = np.concatenate((atca_fluxes, atca_fluxes_2020), axis=0)
    # err_atca_fluxes = np.concatenate((err_atca_fluxes, err_atca_fluxes_2020), axis=0)
    # print(atca_fluxes)
    src_flux = []
    err_src_flux = []
    for i in range(len(epochs)):
        src_flux.append(np.concatenate((mwa_fluxes[i], atca_fluxes[i]), axis=None))
        err_src_flux.append(np.concatenate((err_mwa_fluxes[i], err_atca_fluxes[i]),axis=None))
    mwa_buffer = np.empty((20,))
    mwa_buffer[:] = np.nan
    src_flux.append(np.concatenate((mwa_buffer, atca_fluxes_2020), axis=None))
    err_src_flux.append(np.concatenate((mwa_buffer, err_atca_fluxes_2020), axis=None))
    # print(len(src_flux))
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
