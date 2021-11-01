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
    # "Oct21",
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
                    print("empty pd")
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
    epochs = ("epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6", "2021-10-15")
    for i in range(0, 7):
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


def create_epochcat(directory, tar, gleam_tar, epoch):
    mwa_flux, err_mwa = read_mwa_fluxes("/data/MWA", tar, gleam_tar, epochs[epoch])
    atca_flux_all = read_atca_fluxes(directory, tar, tar)
    atca_flux = atca_flux_all[epoch]
    src_flux = np.hstack((mwa_flux, atca_flux))
    src_errs = np.hstack((err_mwa, atca_flux * 0.05))
    return src_flux, src_errs


def plot_sed(save_dir, data_dir, freq, gleam_tar, tar, colors, model_info, model):
    (
        mwa_flux_yr1,
        err_mwa_flux_yr1,
        mwa_flux_yr2,
        err_mwa_flux_yr2,
        fluxes_extra,
    ) = read_gleam_fluxes("/data/MWA", gleam_tar)
    atca_buffer = np.empty(17)
    atca_buffer[:] = np.nan
    mwa_flux_yr1 = np.hstack((mwa_flux_yr1, atca_buffer))
    err_mwa_flux_yr1 = np.hstack((err_mwa_flux_yr1, atca_buffer))
    mwa_flux_yr2 = np.hstack((mwa_flux_yr2, atca_buffer))
    err_mwa_flux_yr2 = np.hstack((err_mwa_flux_yr2, atca_buffer))
    src_epoch1, err_src_epoch1 = create_epochcat(data_dir, tar, gleam_tar, 0)
    src_epoch2, err_src_epoch2 = create_epochcat(data_dir, tar, gleam_tar, 1)
    src_epoch3, err_src_epoch3 = create_epochcat(data_dir, tar, gleam_tar, 2)
    src_epoch4, err_src_epoch4 = create_epochcat(data_dir, tar, gleam_tar, 3)
    src_epoch5, err_src_epoch5 = create_epochcat(data_dir, tar, gleam_tar, 4)
    src_epoch6, err_src_epoch6 = create_epochcat(data_dir, tar, gleam_tar, 5)
    src_epoch7, err_src_epoch7 = create_epochcat(data_dir, tar, gleam_tar, 6)
    # plotting SED
    f = CF.sed_fig()
    # f.plot_spectrum(
    #     freq,
    #     mwa_flux_yr1,
    #     err_mwa_flux_yr1,
    #     marker="o",
    #     label=epoch_nms[0],
    #     marker_color=colors[0],
    #     s=75,
    # )
    # f.plot_spectrum(
    #     freq,
    #     mwa_flux_yr2,
    #     err_mwa_flux_yr2,
    #     marker="o",
    #     label=epoch_nms[1],
    #     marker_color=colors[1],
    #     s=75,
    # )

    # f.plot_spectrum(
    #     freq,
    #     src_epoch1,
    #     err_src_epoch1,
    #     marker="o",
    #     label=epoch_nms[2],
    #     marker_color=colors[2],
    #     s=75,
    # )
    # f.plot_spectrum(
    #     freq,
    #     src_epoch2,
    #     err_src_epoch2,
    #     marker="o",
    #     label=epoch_nms[3],
    #     marker_color=colors[3],
    #     s=75,
    # )
    # f.plot_spectrum(
    #     freq,
    #     src_epoch3,
    #     err_src_epoch3,
    #     marker="o",
    #     label=epoch_nms[4],
    #     marker_color=colors[4],
    #     s=75,
    # )
    f.plot_spectrum(
        freq,
        src_epoch4,
        err_src_epoch4,
        marker="o",
        # label=epoch_nms[5],
        marker_color=colors[5],
        alpha=0,
        s=75,
    )
    # f.plot_spectrum(
    #     freq,
    #     src_epoch5,
    #     err_src_epoch5,
    #     marker="o",
    #     label=epoch_nms[6],
    #     marker_color=colors[6],
    #     s=75,
    # )
    # f.plot_spectrum(
    #     freq,
    #     src_epoch6,
    #     err_src_epoch6,
    #     marker="o",
    #     label=epoch_nms[7],
    #     marker_color=colors[7],
    #     s=75,
    # )
    # f.plot_spectrum(
    #     freq,
    #     src_epoch7,
    #     err_src_epoch7,
    #     marker="o",
    #     label=epoch_nms[8],
    #     marker_color=colors[8],
    #     s=75,
    # )
    for i in range(len(epoch_nms)):
        # if epoch_nms[i] == "Apr20":
        #     continue
        # else:
        try:
            sampler = open(f"/data/ATCA/analysis/{tar}/{epoch_nms[i]}/{model}/run1/info/results.json")
            results = json.load(sampler)
            params = results["maximum_likelihood"]["point"]
            yvals = model_info[0](freq_cont, *params)
            f.display_model(freq_cont, yvals, colors[i],lw=3)
        except:
            continue
    # f.plt_mcmcfits(sampler, chosen_model, freq_cont, color=colors[5])
    # f.legend(loc="lower center")
    f.title(gleam_tar)
    f.format(xunit="GHz")
    f.save(f"{save_dir}/{tar}_sed_models", ext="png")
    plt.clf()
    plt.close
    return


def plot_epochsed(
    save_dir,
    freq,
    src_flux,
    err_src_flux,
    model_info,
    color,
    tar,
    model,
):
    model_funct = model_info[0]
    # src_flux, err_src_flux = create_epochcat(data_dir, tar, gleam_tar, epoch)
    # plotting SED
    f = CF.sed_fig()
    f.plot_spectrum(
        freq,
        src_flux,
        err_src_flux,
        marker="o",
        marker_color='k',
        s=75,
    )
    epochs = [
        "2013",
        "2014",
        "Jan20",
        "Mar20",
        "Apr20",
        "May20",
        "July20",
        "Oct20",
    ]
    for i in range(len(epochs)):
        sampler = open(f"/data/ATCA/analysis/{tar}/{epochs[i]}/{model}/run1/info/results.json")
        # sequence, final = ultranest.integrator.read_file(
        #     f"/data/ATCA/analysis/{tar}/{epochs[i]}/{model}/run1/",
        #     len(model_info[2:]),
        #     check_insertion_order=False,
        # )
        # band = PredictionBand(freq_cont)
        # for params in final["samples"]:
        #     band.add(model_funct(freq_cont, *params))
        params = sampler["maximum_likelihood"]["point"]
        yvals = model_funct(freq_cont, *params)
        f.display_model(freq_cont, color[i])
    # f.plt_mcmcfits(sampler, chosen_model, freq_cont, color=colors[epoch + 2])
    # f.legend(loc="lower center")
    f.title(f"{tar}")
    # f.format(xunit="GHz")
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


def model_comparison(bayes_name, model1_lnlike, model2_lnlike):
    K_factor = np.exp(model2_lnlike - model1_lnlike)
    print(f"{bayes_name} factor = {K_factor}")
    print(
        f"The first model {K_factor} times more probable than the second model for {bayes_name}"
    )
    return K_factor


def plot_paramswithtime(
    directory,
    target,
    model,
    labels,
    colors=cmr.take_cmap_colors(
        "cmr.rainforest", 12, cmap_range=(0.15, 0.85), return_fmt="hex"
    ),
):
    params = []
    errlo_params = []
    errup_params = []
    months = [-5, -4, 0, 1, 5, 7, 10]
    for epoch in ["2013", "2014", "Jan20", "Mar20", "May20", "July20", "Oct20"]:
        sampler = open(
            f"/data/ATCA/analysis/{target}/{epoch}/{model}/run1/info/results.json"
        )
        results = json.load(sampler)
        paramnames = results["paramnames"]
        parameters = np.array(results["maximum_likelihood"]["point"])
        errlo_arr = np.array(results["posterior"]["errlo"])
        errup_arr = np.array(results["maximum_likelihood"]["point"])
        errlo = parameters - errlo_arr
        errup = errup_arr - parameters
        params.append(parameters)
        errlo_params.append(errlo)
        errup_params.append(errup)
    paramsT = np.transpose(params)
    errloT = np.transpose(errlo_params)
    # errupT = np.transpose(errup_params)
    for i in range(len(paramnames)):
        f = CF.timeseries()
        name = paramnames[i]
        f.plot_params(months, paramsT[i], err_params=errloT[i], s=75)
        f.format()
        f.title(f"{target} {name}")
        f.save(f"{directory}_{name}_{model}", ext="png")
    return


def plot_logzvsmodel(directory, logz, epoch_nms, model_nms):
    num_models = np.arange(np.shape(logz)[1])
    colors = cmr.take_cmap_colors(
        "cmr.bubblegum", len(num_models), cmap_range=(0.15, 0.85), return_fmt="hex"
    )
    print(num_models)
    # plot_logz = np.arange(np.shape(logz)[1])
    for i in range(len(logz)):
        epoch_logz = logz[i]
        max_logz = epoch_logz[5]
        fig = plt.figure(figsize=(10, 8))
        model_logz = np.exp(max_logz - epoch_logz)
        # model_logz[0] = np.nan
        print(model_logz)
        # plot_logz = model_logz
        # print(plot_logz)
        for j in range(5, len(num_models)):
            plt.scatter(
                num_models[j], model_logz[j], label=model_nms[j], color=colors[j]
            )
        plt.legend()
        plt.title(f"Model Comparison of Likelihood {epoch_nms[i]}")
        plt.savefig(f"{directory}{epoch_nms[i]}_logz.png")
    return
