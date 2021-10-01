#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import pandas as pd
import scipy.optimize as opt
import numpy as np
import CFigTools.CustomFigure as CF
import gpscssmodels
import emcee
import corner
import matplotlib.pyplot as plt

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


def create_epochcat(directory, tar, gleam_tar, epoch):
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
    axes[0].set_title(f"{figname}")
    plt.savefig(f"{directory}_paramchains.png", overwrite=True)
    plt.clf()
    plt.close()
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    # Plotting the corner plot of the different variables
    corner.corner(flat_samples, labels=labels, truths=p0)
    plt.title(f"{figname}")
    plt.savefig(f"{directory}_corner.png", overwrite=True)
    plt.clf()
    plt.close()
    return


def plot_sed(save_dir, data_dir, freq_cont, freq, gleam_tar, tar, colors):
    # mwa_flux, model_vals, xtra_values, values = read_gleam_fluxes(
    #     "/data/MWA", gleam_tar
    # )
    src_epoch1, err_src_epoch1 = create_epochcat(data_dir, tar, gleam_tar, 0)
    src_epoch2, err_src_epoch2 = create_epochcat(data_dir, tar, gleam_tar, 1)
    src_epoch3, err_src_epoch3 = create_epochcat(data_dir, tar, gleam_tar, 2)
    src_epoch4, err_src_epoch4 = create_epochcat(data_dir, tar, gleam_tar, 3)
    src_epoch5, err_src_epoch5 = create_epochcat(data_dir, tar, gleam_tar, 4)
    src_epoch6, err_src_epoch6 = create_epochcat(data_dir, tar, gleam_tar, 5)
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
    f.legend(loc="lower center")
    f.title(gleam_tar)
    f.format(xunit="GHz")
    f.save(f"{save_dir}/{tar}_sed", ext="png")
    return


def plot_epochsed(
    save_dir,
    data_dir,
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
    epoch_nm = epoch_nms[epoch+2]
    src_flux, err_src_flux = create_epochcat(data_dir, tar, gleam_tar, epoch)
    # plotting SED
    f = CF.sed_fig()
    f.plot_spectrum(
        freq,
        src_flux,
        err_src_flux,
        marker="o",
        label=f"{model_nm}",
        marker_color=colors[epoch+2],
        s=75,
    )
    f.plt_mcmcfits(sampler, chosen_model, freq_cont, color=colors[epoch+2])
    # f.legend(loc="lower center")
    f.title(f"{gleam_tar} {epoch_nm} {model_nm}")
    f.format(xunit="GHz")
    f.save(f"{save_dir}/{tar}_{epoch_nm}_{model_nm}_sed", ext="png")
    return
