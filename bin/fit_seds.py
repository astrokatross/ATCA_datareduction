#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import fitfuncts
import numpy as np
import gpscssmodels
import cmasher as cmr

NUM_COLORS = 8
colors = cmr.take_cmap_colors(
    "cmr.rainforest", NUM_COLORS, cmap_range=(0.15, 0.85), return_fmt="hex"
)

freq = [
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
freq_cont = np.linspace(0.01, 15, num=10000)
epochs = ["epoch1", "epoch2", "epoch3", "epoch4", "epoch5", "epoch6"]
epoch_nms = ("2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20")

# Initial conditions
nwalkers = 32
ndim = 3
save_dir = "/data/ATCA/analysis/"
data_dir = "/data/ATCA/ATCA_datareduction/"
gleam_tar = "GLEAM J020507-110922"
tar = "J020507"
chosen_model = gpscssmodels.singSSA
epoch = 3
epoch_nm = epoch_nms[epoch + 2]
model_nm = "singSSA"
labels = ["Snorm", "beta", "freqpeak"]

# Calculating a good starting point using scipy.stats.opt.curve_fit or whatever it is
src_flux, err_src_flux = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, epoch)

src_epoch1, err_src_epoch1 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 0)
src_epoch2, err_src_epoch2 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 1)
src_epoch3, err_src_epoch3 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 2)
src_epoch4, err_src_epoch4 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 3)
src_epoch5, err_src_epoch5 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 4)
src_epoch6, err_src_epoch6 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 5)
print(src_flux)

# Uncomment if it's epochs without ATCA
src_epoch5_fit = np.hstack((src_epoch5[0:20], src_epoch4[20:37]))
err_src_epoch5_fit = np.hstack((err_src_epoch5[0:20], err_src_epoch4[20:37]))
src_epoch6_fit = np.hstack((src_epoch6[0:20], src_epoch4[20:37]))
err_src_epoch6_fit = np.hstack((err_src_epoch6[0:20], err_src_epoch4[20:37]))


poptstart, pcovstart = fitfuncts.fit_model(freq, src_flux, err_src_flux, chosen_model)
p0 = poptstart + 1e-4 * np.random.randn(nwalkers, ndim)

sampler = fitfuncts.run_mcmc(
    nwalkers, ndim, p0, freq, src_flux, err_src_flux, chosen_model, niters=5000
)


tau = sampler.get_autocorr_time()
print(tau)
flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)


# Checking how long before it forgets the initial conditions, to be cut kinda like the burn in
fitfuncts.plot_mcmcqualitycheck(
    f"{save_dir}{tar}_{epoch_nm}_{model_nm}", f"{tar} {epoch_nm} {model_nm}", sampler, poptstart, labels
)

fitfuncts.plot_sed(save_dir, data_dir, freq_cont, freq, gleam_tar, tar, colors)
fitfuncts.plot_epochsed(
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
)
