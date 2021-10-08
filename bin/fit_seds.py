#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import fitfuncts
import numpy as np
import gpscssmodels
import priortransfuncts
import cmasher as cmr
import ultranest
from ultranest.plot import PredictionBand
import matplotlib.pylab as plt
import json
import os

num_colors = 12
colors = cmr.take_cmap_colors(
    "cmr.rainforest", num_colors, cmap_range=(0.15, 0.85), return_fmt="hex"
)
freq_cont = np.linspace(0.01, 15, num=10000)
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
epoch_nms = ("2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20")
model_params_dict = {
    "singhomobremss": [
        gpscssmodels.singhomobremss,
        priortransfuncts.singhomobremss,
        "Snorm",
        "alpha",
        "freqpeak",
    ],
    "singinhomobremss": [
        gpscssmodels.singinhomobremss,
        priortransfuncts.singinhomobremss,
        "Snorm",
        "alpha",
        "p",
        "freqpeak",
    ],
    "internalbremss": [
        gpscssmodels.internalbremss,
        priortransfuncts.internalbremss,
        "Snorm",
        "alpha",
        "freqpeak",
    ],
    "singSSA": [
        gpscssmodels.singSSA,
        priortransfuncts.singSSA,
        "Snorm",
        "beta",
        "peakfreq",
    ],
    "singhomobremsscurve": [
        gpscssmodels.singhomobremsscurve,
        priortransfuncts.singhomobremsscurve,
        "Snorm",
        "alpha",
        "freqpeak",
        "q",
    ],
    "singinhomobremsscurve": [
        gpscssmodels.singinhomobremsscurve,
        priortransfuncts.singinhomobremsscurve,
        "Snorm",
        "alpha",
        "p",
        "freqpeak",
        "q",
    ],
    "singinhomobremssbreakexp": [
        gpscssmodels.singinhomobremssbreakexp,
        priortransfuncts.singinhomobremssbreakexp,
        "Snorm",
        "alpha",
        "p",
        "freqpeak",
        "breakfreq",
    ],
    "singhomobremssbreakexp": [
        gpscssmodels.singhomobremssbreakexp,
        priortransfuncts.singhomobremssbreakexp,
        "Snorm",
        "alpha",
        "freqpeak",
        "breakfreq",
    ],
    "singSSAbreakexp": [
        gpscssmodels.singSSAbreakexp,
        priortransfuncts.singSSAbreakexp,
        "Snorm",
        "beta",
        "peakfreq",
        "breakfreq",
    ],
}

# Source/run information
save_dir = "/data/ATCA/analysis/"
data_dir = "/data/ATCA/ATCA_datareduction/"
gleam_tar = "GLEAM J020507-110922"
target = "J020507"
fit_models = [
    "singhomobremss",
    "singinhomobremss",
    "internalbremss",
    "singSSA",
    "singhomobremsscurve",
    "singinhomobremsscurve",
    "singinhomobremssbreakexp",
    "singhomobremssbreakexp",
    "singSSAbreakexp",
]
epochs = [0, 1, 2, 3, 4, 5, 6, 7]


fitfuncts.plot_sed(
    f"{save_dir}/{target}", data_dir, frequency, gleam_tar, target, colors
)
src_epoch1, err_src_epoch1 = fitfuncts.create_epochcat(data_dir, target, gleam_tar, 0)
src_epoch2, err_src_epoch2 = fitfuncts.create_epochcat(data_dir, target, gleam_tar, 1)
src_epoch3, err_src_epoch3 = fitfuncts.create_epochcat(data_dir, target, gleam_tar, 2)
src_epoch4, err_src_epoch4 = fitfuncts.create_epochcat(data_dir, target, gleam_tar, 3)
src_epoch5, err_src_epoch5 = fitfuncts.create_epochcat(data_dir, target, gleam_tar, 4)
src_epoch6, err_src_epoch6 = fitfuncts.create_epochcat(data_dir, target, gleam_tar, 5)

logz = np.zeros((len(fit_models), len(epochs)))
for i in range(len(fit_models)):
    model = fit_models[i]
    model_funct = model_params_dict[model][0]
    model_trans = model_params_dict[model][1]
    labels = model_params_dict[model][2:]

    maxlike_params = []
    err_maxlike_params = []
    model_logz = []
    for j in range(len(epochs)):
        epoch = epochs[j]
        epoch_nm = epoch_nms[epoch]
        color = colors[epoch]
        print(f"Now fitting for epoch {epoch_nm}")
        # Reading in the flux for this epoch (note there accommodations for epochs that don't have atca)
        if epoch == 0:
            (
                mwa_flux_yr1,
                err_mwa_flux_yr1,
                mwa_flux_yr2,
                err_mwa_flux_yr2,
                fluxes_extra,
            ) = fitfuncts.read_gleam_fluxes("/data/MWA", gleam_tar)
            src_flux = np.hstack((mwa_flux_yr1, src_epoch4[20:37]))
            err_src_flux = np.hstack((err_mwa_flux_yr1, err_src_epoch4[20:37]))
        elif epoch == 1:
            (
                mwa_flux_yr1,
                err_mwa_flux_yr1,
                mwa_flux_yr2,
                err_mwa_flux_yr2,
                fluxes_extra,
            ) = fitfuncts.read_gleam_fluxes("/data/MWA", gleam_tar)
            src_flux = np.hstack((mwa_flux_yr2, src_epoch4[20:37]))
            err_src_flux = np.hstack((err_mwa_flux_yr2, err_src_epoch4[20:37]))
        elif epoch == 2:
            (
                mwa_flux_yr1,
                err_mwa_flux_yr1,
                mwa_flux_yr2,
                err_mwa_flux_yr2,
                fluxes_extra,
            ) = fitfuncts.read_gleam_fluxes("/data/MWA", gleam_tar)
            src_flux = np.hstack((mwa_flux_yr1, src_epoch1[20:37]))
            err_src_flux = np.hstack((err_mwa_flux_yr1, err_src_epoch1[20:37]))
        elif epoch == 3:
            (
                mwa_flux_yr1,
                err_mwa_flux_yr1,
                mwa_flux_yr2,
                err_mwa_flux_yr2,
                fluxes_extra,
            ) = fitfuncts.read_gleam_fluxes("/data/MWA", gleam_tar)
            src_flux = np.hstack((mwa_flux_yr1, src_epoch2[20:37]))
            err_src_flux = np.hstack((err_mwa_flux_yr1, err_src_epoch2[20:37]))
        elif epoch == 6:
            src_flux = np.hstack((src_epoch5[0:20], src_epoch4[20:37]))
            err_src_flux = np.hstack((err_src_epoch5[0:20], err_src_epoch4[20:37]))
        elif epoch == 7:
            src_flux = np.hstack((src_epoch6[0:20], src_epoch4[20:37]))
            err_src_flux = np.hstack((err_src_epoch6[0:20], err_src_epoch4[20:37]))
        elif epoch == 4:
            src_flux = src_epoch3
            err_src_flux = err_src_epoch3
        elif epoch == 5: 
            src_flux = src_epoch4
            err_src_flux = err_src_epoch4

        # Making sure there's no nan's in flux
        mask = np.where(~np.isnan(src_flux))
        src_flux = src_flux[mask]
        err_src_flux = err_src_flux[mask]
        freq = frequency[mask]

        try:
            print("Found results of run, continuing with analysis")
            sampler = open(
                f"{save_dir}{target}/{epoch_nm}/{model}/run1/info/results.json"
            )
            print("Found results of run, continuing with analysis")
        except:
            sampler = fitfuncts.run_ultranest_mcmc(
                f"{save_dir}{target}/{epoch_nm}/{model}",
                labels,
                freq,
                src_flux,
                err_src_flux,
                model_funct,
                model_trans,
            )
            sampler.run(max_iters=50000)
            print(f"Finished fitting for {model}")
            sampler.store_tree()
            sampler.plot()
            sampler = open(
                f"/data/ATCA/analysis/{target}/{epoch_nm}/{model}/run1/info/results.json"
            )
        
        if os.path.exists(f"{save_dir}/{target}/seds/{target}_{epoch}_{model}_sed.png")==False:
            sequence, final = ultranest.integrator.read_file(f"{save_dir}{target}/{epoch_nm}/{model}/run1/", len(labels), check_insertion_order=False)

            # TODO: figure out how to generalise the band plotting thing
            band = PredictionBand(freq_cont)
            for params in final["samples"]:
                band.add(model_funct(freq_cont, *params))

            fitfuncts.plot_epochsed(
                f"{save_dir}/{target}/seds/{target}_{epoch_nm}_{model}",
                freq,
                src_flux,
                err_src_flux,
                band,
                color,
                target,
            )

        results = json.load(sampler)
        maxlike_params.append(results["maximum_likelihood"]["point"])
        err_maxlike_params.append(results["posterior"]["stdev"])
        logz[i][j] = results["logz"]
    fitfuncts.plot_paramswithtime(
        f"{save_dir}/{target}/{target}", target, model, labels
    )
        # # homoFFA inhomoFFA
        # K_homoinhomo = fitfuncts.model_comparison(
        #     "homovsinhomoFFA",
        #     sampler_singhomobremss[1]["logz"],
        #     sampler_singinhomobremss[1]["logz"],
        # )


