#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

from numpy.lib.npyio import save
import fitfuncts
import numpy as np
import gpscssmodels
import priortransfuncts
import cmasher as cmr
import ultranest
from ultranest.plot import PredictionBand
import matplotlib.pylab as plt
import json

num_colors = 12
colors = cmr.take_cmap_colors(
    "cmr.rainforest", num_colors, cmap_range=(0.15, 0.85), return_fmt="hex"
)
freq_cont = np.linspace(0.01, 15, num=10000)
freq = np.array(
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
epochs = [0, 1, 2, 3, 4, 5]


fitfuncts.plot_sed(f"{save_dir}/{target}", data_dir, freq, gleam_tar, target, colors)

for model in fit_models:
    model_funct = model_params_dict[model][0]
    model_trans = model_params_dict[model][1]
    labels = model_params_dict[model][2:]


    maxlike_params = []
    err_maxlike_params = []
    for epoch in epochs:
        epoch_nm = epoch_nms[epoch]
        color = colors[epoch]

        # Reading in the flux for this epoch (note there accommodations for epochs that don't have atca)
        src_epoch4, err_src_epoch4 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 3)
        src_epoch5, err_src_epoch5 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 4)
        src_epoch6, err_src_epoch6 = fitfuncts.create_epochcat(data_dir, tar, gleam_tar, 5)
        if epoch == 4:
            src_flux = np.hstack((src_epoch5[0:20], src_epoch4[20:37]))
            err_src_flux = np.hstack((err_src_epoch5[0:20], err_src_epoch4[20:37]))
        elif epoch == 5:
            src_flux = np.hstack((src_epoch6[0:20], src_epoch4[20:37]))
            err_src_flux = np.hstack((err_src_epoch6[0:20], err_src_epoch4[20:37]))
        elif epoch == 0:
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
        else:
            src_flux, err_src_flux = fitfuncts.create_epochcat(
                data_dir, tar, gleam_tar, epoch
            )

        # Making sure there's no nan's in flux
        mask = np.where(~np.isnan(src_flux))
        src_flux = src_flux[mask]
        err_src_flux = err_src_flux[mask]
        freq = freq[mask]

        try:
            sampler = open(
                f"/data/ATCA/analysis/{target}/{epoch_nm}/{model}/run1/info/results.json"
            )
        except (FileNotFoundError, KeyError):
            sampler = fitfuncts.run_ultranest_mcmc(
                f"/data/ATCA/analysis/{target}/{epoch_nm}/{model}",
                labels,
                freq,
                src_flux,
                err_src_flux,
                model,
                model_trans,
            )
            sampler.run(max_iters=50000)
            print(f"Finished fitting for {model}")
            sampler.store_tree()
            sampler.plot()

            # TODO: figure out how to generalise the band plotting thing 
            # band = PredictionBand(freq_cont)
            # for Snorm, alpha, freqpeak in sampler.results["samples"]:
            #     band.add(chosen_model(freq_cont, Snorm, alpha, freqpeak))

            # fitfuncts.plot_epochsed(
            #     f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
            #     freq,
            #     src_flux,
            #     err_src_flux,
            #     band,
            #     color,
            #     tar,
            # )

            sampler = open(
                f"/data/ATCA/analysis/{target}/{epoch_nm}/{model}/run1/info/results.json"
            )

        results = json.load(sampler)
        maxlike_params.append(results["maximum_likelihood"]["point"])
        err_maxlike_params.append(results["posterior"]["stdev"])
        

# ------------------------------------------------------------------------------
#  Running initial no break models: singSSA, singhomobremss, singinhomobremss
# ------------------------------------------------------------------------------

# # Initial conditions of mcmc for singSSA
# chosen_model = gpscssmodels.singSSA
# transform_funct = fitfuncts.priortrans_singSSA
# model_nm = "singSSA"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_singssa = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_singssa[1]["logz"])
# except:
        # sampler_singssa = fitfuncts.run_ultranest_mcmc(
        #     directory,
        #     labels,
        #     freq,
        #     src_flux,
        #     err_src_flux,
        #     chosen_model,
        #     transform_funct,
        #     resume="overwrite",
        # )
        # sampler_singssa.run(max_iters=50000)
        # print(f"Finished fitting for {model_nm}")
        # sampler_singssa.store_tree()
        # sampler_singssa.plot()
        # band = PredictionBand(freq_cont)
        # for Snorm, alpha, freqpeak in sampler_singssa.results["samples"]:
        #     band.add(chosen_model(freq_cont, Snorm, alpha, freqpeak))

        # fitfuncts.plot_epochsed(
        #     f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
        #     freq,
        #     src_flux,
        #     err_src_flux,
        #     band,
        #     color,
        #     tar,
        # )


# # Initial conditions of mcmc for singhomo
# chosen_model = gpscssmodels.singhomobremss
# transform_funct = fitfuncts.priortrans_singhomobremss
# model_nm = "singhomobremss"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_singhomobremss = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_singhomobremss[1]["logz"])
# except:
#     sampler_singhomobremss = fitfuncts.run_ultranest_mcmc(
#         directory,
#         labels,
#         freq,
#         src_flux,
#         err_src_flux,
#         chosen_model,
#         transform_funct,
#         resume="overwrite",
#     )
#     sampler_singhomobremss.run(max_iters=50000)
#     print(f"Finished fitting for {model_nm}")
#     sampler_singhomobremss.store_tree()
#     sampler_singhomobremss.plot()
#     band = PredictionBand(freq_cont)
#     for Snorm, alpha, freqpeak in sampler_singhomobremss.results["samples"]:
#         band.add(chosen_model(freq_cont, Snorm, alpha, freqpeak))

#     fitfuncts.plot_epochsed(
#         f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#         freq,
#         src_flux,
#         err_src_flux,
#         band,
#         color,
#         tar,
#     )


# # Initial conditions of mcmc for singhomobremss
# chosen_model = gpscssmodels.singinhomobremss
# transform_funct = fitfuncts.priortrans_singinhomobremss
# model_nm = "singinhomobremss"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_singinhomobremss = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_singinhomobremss[1]["logz"])
# except:
#     sampler_singinhomobremss = fitfuncts.run_ultranest_mcmc(
#         directory,
#         labels,
#         freq,
#         src_flux,
#         err_src_flux,
#         chosen_model,
#         transform_funct,
#         resume="overwrite",
#     )
#     sampler_singinhomobremss.run(max_iters=50000)
#     print(f"Finished fitting for {model_nm}")
#     sampler_singinhomobremss.store_tree()
#     sampler_singinhomobremss.plot()
#     band = PredictionBand(freq_cont)
#     for Snorm, alpha, p, freqpeak in sampler_singinhomobremss.results["samples"]:
#         band.add(chosen_model(freq_cont, Snorm, alpha, p, freqpeak))

#     fitfuncts.plot_epochsed(
#         f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#         freq,
#         src_flux,
#         err_src_flux,
#         band,
#         color,
#         tar,
#     )


# # ------------------------------------------------------------------------------
# #  Introducing a break at higherfrequencies: singSSA
# # ------------------------------------------------------------------------------


# # Initial conditions of mcmc for singSSAbreakexp
# chosen_model = gpscssmodels.singSSAbreakexp
# transform_funct = fitfuncts.priortrans_singSSAbreakexp
# model_nm = "singSSAbreakexp"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_singssabreak = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_singssabreak[1]["logz"])
# except:
#     sampler_singssabreak = fitfuncts.run_ultranest_mcmc(
#         directory,
#         labels,
#         freq,
#         src_flux,
#         err_src_flux,
#         chosen_model,
#         transform_funct,
#         resume="overwrite",
#     )
#     sampler_singssabreak.run(max_iters=50000)
#     print(f"Finished fitting for {model_nm}")
#     sampler_singssabreak.store_tree()
#     sampler_singssabreak.plot()
#     band = PredictionBand(freq_cont)
#     for Snorm, alpha, freqpeak, breakfreq in sampler_singssabreak.results["samples"]:
#         band.add(chosen_model(freq_cont, Snorm, alpha, freqpeak, breakfreq))

#     fitfuncts.plot_epochsed(
#         f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#         freq,
#         src_flux,
#         err_src_flux,
#         band,
#         color,
#         tar,
#     )


# # ------------------------------------------------------------------------------
# #  Introducing a break at higherfrequencies: homo
# # ------------------------------------------------------------------------------


# # Initial conditions of mcmc for homobreak
# chosen_model = gpscssmodels.singhomobremssbreak
# transform_funct = fitfuncts.priortrans_singhomobremssbreak
# model_nm = "singhomobremssbreak"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_homobreak = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_homobreak[1]["logz"])
# except:
#     sampler_homobreak = fitfuncts.run_ultranest_mcmc(
#         directory,
#         labels,
#         freq,
#         src_flux,
#         err_src_flux,
#         chosen_model,
#         transform_funct,
#         resume="overwrite",
#     )
#     sampler_homobreak.run(max_iters=50000)
#     print(f"Finished fitting for {model_nm}")

#     sampler_homobreak.store_tree()
#     sampler_homobreak.plot()
#     band = PredictionBand(freq_cont)
#     for Snorm, alpha, freqpeak, breakfreq in sampler_homobreak.results["samples"]:
#         band.add(chosen_model(freq_cont, Snorm, alpha, freqpeak, breakfreq))

#     fitfuncts.plot_epochsed(
#         f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#         freq,
#         src_flux,
#         err_src_flux,
#         band,
#         color,
#         tar,
#     )

# # Initial conditions of mcmc for homobreakexp
# chosen_model = gpscssmodels.singhomobremssbreakexp
# transform_funct = fitfuncts.priortrans_singhomobremssbreak
# model_nm = "singhomobremssbreakexp"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")


# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_singhomobremssbreakexp = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_singhomobremssbreakexp[1]["logz"])
# except:
#     sampler_singhomobremssbreakexp = fitfuncts.run_ultranest_mcmc(
#         directory,
#         labels,
#         freq,
#         src_flux,
#         err_src_flux,
#         chosen_model,
#         transform_funct,
#         resume="overwrite",
#     )
#     sampler_singhomobremssbreakexp.run(max_iters=50000)
#     print(f"Finished fitting for {model_nm}")
#     sampler_singhomobremssbreakexp.store_tree()
#     sampler_singhomobremssbreakexp.plot()
#     band = PredictionBand(freq_cont)
#     for Snorm, alpha, freqpeak, breakfreq in sampler_singhomobremssbreakexp.results[
#         "samples"
#     ]:
#         band.add(chosen_model(freq_cont, Snorm, alpha, freqpeak, breakfreq))

#     fitfuncts.plot_epochsed(
#         f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#         freq,
#         src_flux,
#         err_src_flux,
#         band,
#         color,
#         tar,
#     )


# # ------------------------------------------------------------------------------
# #  Introducing a break at higherfrequencies: inhomo
# # ------------------------------------------------------------------------------


# # Initial conditions of mcmc for inhomobreak
# chosen_model = gpscssmodels.singinhomobremssbreak
# transform_funct = fitfuncts.priortrans_singinhomobremssbreak
# model_nm = "singinhomobremssbreak"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_inhomobreak = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_inhomobreak[1]["logz"])
# except:
#     sampler_inhomobreak = fitfuncts.run_ultranest_mcmc(
#         directory,
#         labels,
#         freq,
#         src_flux,
#         err_src_flux,
#         chosen_model,
#         transform_funct,
#         resume="overwrite",
#     )
#     sampler_inhomobreak.run(max_iters=50000)
#     print(f"Finished fitting for {model_nm}")
#     sampler_inhomobreak.store_tree()
#     sampler_inhomobreak.plot()
#     band = PredictionBand(freq_cont)
#     for Snorm, alpha, p, freqpeak, breakfreq in sampler_inhomobreak.results["samples"]:
#         band.add(chosen_model(freq_cont, Snorm, alpha, p, freqpeak, breakfreq))

#     fitfuncts.plot_epochsed(
#         f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#         freq,
#         src_flux,
#         err_src_flux,
#         band,
#         color,
#         tar,
#     )


# # Initial conditions of mcmc for inhomobreakexp
# chosen_model = gpscssmodels.singinhomobremssbreakexp
# transform_funct = fitfuncts.priortrans_singinhomobremssbreak
# model_nm = "singinhomobremssbreakexp"
# labels = model_params_dict[model_nm]
# directory = f"{save_dir}{tar}/{epoch_nm}/{model_nm}/"
# print(f"Fitting for {model_nm}")

# # Trying to see if it's already run the mcmc before and just load it if you have, otherwise runs the mcmc again
# try:
#     sampler_singinhomobremssbreakexp = ultranest.integrator.read_file(
#         f"{directory}/run1/", len(labels), check_insertion_order=False
#     )
#     print(sampler_singinhomobremssbreakexp[1]["logz"])
# except:
# sampler_singinhomobremssbreakexp = fitfuncts.run_ultranest_mcmc(
#     directory,
#     labels,
#     freq,
#     src_flux,
#     err_src_flux,
#     chosen_model,
#     transform_funct,
#     resume="overwrite",
# )
# sampler_singinhomobremssbreakexp.run(max_iters=50000)
# print(f"Finished fitting for {model_nm}")
# sampler_singinhomobremssbreakexp.store_tree()
# sampler_singinhomobremssbreakexp.plot()
# band = PredictionBand(freq_cont)
# for (
#     Snorm,
#     alpha,
#     p,
#     freqpeak,
#     breakfreq,
# ) in sampler_singinhomobremssbreakexp.results["samples"]:
#     band.add(chosen_model(freq_cont, Snorm, alpha, p, freqpeak, breakfreq))

# fitfuncts.plot_epochsed(
#     f"{save_dir}/{tar}_{epoch_nm}_{model_nm}",
#     freq,
#     src_flux,
#     err_src_flux,
#     band,
#     color,
#     tar,
# )


# ------------------------------------------------------------------------------
#  Comparing models using Bayes factor
# ------------------------------------------------------------------------------

# Section of model comparison
# First comparing generic SSA and FFA
K_ssaffa = fitfuncts.model_comparison(
    "SSAvshomoFFA", sampler_singssa[1]["logz"], sampler_singhomobremss[1]["logz"]
)


# homoFFA inhomoFFA
K_homoinhomo = fitfuncts.model_comparison(
    "homovsinhomoFFA",
    sampler_singhomobremss[1]["logz"],
    sampler_singinhomobremss[1]["logz"],
)


# Homo break or not
K_homobreak = fitfuncts.model_comparison(
    "homovsbreak",
    sampler_singhomobremss[1]["logz"],
    sampler_singhomobremssbreakexp[1]["logz"],
)


# inhomo break or not
K_inhomobreak = fitfuncts.model_comparison(
    "inhomovsbreak",
    sampler_singinhomobremss[1]["logz"],
    sampler_singinhomobremssbreakexp[1]["logz"],
)


# homobreak vs inhomobreak
K_homobreakinhomobreak = fitfuncts.model_comparison(
    "homobreakvsinhomobreak",
    sampler_singhomobremssbreakexp[1]["logz"],
    sampler_singinhomobremssbreakexp[1]["logz"],
)


best_model_results = sampler_singinhomobremssbreakexp[1]
print(best_model_results)
params_bestfit = best_model_results["median"]
# TODO: find a better way to incorporate the asymmetric errors (errlo and errup in the results)
err_params_bestfit = best_model_results["stdev"]

chosen_model = gpscssmodels.singinhomobremssbreakexp
model_nm = "singinhomobremssbreakexp"
model = model_nm
target = tar
labels = model_params_dict[model_nm]
print(f"plotting params for {model_nm}")
fitfuncts.plot_paramswithtime(f"{save_dir}/{tar}/{tar}", tar, model_nm, labels)
