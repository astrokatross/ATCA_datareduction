#!/usr/bin/python3
# This script has functions for analysising the variability and returning useful and interpretable parameters
# By K.Ross 29/09/21

import numpy as np
import gpscssmodels
import json
import cmasher as cmr
import matplotlib.pyplot as plt
import math

# from casatasks import uvmodelfit
import datetime
import os

from casacore.tables import table
import fitfuncts
import CFigTools.CustomFigure as CF
import priortransfuncts

num_colors = 8
colors = cmr.take_cmap_colors(
    "cmr.gothic", num_colors, cmap_range=(0.15, 0.8), return_fmt="hex"
)
freq_cont = np.linspace(0.01, 25, num=10000)
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
epoch_nms = [
    "2013",
    "2014",
    # "2020-01",
    # "2020-03",
    "2020-04",
    "2020-05",
    "2020-07",
    "2020-09",
]
model_params_dict = {
    "FFA": [
        gpscssmodels.singhomobremss,
        priortransfuncts.singhomobremss,
        "Snorm",
        "alpha",
        "peak_frequency",
    ],
    "inFFA": [
        gpscssmodels.singinhomobremss,
        priortransfuncts.singinhomobremss,
        "Snorm",
        "alpha",
        "p",
        "peak_frequency",
    ],
    "SSA": [
        gpscssmodels.singSSA,
        priortransfuncts.singSSA,
        "Snorm",
        "beta",
        "peak_frequency",
    ],
    "inFFAb": [
        gpscssmodels.singinhomobremssbreakexp,
        priortransfuncts.singinhomobremssbreakexp,
        "Snorm",
        "alpha",
        "p",
        "peak_frequency",
        "breakfreq",
    ],
    "FFAb": [
        gpscssmodels.singhomobremssbreakexp,
        priortransfuncts.singhomobremssbreakexp,
        "Snorm",
        "alpha",
        "peak_frequency",
        "breakfreq",
    ],
    "SSAb": [
        gpscssmodels.singSSAbreakexp,
        priortransfuncts.singSSAbreakexp,
        "Snorm",
        "beta",
        "peak_frequency",
        "breakfreq",
    ],
}


def calc_yvals(directory, target, epoch):
    freq_cont = np.linspace(0.01, 25, num=10000)
    yvals = []
    logz = []
    model_nms = [
        "SSA",
        "FFA",
        "inFFA",
        "SSAb",
        "FFAb",
        "inFFAb",
    ]
    models = [
        gpscssmodels.singSSA,
        gpscssmodels.singhomobremss,
        gpscssmodels.singinhomobremss,
        gpscssmodels.singSSAbreakexp,
        gpscssmodels.singhomobremssbreakexp,
        gpscssmodels.singinhomobremssbreakexp,
    ]
    for i in range(len(models)):
        try:
            sampler = open(f"{directory}/{epoch}/{model_nms[i]}/run1/info/results.json")
            results = json.load(sampler)
            param_mod = results["maximum_likelihood"]["point"]
            yvals_mod = models[i](freq_cont, *param_mod)
            yvals.append(yvals_mod)
            logz.append(results["logz"])
        except FileNotFoundError:
            yvals.append(np.full(10000, np.nan))
            logz.append(np.nan)
    return yvals, logz


def calc_modelnparams(directory, target, model):
    freq_cont = np.linspace(0.01, 25, num=10000)
    yvals = []
    nu_p = []
    alpha = []
    errlo_nu_p = []
    errup_nu_p = []
    errlo_alpha = []
    errup_alpha = []
    model_nm = model[0]
    model_funct = model[1]
    epochs = [
        "2013",
        "2014",
        "2020-01",
        "2020-03",
        "2020-04",
        "2020-05",
        "2020-07",
        "2020-09",
    ]
    for i in range(len(epochs)):
        try:
            sampler = open(f"{directory}/{epochs[i]}/{model_nm}/run1/info/results.json")
            results = json.load(sampler)
            param_mod = results["maximum_likelihood"]["point"]
            paramnames = results["paramnames"]
            err_lo = np.array(results["posterior"]["errlo"])
            err_hi = np.array(results["posterior"]["errup"])
            errs_low = param_mod - err_lo
            errs_up = err_hi - param_mod
            yvals_mod = model_funct(freq_cont, *param_mod)
            yvals.append(yvals_mod)
            nu_p.append(param_mod[paramnames.index("peak_frequency")])
            errlo_nu_p.append(abs(errs_low[paramnames.index("peak_frequency")]))
            errup_nu_p.append(abs(errs_up[paramnames.index("peak_frequency")]))
            # try:
            if model_nm in ["SSAb", "SSA"]:
                alpha.append(param_mod[paramnames.index("beta")])
                errlo_alpha.append(abs(errs_low[paramnames.index("beta")]))
                errup_alpha.append(abs(errs_up[paramnames.index("beta")]))
            else:
                alpha.append(param_mod[paramnames.index("alpha")])
                errlo_alpha.append(abs(errs_low[paramnames.index("alpha")]))
                errup_alpha.append(abs(errs_up[paramnames.index("alpha")]))
            # except:
        except FileNotFoundError:
            yvals.append(np.full(len(freq_cont), np.nan))
            nu_p.append(np.nan)
            errlo_nu_p.append(np.nan)
            errup_nu_p.append(np.nan)
            alpha.append(np.nan)
            errlo_alpha.append(np.nan)
            errup_alpha.append(np.nan)
            print("Can't find sampler....")
        err_nu_p = np.stack((errlo_nu_p, errup_nu_p))
        err_alpha = np.stack((errlo_alpha, errup_alpha))
    return yvals, nu_p, err_nu_p, alpha, err_alpha


def read_results(dir):
    try:
        sampler = open(f"{dir}/run1/info/results.json")
        results = json.load(sampler)
        logz = results["logz"]
        return logz
    except:
        print("Can't find sampler")
        return np.nan


def read_alllogz(target, epoch):
    ssa_logz = read_results(f"/data/ATCA/analysis/{target}/{epoch}/SSA")
    ffa_logz = read_results(f"/data/ATCA/analysis/{target}/{epoch}/FFA")
    inffa_logz = read_results(f"/data/ATCA/analysis/{target}/{epoch}/inFFA")
    ssab_logz = read_results(f"/data/ATCA/analysis/{target}/{epoch}/SSAb")
    ffab_logz = read_results(f"/data/ATCA/analysis/{target}/{epoch}/FFAb")
    inffab_logz = read_results(f"/data/ATCA/analysis/{target}/{epoch}/inFFAb")
    return ssa_logz, ffa_logz, inffa_logz, ssab_logz, ffab_logz, inffab_logz


def calc_logzmod(target):
    model_nms = [
        "SSA",
        "FFA",
        "inFFA",
        "SSAb",
        "FFAb",
        "inFFAb",
    ]
    epochs = [
        "2013",
        "2014",
        # "2020-01",
        # "2020-03",
        "2020-04",
        "2020-05",
        "2020-07",
        "2020-09",
    ]
    logz = []
    for i in range(len(epochs)):
        try:
            logz_epoch = read_alllogz(target, epochs[i])
            logz.append(logz_epoch)
        except:
            print("Couldn't read logz")
            logz.append(np.nan)
    logz = np.array(logz)
    sum_logz = np.zeros(len(logz[0]))
    for i in range(len(logz[0])):
        sum_logz[i] = np.nansum(logz[:, i])
    if target in ["J020507", "J024838"]:
        avg_logz = sum_logz / 5
    else:
        avg_logz = sum_logz / 6
    return avg_logz


def calc_secbayes(avg_logz):
    rest_logz = np.delete(avg_logz, np.where(avg_logz == np.max(avg_logz)))

    K = np.exp(np.max((avg_logz) - np.max(rest_logz)))
    print(f"K (best compared to second) = {K:.4f}")
    return


def read_timeranges(start_times, end_times):
    timeranges = []
    for i in range(len(start_times)):
        start_time = datetime.datetime.strptime(start_times[i], "%H:%M:%S")
        end_time = datetime.datetime.strptime(end_times[i], "%H:%M:%S")
        new_time = start_time
        current_time = start_time
        timeranges.append(start_times[i])
        while new_time < end_time:
            new_time = new_time + datetime.timedelta(0, 30)
            current_str = new_time.strftime("%H:%M:%S")
            timeranges.append(current_str)
    return timeranges


def plt_alphatime(
    save_dir,
    alpha,
    err_alpha,
    plot_title,
    ext="pdf",
    colors=cmr.take_cmap_colors(
        "cmr.gothic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
):
    figsize = (3.5, 1.25)
    # fig = plt.figure(1, figsize=figsize, facecolor="white")
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, wspace=0.05, width_ratios=[1, 3])
    ax = gs.subplots()
    fig.suptitle(plot_title, fontsize=8)
    ax[0].set_ylabel(r"$\alpha$", fontsize=8)
    ax[0].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    ax[0].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    ax[0].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    ax[1].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=1.5, pad=5
    )
    ax[1].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    ax[1].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    ax[0].spines["right"].set_visible(False)
    ax[1].spines["left"].set_visible(False)
    ax[0].yaxis.set_ticks_position("left")
    ax[1].yaxis.set_ticks_position("right")
    ax[1].yaxis.set_ticklabels([])
    d = 1
    kwargs = dict(
        marker=[(-1, -d), (1, d)],
        markersize=5,
        linestyle="none",
        color="k",
        mec="k",
        mew=1,
        clip_on=False,
    )
    ax[0].plot([1, 1], [0, 1], transform=ax[0].transAxes, **kwargs)
    ax[1].plot([0, 0], [0, 1], transform=ax[1].transAxes, **kwargs)

    months = [-5, -4, 4, 5, 7, 10]
    for i in range(len(months)):
        ax[0].errorbar(
            months[i],
            alpha[i],
            yerr=np.array([[err_alpha[0, i], err_alpha[1, i]]]).T,
            fmt="o",
            color=colors[i],
            capsize=3,
            markersize=5,
            elinewidth=0.5,
        )
        ax[1].errorbar(
            months[i],
            alpha[i],
            yerr=np.array([[err_alpha[0, i], err_alpha[1, i]]]).T,
            fmt="o",
            color=colors[i],
            capsize=3,
            markersize=5,
            elinewidth=0.5,
        )
    # for axis in ["top", "bottom", "left", "right"]:
    #     ax[0].spines[axis].set_linewidth(2)
    #     ax[1].spines[axis].set_linewidth(2)
    ax[0].set_xticks(months)
    ax[1].set_xticks(months)
    ax[0].set_xticklabels(["2013", "2014", "Apr20", "May20", "Jul20", "Sept20"])
    ax[1].set_xticklabels(["2013", "2014", "Apr20", "May20", "Jul20", "Sept20"])
    ax[0].set_xlim([-6, -3])
    ax[1].set_xlim([0, 12])
    plt.tight_layout()

    plt.savefig(f"{save_dir}_alphavstime.{ext}", overwrite=True)
    plt.clf()
    plt.close()
    return


def plt_peakftime(
    save_dir,
    nu_p,
    err_nu_p,
    plot_title,
    ext="pdf",
    colors=cmr.take_cmap_colors(
        "cmr.gothic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
):
    figsize = (3.5, 1.5)
    # fig = plt.figure(1, figsize=figsize, facecolor="white")
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, wspace=0.05, width_ratios=[1, 3])
    ax = gs.subplots()
    fig.suptitle(plot_title, fontsize=8)
    ax[0].set_ylabel(r"$\nu_p$ (GHz)", fontsize=8)
    # ax[1].set_xlabel(r"Date", fontsize=8)
    fig.text(
        0.5, -0.06, f"Date", ha="center", fontsize=8
    )
    ax[0].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    ax[0].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    ax[0].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    ax[1].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    ax[1].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    ax[1].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    ax[0].spines["right"].set_visible(False)
    ax[1].spines["left"].set_visible(False)
    ax[0].yaxis.set_ticks_position("left")
    ax[1].yaxis.set_ticks_position("right")
    ax[1].yaxis.set_ticklabels([])
    d = 1
    kwargs = dict(
        marker=[(-1, -d), (1, d)],
        markersize=3,
        linestyle="none",
        color="k",
        mec="k",
        mew=1,
        clip_on=False,
    )
    ax[0].plot([1, 1], [0, 1], transform=ax[0].transAxes, **kwargs)
    ax[1].plot([0, 0], [0, 1], transform=ax[1].transAxes, **kwargs)

    months = [-5, -4, 1, 3, 4, 5, 7, 9]
    for i in range(len(months)):
        ax[0].errorbar(
            months[i],
            nu_p[i],
            yerr=np.array([[err_nu_p[0, i], err_nu_p[1, i]]]).T,
            fmt="o",
            color=colors[i],
            capsize=1,
            markersize=2,
            elinewidth=0.5,
        )
        ax[1].errorbar(
            months[i],
            nu_p[i],
            yerr=np.array([[err_nu_p[0, i], err_nu_p[1, i]]]).T,
            fmt="o",
            color=colors[i],
            capsize=1,
            markersize=2,
            elinewidth=0.5,
        )
    # for axis in ["top", "bottom", "left", "right"]:
    #     ax[0].spines[axis].set_linewidth(2)
    #     ax[1].spines[axis].set_linewidth(2)
    ax[0].set_xticks([-5, -4, 4, 5, 7, 9])
    ax[1].set_xticks([-5, -4, 4, 5, 7, 9])
    ax[0].set_xticklabels(["'13", "'14", "Apr20", "May20", "Jul20", "Sept20"])
    ax[1].set_xticklabels(["'13", "'14", "Apr20", "May20", "Jul20", "Sept20"])
    ax[0].set_xlim([-6, -3])
    ax[1].set_xlim([3.5, 9.5])
    # plt.tight_layout()

    plt.savefig(f"{save_dir}_peakfreqvstime.{ext}", overwrite=True, bbox_inches='tight')
    plt.clf()
    plt.close()
    return


# def read_lightcurveflux(data_dir, outfile_dir, timeranges):
#     fluxes = []
#     err_fluxes = []
#     for i in range(len(timeranges)):
#         timerange = f"{timeranges[i]}+00:00:30"
#         outfile = f"{outfile_dir}_{timerange}.cl"
#         tar_ms = f"{data_dir}_selfcal.ms"
#         print(tar_ms)
#         if (os.path.exists(outfile)) is False:
#             uvmodelfit(
#                 vis=tar_ms,
#                 niter=10,
#                 comptype="P",
#                 outfile=outfile,
#                 field="0",
#                 selectdata=True,
#                 timerange=timerange,
#             )
#             tbl = table(outfile)
#             flux = tbl.getcell("Flux", 0)[0].astype("float64")
#             err = np.sqrt((flux * 0.05) ** 2 + (0.002 ** 2))
#             fluxes.append(flux)
#             err_fluxes.append(err)
#         else:
#             tbl = table(outfile)
#             flux = tbl.getcell("Flux", 0)[0].astype("float64")
#             err = np.sqrt((0.002 ** 2))# + (flux * 0.05) ** 2)
#             fluxes.append(flux)
#             err_fluxes.append(err)
#     print(np.std(fluxes)/np.median(fluxes))
#     mod = round(np.std(fluxes)/np.median(fluxes), -int(math.floor(math.log10(abs(np.std(fluxes)/np.median(fluxes))))))
#     err_fluxes = (err_fluxes/np.median(fluxes))*100
#     fluxes = ((fluxes/np.median(fluxes)) - 1)*100
#     return fluxes, err_fluxes, mod


def plt_lightcurve_continual(
    save_dir,
    scan_times_src1,
    src1_fluxes,
    err_src1_fluxes,
    mod_src1,
    scan_times_src2,
    src2_fluxes,
    err_src2_fluxes,
    mod_src2,
    ext="pdf",
    src1_nm="J215436",
    src2_nm="2211-388",
):
    src1_fluxes_c = src1_fluxes[0]
    src1_fluxes_x = src1_fluxes[1]
    err_src1_fluxes_c = err_src1_fluxes[0]
    err_src1_fluxes_x = err_src1_fluxes[1]

    mod_src1_c = mod_src1[0]
    mod_src1_x = mod_src1[1]

    src2_fluxes_c = src2_fluxes[0]
    src2_fluxes_x = src2_fluxes[1]
    err_src2_fluxes_c = err_src2_fluxes[0]
    err_src2_fluxes_x = err_src2_fluxes[1]

    mod_src2_c = mod_src2[0]
    mod_src2_x = mod_src2[1]

    fig = plt.figure(figsize=(3.5, 2.2))
    gs = fig.add_gridspec(2, 1, hspace=0, wspace=0.05)
    ax = gs.subplots(sharey=True)
    fig.suptitle(f"Light Curve 2021-07-24", fontsize=8)

    axc = ax[0]
    axx = ax[1]
    axc.errorbar(
        scan_times_src1,
        src1_fluxes_c,
        err_src1_fluxes_c,
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_c}",
        linestyle="None",
        marker=".",
        markersize=5,
    )

    axx.errorbar(
        scan_times_src1,
        src1_fluxes_x,
        err_src1_fluxes_x,
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_x}",
        linestyle="None",
        marker=".",
        markersize=5,
    )
    axc.errorbar(
        scan_times_src2,
        src2_fluxes_c,
        err_src2_fluxes_c,
        color="C6",
        label=f"{src2_nm}, $m=$ {mod_src2_c}",
        linestyle="None",
        marker=".",
        markersize=5,
    )

    axx.errorbar(
        scan_times_src2,
        src2_fluxes_x,
        err_src2_fluxes_x,
        color="C6",
        label=f"{src2_nm}, $m=$ {mod_src2_x}",
        linestyle="None",
        marker=".",
        markersize=5,
    )

    axc.set_ylabel(r"$S_{5.5\mathrm{GHz}}$ Offset (\%)", fontsize=8)
    axx.set_ylabel(r"$S_{9\mathrm{GHz}}$ Offset (\%)", fontsize=8)
    fig.text(
        0.5, 0.04, f"Time since start of first scan (minutes)", ha="center", fontsize=8
    )
    # formating!
    # customize tick directions and lengths
    axc.tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    axc.tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    axc.tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    axc.tick_params(axis="both", which="both", labelsize="8", right=True, top=True)

    axx.tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    axx.tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    axx.tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    axc.legend(loc="lower center", fontsize=8)
    axx.legend(loc="lower center", fontsize=8)

    # axc.spines["right"].set_visible(False)
    # axc[1].spines["left"].set_visible(False)
    # axc[0].yaxis.set_ticks_position("left")
    # axc[1].yaxis.set_ticks_position("right")
    # axx[0].spines["right"].set_visible(False)
    # axx[1].spines["left"].set_visible(False)
    # axx[0].yaxis.set_ticks_position("left")
    # axx[1].yaxis.set_ticks_position("right")

    axc.xaxis.set_ticklabels([])

    # axc.set_xlim([-2.5, 27.5])
    # axx.set_xlim([-2.5, 27.5])
    # d = 1

    # kwargs = dict(
    #     marker=[(-1, -d), (1, d)],
    #     markersize=12,
    #     linestyle="none",
    #     color="k",
    #     mec="k",
    #     mew=1,
    #     clip_on=False,
    # )

    # axc[0].plot([1, 1], [0, 1], transform=axc[0].transAxes, **kwargs)
    # axc[1].plot([0, 0], [0, 1], transform=axc[1].transAxes, **kwargs)

    # axx[0].plot([1, 1], [0, 1], transform=axx[0].transAxes, **kwargs)
    # axx[1].plot([0, 0], [0, 1], transform=axx[1].transAxes, **kwargs)

    axc.axhline(y=0, color="k", alpha=0.5, linestyle="--")
    axx.axhline(y=0, color="k", alpha=0.5, linestyle="--")
    plt.tight_layout()
    print("saving figure")
    plt.savefig(f"{save_dir}/lightcurve_j215436.{ext}", overwrite=True, bbox_inches='tight')
    plt.clf()
    plt.close()
    return


def plt_lightcurve(
    save_dir,
    scan_times_src1,
    scan_times_src2,
    src1_fluxes,
    src2_fluxes,
    err_src1_fluxes,
    err_src2_fluxes,
    mod_src1,
    mod_src2,
    ext="pdf",
    src1_nm="J001513",
    src2_nm="J020507",
):
    src1_fluxes_c = src1_fluxes[0]
    src1_fluxes_x = src1_fluxes[1]
    err_src1_fluxes_c = err_src1_fluxes[0]
    err_src1_fluxes_x = err_src1_fluxes[1]

    src2_fluxes_c = src2_fluxes[0]
    src2_fluxes_x = src2_fluxes[1]
    err_src2_fluxes_c = err_src2_fluxes[0]
    err_src2_fluxes_x = err_src2_fluxes[1]

    mod_src1_c = mod_src1[0]
    mod_src1_x = mod_src1[1]
    mod_src2_c = mod_src2[0]
    mod_src2_x = mod_src2[1]

    fig = plt.figure(figsize=(3.5, 2.2))
    gs = fig.add_gridspec(2, 2, hspace=0, wspace=0.05)
    ax = gs.subplots(sharey=True)
    fig.suptitle("Light Curves October 2021", fontsize=8)

    axc = ax[0]
    axx = ax[1]
    axc[0].errorbar(
        scan_times_src1[0],
        src1_fluxes_c[0:20],
        yerr=err_src1_fluxes_c[0:20],
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_c}",
        linestyle="None",
        marker=".",
        elinewidth=0.5,
        markersize=2,
    )

    axc[1].errorbar(
        scan_times_src1[1],
        src1_fluxes_c[20:40],
        yerr=err_src1_fluxes_c[20:40],
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_c}",
        linestyle="None",
        marker=".",
        elinewidth=0.5,
        markersize=2,
    )

    axx[0].errorbar(
        scan_times_src1[0],
        src1_fluxes_x[0:20],
        yerr=err_src1_fluxes_x[20:40],
        color="C4",
        elinewidth=0.5,
        label=f"{src1_nm}, $m=$ {mod_src1_x}",
        linestyle="None",
        marker=".",
        markersize=2,
    )

    axx[1].errorbar(
        scan_times_src1[1],
        src1_fluxes_x[20:40],
        yerr=err_src1_fluxes_x[0:20],
        color="C4",
        elinewidth=0.5,
        label=f"{src1_nm}, $m=$ {mod_src1_x}",
        linestyle="None",
        marker=".",
        markersize=2,
    )

    axc[0].errorbar(
        scan_times_src2[0],
        src2_fluxes_c[0:20],
        yerr=err_src2_fluxes_c[0:20],
        color="C6",
        elinewidth=0.5,
        label=f"{src2_nm}, $m=$ {mod_src2_c}",
        linestyle="None",
        marker=".",
        markersize=2,
    )

    axc[1].errorbar(
        scan_times_src2[1],
        src2_fluxes_c[20:40],
        yerr=err_src2_fluxes_c[20:40],
        color="C6",
        elinewidth=0.5,
        label=f"{src2_nm}, $m=$ {mod_src2_c}",
        linestyle="None",
        marker=".",
        markersize=2,
    )

    axx[0].errorbar(
        scan_times_src2[0],
        src2_fluxes_x[0:20],
        yerr=err_src2_fluxes_x[20:40],
        color="C6",
        elinewidth=0.5,
        label=f"{src2_nm}, $m=$ {mod_src2_x}",
        linestyle="None",
        marker=".",
        markersize=2,
    )

    axx[1].errorbar(
        scan_times_src2[1],
        src2_fluxes_x[20:40],
        yerr=err_src2_fluxes_x[0:20],
        color="C6",
        elinewidth=0.5,
        label=f"{src2_nm}, $m=$ {mod_src2_x}",
        linestyle="None",
        marker=".",
        markersize=2,
    )

    axc[0].set_ylabel(r"$S_{5.5\mathrm{GHz}}$ Offset (\%)", fontsize=6)
    axx[0].set_ylabel(r"$S_{9\mathrm{GHz}}$ Offset (\%)", fontsize=6)
    fig.text(
        0.5, -0.04, f"Time since start of first scan (minutes)", ha="center", fontsize=6
    )
    # formating!
    # customize tick directions and lengths
    axc[0].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    axc[0].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    axc[1].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    axc[1].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    axc[0].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    axc[1].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)

    axx[0].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    axx[0].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    axx[1].tick_params(
        axis="both", which="major", direction="in", length=2.5, width=0.5, pad=5
    )
    axx[1].tick_params(axis="both", which="minor", direction="in", length=2, width=0.5)
    axx[0].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)
    axx[1].tick_params(axis="both", which="both", labelsize="8", right=True, top=True)

    axc[0].legend(loc="lower right", fontsize=6)
    axx[0].legend(loc="lower right", fontsize=6)

    axc[0].spines["right"].set_visible(False)
    axc[1].spines["left"].set_visible(False)
    axc[0].yaxis.set_ticks_position("left")
    axc[1].yaxis.set_ticks_position("right")
    axx[0].spines["right"].set_visible(False)
    axx[1].spines["left"].set_visible(False)
    axx[0].yaxis.set_ticks_position("left")
    axx[1].yaxis.set_ticks_position("right")

    axc[0].xaxis.set_ticklabels([])
    axc[1].xaxis.set_ticklabels([])

    axc[0].set_xlim([-2.5, 27.5])
    axx[0].set_xlim([-2.5, 27.5])

    axc[1].set_xlim([98.5 - 2.5, 98.5 + 27.5])
    axx[1].set_xlim([98.5 - 2.5, 98.5 + 27.5])
    d = 1

    kwargs = dict(
        marker=[(-1, -d), (1, d)],
        markersize=5,
        linestyle="none",
        color="k",
        mec="k",
        mew=1,
        clip_on=False,
    )

    axc[0].plot([1, 1], [0, 1], transform=axc[0].transAxes, **kwargs)
    axc[1].plot([0, 0], [0, 1], transform=axc[1].transAxes, **kwargs)

    axx[0].plot([1, 1], [0, 1], transform=axx[0].transAxes, **kwargs)
    axx[1].plot([0, 0], [0, 1], transform=axx[1].transAxes, **kwargs)

    axc[0].axhline(y=0, color="k", alpha=0.5, linestyle="--", linewidth=0.5)
    axc[1].axhline(y=0, color="k", alpha=0.5, linestyle="--", linewidth=0.5)
    axx[0].axhline(y=0, color="k", alpha=0.5, linestyle="--", linewidth=0.5)
    axx[1].axhline(y=0, color="k", alpha=0.5, linestyle="--", linewidth=0.5)
    # plt.tight_layout()
    plt.savefig(f"{save_dir}/lightcurve.{ext}", overwrite=True, bbox_inches='tight')
    plt.clf()
    plt.close()
    return


def plt_mwa_sed(
    data_dir,
    save_dir,
    gleam_target,
    model,
    ext="pdf",
    epochs=[
        "2013",
        "2014",
        "2020-01",
        "2020-03",
        "2020-04",
        "2020-05",
        "2020-07",
        "2020-09",
    ],
    colors=cmr.take_cmap_colors(
        "cmr.gothic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
    epoch_nms=["2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "Jul20", "Sept20"],
):
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
        ]
    )

    mwa_fluxes, err_mwa_fluxes = fitfuncts.read_mwa_fluxes(
        "/data/MWA", target, gleam_target, epochs
    )
    apr_normfact = np.nanmax(mwa_fluxes[4])
    may_normfact = np.nanmedian(mwa_fluxes[4][10:20] - mwa_fluxes[5][10:20])
    jul_normfact = np.nanmedian(mwa_fluxes[4][10:20] - mwa_fluxes[6][10:20])
    oct_normfact = np.nanmedian(mwa_fluxes[4][10:20] - mwa_fluxes[7][10:20])
    yvals, nu_p, err_nu_p, alpha, err_alpha = calc_modelnparams(data_dir, target, model)
    if target == "J223933":
        ylabel = "Relative Flux Density"
        yunit = ""
        yvals = yvals[4:8]
        colors = colors[4:8]
        epochnms = epoch_nms[4:8]
        mwaflux_2020 = mwa_fluxes[4:8]
        errmwaflux_2020 = err_mwa_fluxes[4:8]
        errmwaflux_2020 = errmwaflux_2020 / apr_normfact
        mwaflux_2020[0] = mwaflux_2020[0] / apr_normfact
        mwaflux_2020[1] = (mwaflux_2020[1] + may_normfact) / apr_normfact
        mwaflux_2020[2] = (mwaflux_2020[2] + jul_normfact) / apr_normfact
        mwaflux_2020[3] = (mwaflux_2020[3] + oct_normfact) / apr_normfact
        yvals[0] = yvals[0] / apr_normfact
        yvals[1] = (yvals[1] + may_normfact) / apr_normfact
        yvals[2] = (yvals[2] + jul_normfact) / apr_normfact
        yvals[3] = (yvals[3] + oct_normfact) / apr_normfact
        yvals_2020 = yvals
    else:
        ylabel = "Flux Density"
        yunit = "Jy"
        yvals_2020 = yvals
        epochnms = epoch_nms
        mwaflux_2020 = mwa_fluxes
        errmwaflux_2020 = err_mwa_fluxes
    f = CF.sed_fig()
    for i in range(len(epochnms)):
        if epochnms[i] in ["2013", "2014", "Apr20", "May20", "Jul20", "Sept20"]:
            f.plot_spectrum(
                frequency,
                mwaflux_2020[i],
                errmwaflux_2020[i],
                marker="o",
                label=epochnms[i],
                marker_color=colors[i],
                s=5,
            )
            f.display_model(
                np.linspace(0.01, 25, num=10000), yvals_2020[i], colors[i], lw=1
            )
    f.legend(loc="lower right")#, fontsize=20)
    # f.format(xunit="GHz",xlabelfontsize=25, ylabelfontsize=25, ticklabelfontsize=20,majorticklength=7.5, minorticklength=6, tickwidth=1, ylabel=ylabel, yunit=yunit)
    f.format(xunit="GHz", ylabel=ylabel, yunit=yunit)
    f.title(f"{gleam_target}")#, fontsize=40)
    f.save(f"{save_dir}{target}_mwa", ext=ext)
    plt.clf()
    plt.close()
    return


def plt_modelsonly(
    save_dir,
    gleam_target,
    fit_freq,
    src_flux,
    err_src_flux,
    yvals,
    logz,
    epoch_nm,
    ext="pdf",
    colors=cmr.take_cmap_colors(
        "cmr.cosmic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
):

    model_nms = [
        "SSA",
        "FFA",
        "inFFA",
        "SSAb",
        "FFAb",
        "inFFAb",
    ]
    target = gleam_target.strip("GLEAM ")[0:7]
    f = CF.sed_fig()
    f.plot_spectrum(
        fit_freq,
        src_flux,
        err_src_flux,
        marker="o",
        marker_color="k",
        s=5,
    )
    for i in range((len(yvals))):
        f.display_model(
            np.linspace(0.01, 25, num=10000),
            yvals[i],
            colors[i],
            lw=0.5,
            label=f"{model_nms[i]}, logz: {logz[i]:.2f}",
        )
    f.legend(loc="lower left", ncol=1)
    f.format(xunit="GHz")
    f.title(f"{gleam_target}: {epoch_nm} Models")
    f.save(f"{save_dir}{target}_{epoch_nm}_models", ext=ext)
    plt.close()
    plt.clf()


def plt_sed(
    data_dir,
    save_dir,
    gleam_target,
    src_flux,
    err_src_flux,
    extra_surveys,
    yvals,
    ext="pdf",
    colors=cmr.take_cmap_colors(
        "cmr.gothic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
    epochnms=[
        "2013",
        "2014",
        "Jan20",
        "Mar20",
        "Apr20",
        "May20",
        "Jul20",
        "Sept20",
        "2020",
    ],
    frequency=np.array(
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
    ),
):
    extra_frequencies = [
        0.150,
        0.408,
        0.843,
        1.400,
        0.074,
        20,
        8.6,
        4.8,
        0.8875,
    ]
    target = gleam_target.strip("GLEAM ")[0:7]
    f = CF.sed_fig()
    extra_markers = ["X", "s", "*", "p", "D", ">", "<", "^", "1"]
    print(np.shape(yvals))
    print(np.shape(src_flux))
    for i in range((len(yvals))):
        f.display_model(np.linspace(0.01, 25, num=10000), yvals[i], colors[i])
    for i in range(len(extra_surveys)):
        if extra_surveys[i] == 0.0:
            extra_surveys[i] = np.nan
        f.plot_point(
            extra_frequencies[i],
            extra_surveys[i],
            marker=extra_markers[i],
            marker_color="k",
            s=10,
            alpha=0.5,
        )
    for i in range((8)):
        f.plot_spectrum(
            frequency,
            src_flux[i],
            err_src_flux[i],
            marker="o",
            label=epochnms[i],
            marker_color=colors[i],
            s=5,
        )

    # f.plot_spectrum(frequency, src_flux[-1], err_src_flux[-1], marker="X", label="2020", marker_color="k",s=5)
    f.legend(loc="lower center")#, fontsize=20)
    f.format(xunit="GHz")#,xlabelfontsize=25, ylabelfontsize=25, ticklabelfontsize=20,majorticklength=7.5, minorticklength=6, tickwidth=1)
    f.title(f"{gleam_target}")#, fontsize=40)
    f.save(f"{save_dir}{target}_sed", ext=ext)
    plt.close()
    plt.clf()
    return


def run_everything(
    save_dir,
    data_dir,
    gleam_tar,
    epochs=[
        "2013",
        "2014",
        "2020-01",
        "2020-03",
        "2020-04",
        "2020-05",
        "2020-07",
        "2020-09",
        "2020",
    ],
):
    target = gleam_tar.strip("GLEAM ")[0:7]
    fit_models = [
        "SSA",
        "FFA",
        "inFFA",
        "SSAb",
        "FFAb",
        "inFFAb",
    ]
    fit_flux, err_fit_flux, fit_freq = fitfuncts.createfitflux(data_dir, gleam_tar)
    print(np.shape(fit_flux))
    src_flux, err_src_flux = fitfuncts.createsrcflux(data_dir, gleam_tar)
    # print(np.shape(src_flux))
    diff = src_flux[0][4:20] - src_flux[1][4:20]
    diff_err = np.sqrt((err_src_flux[0][4:20])**2+(err_src_flux[1][4:20])**2)
    var_param = (np.sum(((src_flux[0][4:20]-src_flux[1][4:20])**2)/(diff_err)**2))
    print(var_param)
    extra_fluxes = fitfuncts.read_extra_fluxes("/data/MWA", gleam_tar)
    # moss = np.zeros(len(flux_yr1[0]))
    # for i in range(len(flux_yr1[0])):
    mean_resid = np.mean(diff)
    diff_median = np.median(diff)
    moss = np.sum((diff_median-diff)**2/(diff_err)**2)
    print(moss)
    extra_fluxes = fitfuncts.read_extra_fluxes("/data/MWA", gleam_tar)
    print(f"Running for {target}")
    for i in range(len(fit_flux)):
        if epoch_nms[i] == "2020-04" and target == "J020507":
            print(f"{target} {epoch_nms[i]}, skipping .... ")
        elif epoch_nms[i] == "2020-04" and target == "J024838":
            print(f"{target} {epoch_nms[i]}, skipping .... ")
        else:
            for j in range(len(fit_models)):
                model = fit_models[j]
                model_funct = model_params_dict[model][0]
                model_trans = model_params_dict[model][1]
                labels = model_params_dict[model][2:]
                try:
                    sampler = open(
                        f"{save_dir}{target}/{epoch_nms[i]}/{model}/run1/info/results.json"
                    )
                    # print("Found results of run, continuing with analysis")
                except FileNotFoundError:
                    sampler = fitfuncts.run_ultranest_mcmc(
                        f"{save_dir}{target}/{epoch_nms[i]}/{model}",
                        labels,
                        fit_freq[i],
                        fit_flux[i],
                        err_fit_flux[i],
                        model_funct,
                        model_trans,
                        run_num=1,
                    )
                    sampler.run(max_iters=50000, show_status=False, viz_callback=False)
                    print(f"Finished fitting for {model} {epoch_nms[i]}")
                    sampler.store_tree()
                    sampler = open(
                        f"{save_dir}{target}/{epoch_nms[i]}/{model}/run1/info/results.json"
                    )
        yvals_epoch, logz_epoch = calc_yvals(
            f"{save_dir}/{target}", target, epoch_nms[i]
        )
        plt_modelsonly(
            f"{save_dir}",
            gleam_tar,
            fit_freq[i],
            fit_flux[i],
            err_fit_flux[i],
            yvals_epoch,
            logz_epoch,
            epoch_nms[i],
            ext="pdf",
        )

    # Plotting section
    avg_logz = calc_logzmod(target)
    calc_secbayes(avg_logz)
    model_nm = fit_models[np.argmax(avg_logz)]
    print(f"Best model: {model_nm}")
    model_funct = model_params_dict[model_nm][0]
    model = [model_nm, model_funct]
    yvals, nu_p, err_nu_p, alpha, err_alpha = calc_modelnparams(
        f"{save_dir}/{target}", target, model
    )
    plt.clf()
    plt.close()
    if target in ["J015445", "J020507", "J024838", "J223933", "J215436"]:
        plt_mwa_sed(f"{save_dir}/{target}", f"{save_dir}Plots/", gleam_tar, model)
        plt.clf()
        plt.close()

    plt_sed(
        data_dir,
        f"{save_dir}Plots/",
        gleam_tar,
        src_flux,
        err_src_flux,
        extra_fluxes,
        yvals,
        ext="pdf",
    )
    if target in ["J015445", "J020507", "J024838"]:
        plt_peakftime(
            f"{save_dir}Plots/{target}",
            nu_p,
            err_nu_p,
            f"{gleam_tar} Peak Frequency",
            # ext="pdf"
        )
    plt_alphatime(
        f"{save_dir}{target}",
        alpha,
        err_alpha,
        f"{gleam_tar} Spectral Index",
        ext="pdf",
    )
    if target == "J020507":
        yvals_oct, logz_oct = calc_yvals(f"{save_dir}/{target}", target, "2020-09")
        oct_flux = fit_flux[-1]
        err_oct_flux = err_fit_flux[-1]
        # oct_flux2 = src_flux[4][20:37]
        # err_oct_flux2 = err_src_flux[4][20:37]

        # oct_flux = np.concatenate((oct_flux1, oct_flux2))
        # err_oct_flux = np.concatenate((err_oct_flux1, err_oct_flux2))
        plt_modelsonly(
            f"{save_dir}Plots/",
            gleam_tar,
            fit_freq[0],
            oct_flux,
            err_oct_flux,
            yvals_oct,
            logz_oct,
            "Sept20",
        )
    return avg_logz
