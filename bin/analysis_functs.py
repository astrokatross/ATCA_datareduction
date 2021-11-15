#!/usr/bin/python3
# This script has functions for analysising the variability and returning useful and interpretable parameters
# By K.Ross 29/09/21

import numpy as np
import gpscssmodels
import json
import cmasher as cmr
import matplotlib.pyplot as plt
import math
from casatasks import uvmodelfit
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
    "2020-01",
    "2020-03",
    "2020-04",
    "2020-05",
    "2020-07",
    "2020-10",
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


def calc_S_v(model, params, vp):
    # Firstly for 150MHz:
    S_150mhz = model(0.15, *params)

    # Calculating at 2.1GHz
    S_2Ghz = model(2.1, *params)

    # 5.5GHz
    S_5Ghz = model(5.5, *params)

    # 9.5GHz
    S_9Ghz = model(9.5, *params)

    # At Peak Frequency caluclated from fitting previously
    S_vp = model(vp, *params)
    return np.array([S_150mhz, S_2Ghz, S_5Ghz, S_9Ghz, S_vp])


def plot_S_v(directory, target, epochs, model_nm, model):
    S_v_temp = []
    vp = []
    for epoch in epochs:
        sampler = open(
            f"/data/ATCA/analysis/{target}/{epoch}/{model_nm}/run1/info/results.json"
        )
        results = json.load(sampler)
        paramnames = np.array(results["paramnames"])
        parameters = np.array(results["maximum_likelihood"]["point"])
        vp_temp = np.squeeze(np.array(parameters[paramnames == "freqpeak"]))
        S_v_temp.append(calc_S_v(model, parameters, vp_temp))
        vp.append(vp_temp)
    S_v = list(zip(*S_v_temp))
    delta_S150 = max(S_v[0]) - min(S_v[0])
    delta_S2 = max(S_v[1]) - min(S_v[1])
    delta_S5 = max(S_v[2]) - min(S_v[2])
    delta_S9 = max(S_v[3]) - min(S_v[3])
    delta_Svp = max(S_v[4]) - min(S_v[4])
    colors = cmr.take_cmap_colors(
        "cmr.bubblegum", len(S_v[0]), cmap_range=(0.15, 0.85), return_fmt="hex"
    )
    num_epochs = np.array([-10, -5, 1, 4, 5, 7, 10])
    fig = plt.figure(1, figsize=(40, 18))
    gs = fig.add_gridspec(5, hspace=0)
    axes = gs.subplots(sharex=True)
    # fig, axes = plt.subplots(5, 1)
    fig.suptitle("Flux Density With Time", fontsize=30)
    axes[0].set_title(
        f"$\Delta$S\_150MHz={delta_S150:.3f}Jy, $\Delta$S\_2.1GHz={delta_S2:.3f}Jy, $\Delta$S\_5.5GHz={delta_S5:.3f}Jy, $\Delta$S\_9.5GHz={delta_S9:.3f}, $\Delta$S\_vp={delta_Svp:.3f}",
        fontsize=20,
    )
    for i in range(len(S_v[0])):
        axes[0].scatter(num_epochs[i], S_v[0][i], color=colors[i])
        axes[1].scatter(num_epochs[i], S_v[1][i], color=colors[i])
        axes[2].scatter(num_epochs[i], S_v[2][i], color=colors[i])
        axes[3].scatter(num_epochs[i], S_v[3][i], color=colors[i])
        axes[4].scatter(
            num_epochs[i], S_v[4][i], color=colors[i], label=f"vp={vp[i]:.2f}GHz"
        )
    axes[0].set_ylabel("S\_150MHz (Jy)", fontsize=20)
    axes[1].set_ylabel("S\_2.1GHz (Jy)", fontsize=20)
    axes[2].set_ylabel("S\_5.5GHz (Jy)", fontsize=20)
    axes[3].set_ylabel("S\_9.5GHz (Jy)", fontsize=20)
    axes[4].set_ylabel("S\_vp (Jy)", fontsize=20)

    axes[0].set_ylim([0.9 * np.min(S_v[0]), 1.1 * np.max(S_v[0])])
    axes[1].set_ylim([0.9 * np.min(S_v[1]), 1.1 * np.max(S_v[1])])
    axes[2].set_ylim([0.9 * np.min(S_v[2]), 1.1 * np.max(S_v[2])])
    axes[3].set_ylim([0.9 * np.min(S_v[3]), 1.1 * np.max(S_v[3])])
    axes[4].set_ylim([0.9 * np.min(S_v[4]), 1.1 * np.max(S_v[4])])

    axes[4].legend(loc="upper center", bbox_to_anchor=(0.48, -0.2), ncol=7, fontsize=14)

    axes[4].set_xticks(num_epochs)
    axes[4].set_xticklabels(epochs, fontsize=20)
    axes[0].tick_params(axis="y", labelsize=20)
    axes[1].tick_params(axis="y", labelsize=20)
    axes[2].tick_params(axis="y", labelsize=20)
    axes[3].tick_params(axis="y", labelsize=20)
    axes[4].tick_params(axis="y", labelsize=20)
    plt.savefig(f"{directory}/{target}_Sv_vs_time.png", overwrite=True)
    return


def calc_modelnparams(directory, target, model):
    freq_cont = np.linspace(0.01, 25, num=10000)
    yvals = []
    nu_p = []
    errlo_nu_p = []
    errup_nu_p = []
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
        "2020-10",
    ]
    for i in range(len(epochs)):
        try:
            sampler = open(f"{directory}/{epochs[i]}/{model_nm}/run2/info/results.json")
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
        except FileNotFoundError:
            yvals.append(np.full(len(freq_cont), np.nan))
            nu_p.append(np.nan)
            errlo_nu_p.append(np.nan)
            errup_nu_p.append(np.nan)
            print("Can't find sampler....")
        err_nu_p = np.stack((errlo_nu_p, errup_nu_p))
    return yvals, nu_p, err_nu_p


def read_results(dir):
    try:
        sampler = open(f"{dir}/run2/info/results.json")
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
        "2020-01",
        "2020-03",
        "2020-04",
        "2020-05",
        "2020-07",
        "2020-10",
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
        avg_logz = sum_logz / 7
    else:
        avg_logz = sum_logz / 8
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
            new_time = (new_time + datetime.timedelta(0,30))
            current_str = new_time.strftime("%H:%M:%S")
            timeranges.append(current_str)
    return timeranges


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
    figsize = (20, 10)
    # fig = plt.figure(1, figsize=figsize, facecolor="white")
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, wspace=0.05, width_ratios=[1, 3])
    ax = gs.subplots()
    fig.suptitle(plot_title, fontsize=40)
    ax[0].set_ylabel(r"$\nu_p$", fontsize=30)
    ax[0].tick_params(
        axis="both", which="major", direction="in", length=6, width=1.5, pad=5
    )
    ax[0].tick_params(axis="both", which="minor", direction="in", length=4, width=1.5)
    ax[0].tick_params(axis="both", which="both", labelsize="25", right=True, top=True)
    ax[1].tick_params(
        axis="both", which="major", direction="in", length=6, width=1.5, pad=5
    )
    ax[1].tick_params(axis="both", which="minor", direction="in", length=4, width=1.5)
    ax[1].tick_params(axis="both", which="both", labelsize="25", right=True, top=True)
    ax[0].spines["right"].set_visible(False)
    ax[1].spines["left"].set_visible(False)
    ax[0].yaxis.set_ticks_position("left")
    ax[1].yaxis.set_ticks_position("right")
    ax[1].yaxis.set_ticklabels([])
    d = 1
    kwargs = dict(
        marker=[(-1, -d), (1, d)],
        markersize=12,
        linestyle="none",
        color="k",
        mec="k",
        mew=1,
        clip_on=False,
    )
    ax[0].plot([1, 1], [0, 1], transform=ax[0].transAxes, **kwargs)
    ax[1].plot([0, 0], [0, 1], transform=ax[1].transAxes, **kwargs)

    months = [-5, -4, 1, 3, 4, 5, 7, 10]
    for i in range(len(months)):
        ax[0].errorbar(
            months[i],
            nu_p[i],
            yerr=np.array([[err_nu_p[0, i], err_nu_p[1, i]]]).T,
            fmt="o",
            color=colors[i],
            capsize=3,
            markersize=10,
        )
        ax[1].errorbar(
            months[i],
            nu_p[i],
            yerr=np.array([[err_nu_p[0, i], err_nu_p[1, i]]]).T,
            fmt="o",
            color=colors[i],
            capsize=3,
            markersize=10,
        )
    for axis in ["top", "bottom", "left", "right"]:
        ax[0].spines[axis].set_linewidth(2)
        ax[1].spines[axis].set_linewidth(2)
    ax[0].set_xticks(months)
    ax[1].set_xticks(months)
    ax[0].set_xticklabels(
        ["2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20"]
    )
    ax[1].set_xticklabels(
        ["2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20"]
    )
    ax[0].set_xlim([-6, -3])
    ax[1].set_xlim([0, 12])
    plt.tight_layout()
    
    plt.savefig(f"{save_dir}_peakfreqvstime.{ext}", overwrite=True)
    plt.clf()
    plt.close()
    return


def read_lightcurveflux(data_dir, outfile_dir, timeranges):
    fluxes = []
    err_fluxes = []
    for i in range(len(timeranges)):
        timerange = f"{timeranges[i]}+00:00:30"
        outfile = f"{outfile_dir}_{timerange}.cl"
        tar_ms = f"{data_dir}_selfcaltime.ms"
        print(tar_ms)
        if (os.path.exists(outfile)) is False:
            uvmodelfit(
                vis=tar_ms,
                niter=10,
                comptype="P",
                outfile=outfile,
                field="0",
                selectdata=True,
                timerange=timerange,
            )
            tbl = table(outfile)
            flux = tbl.getcell("Flux", 0)[0].astype("float64")
            err = np.sqrt((flux * 0.05) ** 2 + (0.002 ** 2))
            fluxes.append(flux)
            err_fluxes.append(err)
        else:
            tbl = table(outfile)
            flux = tbl.getcell("Flux", 0)[0].astype("float64")
            err = np.sqrt((flux * 0.05) ** 2 + (0.002 ** 2))
            fluxes.append(flux)
            err_fluxes.append(err)
    mod = round(np.std(fluxes)/np.median(fluxes), -int(math.floor(math.log10(abs(np.std(fluxes)/np.median(fluxes))))))
    err_fluxes = (err_fluxes/np.median(fluxes))
    fluxes = (fluxes/np.median(fluxes)) - 1
    return fluxes, err_fluxes, mod


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
    
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2, hspace=0, wspace=0.05)
    ax = gs.subplots(sharey=True)
    fig.suptitle("Light Curves 2021-10-15", fontsize=40)

    axc = ax[0]
    axx = ax[1]
    axc[0].errorbar(
        scan_times_src1[0],
        src1_fluxes_c[0:20],
        yerr=err_src1_fluxes_c[0:20],
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_c}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axc[1].errorbar(
        scan_times_src1[1],
        src1_fluxes_c[20:40],
        yerr=err_src1_fluxes_c[20:40],
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_c}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axx[0].errorbar(
        scan_times_src1[0],
        src1_fluxes_x[0:20],
        yerr=err_src1_fluxes_x[20:40],
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_x}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axx[1].errorbar(
        scan_times_src1[1],
        src1_fluxes_x[20:40],
        yerr=err_src1_fluxes_x[0:20],
        color="C4",
        label=f"{src1_nm}, $m=$ {mod_src1_x}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axc[0].errorbar(
        scan_times_src2[0],
        src2_fluxes_c[0:20],
        yerr=err_src2_fluxes_c[0:20],
        color="C6",
        label=f"{src2_nm}, $m=$ {mod_src2_c}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axc[1].errorbar(
        scan_times_src2[1],
        src2_fluxes_c[20:40],
        yerr=err_src2_fluxes_c[20:40],
        color="C6",
        label=f"{src2_nm}, $m=$ {mod_src2_c}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axx[0].errorbar(
        scan_times_src2[0],
        src2_fluxes_x[0:20],
        yerr=err_src2_fluxes_x[20:40],
        color="C6",
        label=f"{src2_nm}, $m=$ {mod_src2_x}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axx[1].errorbar(
        scan_times_src2[1],
        src2_fluxes_x[20:40],
        yerr=err_src2_fluxes_x[0:20],
        color="C6",
        label=f"{src2_nm}, $m=$ {mod_src2_x}",
        linestyle="None",
        marker="o",
        markersize=5,
    )

    axc[0].set_ylabel(r"$S_{5.5\mathrm{GHz}}$ Offset (\%)", fontsize=30)
    axx[0].set_ylabel(r"$S_{9\mathrm{GHz}}$ Offset (\%)", fontsize=30)
    fig.text(0.5, 0.04, f"Time since start of first scan (minutes)", ha='center', fontsize=30)
    # formating!
    # customize tick directions and lengths
    axc[0].tick_params(
        axis="both", which="major", direction="in", length=6, width=1.5, pad=5
    )
    axc[0].tick_params(axis="both", which="minor", direction="in", length=4, width=1.5)
    axc[1].tick_params(
        axis="both", which="major", direction="in", length=6, width=1.5, pad=5
    )
    axc[1].tick_params(axis="both", which="minor", direction="in", length=4, width=1.5)
    axc[0].tick_params(axis="both", which="both", labelsize="20", right=True, top=True)
    axc[1].tick_params(axis="both", which="both", labelsize="20", right=True, top=True)

    axx[0].tick_params(
        axis="both", which="major", direction="in", length=6, width=1.5, pad=5
    )
    axx[0].tick_params(axis="both", which="minor", direction="in", length=4, width=1.5)
    axx[1].tick_params(
        axis="both", which="major", direction="in", length=6, width=1.5, pad=5
    )
    axx[1].tick_params(axis="both", which="minor", direction="in", length=4, width=1.5)
    axx[0].tick_params(axis="both", which="both", labelsize="20", right=True, top=True)
    axx[1].tick_params(axis="both", which="both", labelsize="20", right=True, top=True)

    axc[0].legend(loc="lower right", fontsize=20)
    axx[0].legend(loc="lower right", fontsize=20)

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
        markersize=12,
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

    axc[0].axhline(y=0, color="k", alpha=0.5, linestyle="--")
    axc[1].axhline(y=0, color="k", alpha=0.5, linestyle="--")
    axx[0].axhline(y=0, color="k", alpha=0.5, linestyle="--")
    axx[1].axhline(y=0, color="k", alpha=0.5, linestyle="--")
    plt.tight_layout()
    plt.savefig(f"{save_dir}/lightcurve.{ext}", overwrite=True)
    plt.clf()
    plt.close()
    return


def plt_mwa_sed(
    data_dir,
    save_dir,
    gleam_target,
    model,
    ext="pdf",
    colors=cmr.take_cmap_colors(
        "cmr.gothic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
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
    (
        mwa_flux_2013,
        err_mwa_flux_2013,
        mwa_flux_2014,
        err_mwa_flux_2014,
        fluxes_extra,
    ) = fitfuncts.read_gleam_fluxes("/data/MWA", gleam_target)
    mwaflux_2020 = []
    errmwaflux_2020 = []
    mwaflux_2020.append(mwa_flux_2013)
    mwaflux_2020.append(mwa_flux_2014)
    errmwaflux_2020.append(err_mwa_flux_2013)
    errmwaflux_2020.append(err_mwa_flux_2014)
    empty_buffer = np.full(20, np.nan)
    mwaflux_2020.append(empty_buffer)
    errmwaflux_2020.append(empty_buffer)
    mwaflux_2020.append(empty_buffer)
    errmwaflux_2020.append(empty_buffer)
    epoch_nms = ["2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "Jul20", "Oct20"]
    for epochs in ["epoch3", "epoch4", "epoch5", "epoch6"]:
        mwa_flux, err_mwa = fitfuncts.read_mwa_fluxes(
            "/data/MWA", target, gleam_target, epochs
        )
        mwaflux_2020.append(mwa_flux)
        errmwaflux_2020.append(err_mwa)
    apr_normfact = np.nanmax(mwaflux_2020[4])
    may_normfact = np.nanmedian(mwaflux_2020[4][10:20] - mwaflux_2020[5][10:20])
    jul_normfact = np.nanmedian(mwaflux_2020[4][10:20] - mwaflux_2020[6][10:20])
    oct_normfact = np.nanmedian(mwaflux_2020[4][10:20] - mwaflux_2020[7][10:20])
    yvals, nu_p, err_nu_p = calc_modelnparams(data_dir, target, model)
    if target == "J223933":
        ylabel = "Relative Flux Density"
        yunit = ""
        yvals = yvals[4:8]
        colors = colors[4:8]
        epochnms = epoch_nms[4:8]
        mwaflux_2020 = mwaflux_2020[4:8]
        errmwaflux_2020 = errmwaflux_2020[4:8]
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
    f = CF.sed_fig()
    for i in range(len(epochnms)):
        if epochnms[i] in ["2013", "2014", "Apr20", "May20", "Jul20", "Oct20"]:
            f.plot_spectrum(
                frequency,
                mwaflux_2020[i],
                errmwaflux_2020[i],
                marker="o",
                label=epochnms[i],
                marker_color=colors[i],
                s=60,
            )
            f.display_model(
                np.linspace(0.01, 25, num=10000), yvals_2020[i], colors[i], lw=3
            )
    f.legend(loc="lower center")
    f.format(xunit="GHz", ylabel=ylabel, yunit=yunit)
    f.title(f"{gleam_target}")
    f.save(f"{save_dir}{target}_mwa", ext=ext)
    plt.clf()
    plt.close()
    return


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
):

    epochnms = ["2013", "2014", "Jan20", "Mar20", "Apr20", "May20", "July20", "Oct20"]
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
    extra_markers = ["X", "s", "*", "p", "D", "1", "2", "3", "4"]
    for i in range(len(extra_surveys)):
        if extra_surveys[i] == 0.:
            extra_surveys[i] = np.nan
        f.plot_point(
            extra_frequencies[i],
            extra_surveys[i],
            marker=extra_markers[i],
            marker_color='k',
            s=30,
        )
    for i in range((8)):
        f.plot_spectrum(
            frequency,
            src_flux[i],
            err_src_flux[i],
            marker="o",
            label=epochnms[i],
            marker_color=colors[i],
            s=60,
        )
        f.display_model(np.linspace(0.01, 25, num=10000), yvals[i], colors[i], lw=1)
    f.legend(loc="lower center")
    f.format(xunit="GHz")
    f.title(f"{gleam_target}")
    f.save(f"{save_dir}{target}_sed", ext=ext)
    plt.close()
    plt.clf()
    return


def run_everything(save_dir, data_dir, gleam_tar):
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
    src_flux, err_src_flux = fitfuncts.createsrcflux(data_dir, gleam_tar)
    extra_fluxes = fitfuncts.read_extra_fluxes("/data/MWA", gleam_tar)
    for i in range(len(epoch_nms)):
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
                        f"{save_dir}{target}/{epoch_nms[i]}/{model}/run2/info/results.json"
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
                        run_num=2,
                    )
                    sampler.run(max_iters=50000)
                    print(f"Finished fitting for {model} {epoch_nms[i]}")
                    sampler.store_tree()
                    sampler = open(
                        f"{save_dir}{target}/{epoch_nms[i]}/{model}/run2/info/results.json"
                    )
                # if (
                #     os.path.exists(
                #         f"{save_dir}/{target}/seds/{target}_{epoch_nm}_{model}_sed.png"
                #     )
                #     is False
                # ):
                #     sequence, final = ultranest.integrator.read_file(
                #         f"{save_dir}{target}/{epoch_nm}/{model}/run2/",
                #         len(labels),
                #         check_insertion_order=False,
                #     )

                # band = PredictionBand(freq_cont)
                # for params in final["samples"]:
                #     band.add(model_funct(freq_cont, *params))

                # fitfuncts.plot_epochsed(
                #     f"{save_dir}/{target}/seds/{target}_{epoch_nm}_{model}",
                #     freq,
                #     src_flux,
                #     err_src_flux,
                #     band,
                #     color,
                #     target,
                #     model,
                #     epoch,
                # )

    # Plotting section
    avg_logz = calc_logzmod(target)
    calc_secbayes(avg_logz)
    model_nm = fit_models[np.argmax(avg_logz)]
    print(f"Best model: {model_nm}")
    model_funct = model_params_dict[model_nm][0]
    model = [model_nm, model_funct]
    yvals, nu_p, err_nu_p = calc_modelnparams(
        f"{save_dir}/{target}", target, model
    )
    plt.clf()
    plt.close()
    if target in ["J015445", "J020507", "J024838", "J223933"]:
        plt_mwa_sed(f"{save_dir}/{target}", f"{save_dir}Plots/", gleam_tar, model)
        plt.clf()
        plt.close()
    plt_sed(data_dir, f"{save_dir}Plots/", gleam_tar, src_flux, err_fit_flux, extra_fluxes, yvals)
    if target in ["J015445", "J020507", "J024838"]:
        plt_peakftime(
                f"{save_dir}Plots/{target}",
                nu_p,
                err_nu_p,
                f"{gleam_tar} Peak Frequency",
            )
    return nu_p, err_nu_p
