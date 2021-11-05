#!/usr/bin/python3
# This script has functions for analysising the variability and returning useful and interpretable parameters
# By K.Ross 29/09/21

import numpy as np
import gpscssmodels
import json
import cmasher as cmr
import matplotlib.pyplot as plt


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
        vp_temp = np.squeeze(np.array(parameters[paramnames == 'freqpeak']))
        S_v_temp.append(calc_S_v(model, parameters, vp_temp))
        vp.append(vp_temp)
    S_v = list(zip(*S_v_temp))
    delta_S150 = max(S_v[0])-min(S_v[0])
    delta_S2 = max(S_v[1])-min(S_v[1])
    delta_S5 = max(S_v[2])-min(S_v[2])
    delta_S9 = max(S_v[3])-min(S_v[3])
    delta_Svp = max(S_v[4])-min(S_v[4])
    colors = cmr.take_cmap_colors(
        "cmr.bubblegum", len(S_v[0]), cmap_range=(0.15, 0.85), return_fmt="hex"
    )
    num_epochs = np.array([-10, -5, 1, 4, 5, 7, 10])
    fig = plt.figure(1, figsize=(40,18))
    gs = fig.add_gridspec(5, hspace=0)
    axes = gs.subplots(sharex=True)
    # fig, axes = plt.subplots(5, 1)
    fig.suptitle("Flux Density With Time",fontsize=30)
    axes[0].set_title(f"$\Delta$S\_150MHz={delta_S150:.3f}Jy, $\Delta$S\_2.1GHz={delta_S2:.3f}Jy, $\Delta$S\_5.5GHz={delta_S5:.3f}Jy, $\Delta$S\_9.5GHz={delta_S9:.3f}, $\Delta$S\_vp={delta_Svp:.3f}",fontsize=20)
    for i in range(len(S_v[0])):
        axes[0].scatter(num_epochs[i], S_v[0][i], color=colors[i])
        axes[1].scatter(num_epochs[i], S_v[1][i], color=colors[i])
        axes[2].scatter(num_epochs[i], S_v[2][i], color=colors[i])
        axes[3].scatter(num_epochs[i], S_v[3][i], color=colors[i])
        axes[4].scatter(num_epochs[i], S_v[4][i], color=colors[i], label=f"vp={vp[i]:.2f}GHz")
    axes[0].set_ylabel("S\_150MHz (Jy)",fontsize=20)
    axes[1].set_ylabel("S\_2.1GHz (Jy)",fontsize=20)
    axes[2].set_ylabel("S\_5.5GHz (Jy)",fontsize=20)
    axes[3].set_ylabel("S\_9.5GHz (Jy)",fontsize=20)
    axes[4].set_ylabel("S\_vp (Jy)",fontsize=20)

    axes[0].set_ylim([0.9*np.min(S_v[0]),1.1*np.max(S_v[0])])
    axes[1].set_ylim([0.9*np.min(S_v[1]),1.1*np.max(S_v[1])])
    axes[2].set_ylim([0.9*np.min(S_v[2]),1.1*np.max(S_v[2])])
    axes[3].set_ylim([0.9*np.min(S_v[3]),1.1*np.max(S_v[3])])
    axes[4].set_ylim([0.9*np.min(S_v[4]),1.1*np.max(S_v[4])])

    axes[4].legend(loc='upper center', bbox_to_anchor=(0.48,-0.2), ncol=7, fontsize=14)

    axes[4].set_xticks(num_epochs)
    axes[4].set_xticklabels(epochs,fontsize=20)
    axes[0].tick_params(axis='y', labelsize=20)
    axes[1].tick_params(axis='y', labelsize=20)
    axes[2].tick_params(axis='y', labelsize=20)
    axes[3].tick_params(axis='y', labelsize=20)
    axes[4].tick_params(axis='y', labelsize=20)
    plt.savefig(f"{directory}/{target}_Sv_vs_time.png", overwrite=True)
    return
