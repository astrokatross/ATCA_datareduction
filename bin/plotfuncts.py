#!/usr/bin/python
# This script is to plot and analyse the ATCA data reduction
# By K.Ross 20/1/21


import matplotlib.pyplot as plt

plt.rcParams["axes.grid"] = False
plt.rcParams["font.family"] = "serif"


def plt_sn(directory, sn, self_cal_round):
    fig = plt.figure(1, figsize=(15, 10))
    gs = plt.GridSpec(1, 1)
    ax = fig.subplot(gs[0])

    ax.scatter(self_cal_round, sn, color="k")
    plt.xlabel("Self Cal Round", fontsize=20)
    plt.ylabel("local_rms", fontsize=20)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)

    ax.tick_params(
        axis="both", which="both", direction="in", labelsize=20, top="on", right="on"
    )
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=8,
        width=1.5,
        top="on",
        right="on",
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=5,
        width=1.5,
        top="on",
        right="on",
    )
    plt.savefig(directory, bbox_inches="tight")
    plt.clf()
