#!/usr/bin/python3
# This script plots SEDs for sources after. Assumes run_process.py has run for
# each epoch/band successfully

# By K.Ross 25/5/21

# Importing relevant packages
# import plotfuncts
import matplotlib.pyplot as plt
import os
import numpy as np
import gpscssmodels
import cmasher as cmr
import plotfuncts

plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False



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
nwalkers = 50
ndim = 6
save_dir = "/data/ATCA/analysis/"
data_dir = "data/ATCA/ATCA_datareduction/"
gleam_tar = "GLEAM J020507-110922"
tar = "J020507"
chosen_model = gpscssmodels.singhomobremss
epoch = 3
epoch_nm = epoch_nms[epoch]
model_nm = "singhomobremss"
labels = ["Snorm", "alpha", "freqpeak"]


# directory = str(os.environ["PROJECT"])
# tarname = str(os.environ["TARGET"])
# calibrator = str(os.environ["CALIBRATOR"])

# plotfuncts.plt_sed(directory, tarname)
# plotfuncts.plt_secondary(directory, tarname)
# plotfuncts.plt_nearby(directory, tarname)

# if calibrator == "PRI":
#     plotfuncts.plt_secondary(directory, tarname)
