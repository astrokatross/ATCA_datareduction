#!/usr/bin/python3
# This script plots SEDs for sources after. Assumes run_process.py has run for
# each epoch/band successfully

# By K.Ross 25/5/21

# Importing relevant packages
import plotfuncts
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False

directory = str(os.environ["PROJECT"])
tarname = str(os.environ["TARGET"])
calibrator = str(os.environ["CALIBRATOR"])

plotfuncts.plt_sed(directory, tarname)
plotfuncts.plt_secondary(directory, tarname)
plotfuncts.plt_nearby(directory, tarname)

# if calibrator == "PRI":
#     plotfuncts.plt_secondary(directory, tarname)
