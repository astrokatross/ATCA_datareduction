#!/usr/bin/python3
# This script plots SEDs for sources after. Assumes run_process.py has run for
# each epoch/band successfully

# By K.Ross 25/5/21

# Importing relevant packages
import plotfuncts
import matplotlib.pyplot as plt


plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False

# constants:
tarname = "J032213"
directory = "/data/ATCA/ATCA_datareduction/"

plotfuncts.plt_sed(directory, tarname, models=True, extra_surveys=True)
