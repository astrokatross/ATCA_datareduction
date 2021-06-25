#!/usr/bin/python3
# This script plots SEDs for sources after. Assumes run_process.py has run for each epoch/band successfully

# By K.Ross 25/5/21

# Importing relevant packages
import matplotlib.pyplot as plt
import plotfuncts

plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False

# constants:
tarname = "J015445"
directory = "/data/var_analysis/ATCA/Code/"

plotfuncts.plt_sed(directory, tarname, models=True, extra_surveys=True)
