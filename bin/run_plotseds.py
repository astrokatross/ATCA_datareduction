#!/usr/bin/python
# This script plots SEDs for sources after. Assumes run_process.py has run for each epoch/band successfully

# By K.Ross 25/5/21

# Importing relevant packages
import matplotlib.pyplot as plt
import plotfuncts

plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False

# constants:
tarname = "J001513"
directory = "/data/var_analysis/ATCA/Code/"

plotfuncts.plt_sed(directory, tarname, epoch3=False)
