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

# constants:

source_dict = [
    # "J001513",
    # "J015445",
    # "J020507",
    # "J021246",
    # "J022744",
    # "J022802",
    # "J024838",
    # "J025040",
    # "J032213",
    # "J032836",
    # "J033023",
    # "J042502",
    # "J044033",
    # "J044737",
    # "J052824",
    # "J223933",
    "J224408",
    # "J215436",
    # "J020510",
    # "J001530",
    # "J032119",
    # "J032843",
    # "J033147",
    # "J044149",
    # "J224000",
    # "J224558",
]
sec_dict = ["2240-260"]

directory = str(os.environ["PROJECT"])
tarname = str(os.environ["TARGET"])

plotfuncts.plt_sed(directory, tarname)
plotfuncts.plt_secondary(directory, tarname)
plotfuncts.plt_nearby(directory, tarname)
