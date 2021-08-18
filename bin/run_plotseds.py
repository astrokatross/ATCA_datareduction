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
    "J001513",
    "J015445",
    "J020507",
    "J021246",
    "J022744",
    "J024838",
    "J032213",
    "J032836",
    "J033023",
    "J042502",
    "J044033",
    "J044737",
    "J052824",
    "J223933",
    "J224408",
    "J215436",
]

directory = str(os.environ["PROJECT"])
tarname = str(os.environ["TARGET"])

# tarname = "J224408"
# directory = "/data/ATCA/ATCA_datareduction/"
# for tarname in source_dict:

plotfuncts.plt_sed(directory, tarname)
