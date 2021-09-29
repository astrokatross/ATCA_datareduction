#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds 
# By K.Ross 29/09/21

import scipy.optimize as opt
import numpy as np


# Originally in Tingay, de Kool 2003
def singSSA(freq, S_norm, beta, peak_freq):  # Single SSA model
    return (
        S_norm
        * ((freq / peak_freq) ** (-(beta - 1) / 2))
        * (1 - np.exp(-((freq / peak_freq) ** (-(beta + 4) / 2))))
        / ((freq / peak_freq) ** (-(beta + 4) / 2))
    )


    poptpowlaw, pcovpowlaw = opt.curve_fits(singSSA(), freqplot, fluxplot, p0 = p0pow, sigma = flux_errplot, maxfev = 10000)