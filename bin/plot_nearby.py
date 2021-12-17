# Script just to plot the nearby MWA sources because it's too messy to do somewhere else

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd


subchans_dict = {
    "69": ["072-080", "080-088", "088-095", "095-103"],
    "93": ["103-111", "111-118", "118-126", "126-134"],
    "121": ["139-147", "147-154", "154-162", "162-170"],
    "145": ["170-177", "177-185", "185-193", "193-200"],
    "169": ["200-208", "208-216", "216-223", "223-231"],
}
epochs = ["2020-04", "2020-05", "2020-07", "2020-09"]
chans = ["69", "93", "121", "145", "169"]
exts = [
    "107",
    "115",
    "123",
    "130",
    "143",
    "150",
    "158",
    "166",
    "174",
    "181",
    "189",
    "197",
    "204",
    "212",
    "220",
    "227",
]
mwa_2013_fluxes = ["S_076", "S_084", "S_092", "S_099"]
mwa_2014_fluxes = []
mwa_2013_errors = ["S_076_err", "S_084_err", "S_092_err", "S_099_err"]
mwa_2014_errors = []
for ext in exts:
    mwa_2013_fluxes.append(f"S_{ext}_yr1")
    mwa_2014_fluxes.append(f"S_{ext}_yr2")
    mwa_2013_errors.append(f"local_rms_{ext}_yr1")
    mwa_2014_errors.append(f"local_rms_{ext}_yr2")


def read_MWA_fluxes(directory, cat_nm, name):
    master_pop_pd = pd.read_csv(f"{directory}/master_pop_extended.csv")
    print(name)
    mask = master_pop_pd["Name"] == name
    src_pd = master_pop_pd[mask]

    mwa_flux = np.full([len(epochs), 20], np.nan)
    err_mwa = np.full([len(epochs), 20], np.nan)

    for j in range(len(epochs)):
        epoch = epochs[j]
        chan_flux = []
        err_chan_flux = []
        if epoch == "2013":
            chan_flux = np.squeeze(src_pd[mwa_2013_fluxes].values)
            err_chan_flux = np.squeeze(
                np.sqrt(src_pd[mwa_2013_errors]) ** 2 + (0.02 * chan_flux) ** 2
            )
        elif epoch == "2014":
            buffer = np.full(4, np.nan)

            chan_flux_temp = np.squeeze(src_pd[mwa_2014_fluxes].values)
            err_chan_flux_temp = np.squeeze(
                np.sqrt(src_pd[mwa_2014_errors]) ** 2 + (0.02 * chan_flux_temp) ** 2
            )
            chan_flux = np.hstack((buffer, chan_flux_temp))
            err_chan_flux = np.hstack((buffer, err_chan_flux_temp))
        else:
            for i in range(len(chans)):
                subchans = subchans_dict[chans[i]]
                chan = chans[i]
                if epoch == "2020-04":
                    percentage = 0.05
                else:
                    percentage = 0.02
                for subchan in subchans:
                    try:
                        src_mwa_pd = pd.read_csv(
                            f"{directory}/{epoch}/{chan}/minimosaic/{cat_nm}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv"
                        )
                        mask = src_mwa_pd["Name"] == name
                        src_pd = src_mwa_pd[mask]
                        if src_pd.empty is False:
                            mwa_flux_chan = np.squeeze(src_pd["int_flux"].values)
                            mwa_errs_chan = np.squeeze(
                                np.sqrt(src_pd["local_rms"]) ** 2
                                + (percentage * mwa_flux_chan) ** 2
                            )
                            chan_flux.append(mwa_flux_chan)
                            err_chan_flux.append(mwa_errs_chan)
                        else:
                            chan_flux.append(np.nan)
                            err_chan_flux.append(np.nan)
                            pass
                    except (FileNotFoundError, KeyError):
                        chan_flux.append(np.nan)
                        err_chan_flux.append(np.nan)
                        pass

        # mwa_flux[j] = chan_flux
        # err_mwa[j] = err_chan_flux
        mwa_flux = 0
    return mwa_flux
