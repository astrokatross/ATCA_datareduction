#!/usr/bin/python3
# This is a script with functions needed for the data reduction, use run_process.py to actually analyse
# Updated from B. Quici script By K.Ross 19/5/21


# TODO: Make a try/else situation for if it tries to make a psf but can't, then try calibrating with a lower snr. Not sure about exact implementation though
import os
from casacore.tables import table
from casatasks import (
    flagmanager,
    flagdata,
    mstransform,
    listobs,
    setjy,
    gaincal,
    bandpass,
    fluxscale,
    applycal,
    tclean,
    rmtables,
    impbcor,
    split,
    uvmodelfit,
    exportfits,
)
import numpy as np
from casaplotms import plotms
import matplotlib.pyplot as plt
from casatools import image as IA
from astropy.wcs import WCS
from astropy.visualization import simple_norm

ia = IA()
plt.rcParams["font.family"] = "serif"


def buildImage(imname="", chan=0):
    ia.open(imname)
    pix = ia.getchunk()[:, :, 0, chan]
    csys = ia.coordsys()
    ia.close()

    rad_to_deg = 180 / np.pi
    w = WCS(naxis=2)
    w.wcs.crpix = csys.referencepixel()["numeric"][0:2]
    w.wcs.cdelt = csys.increment()["numeric"][0:2] * rad_to_deg
    w.wcs.crval = csys.referencevalue()["numeric"][0:2] * rad_to_deg
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]

    return pix, w


def flag_ms(visname):  # , rawname1, rawname2, rawname3):
    print("Flagging antennas affected by shadowing...")
    # importatca(
    #     vis=visname, files=[rawname1, rawname2, rawname3], options="birdie,noac", edge=4
    # )
    flagmanager(vis=visname, mode="save", versionname="before_online_flagging")
    print("Flagging antennas affected by shadowing...")
    flagdata(vis=visname, mode="shadow", tolerance=0.0, flagbackup=False)
    print("Flagging visibilities with zero amplitudes...")
    flagdata(vis=visname, mode="clip", clipzeros=True, flagbackup=False)
    print("Quacking visibilities ...")
    flagdata(
        vis=visname, mode="quack", quackinterval=5.0, quackmode="beg", flagbackup=False
    )
    flagmanager(vis=visname, mode="save", versionname="after_online_flagging")
    flagdata(
        vis=visname,
        mode="tfcrop",
        datacolumn="data",
        action="apply",
        display="report",
        flagbackup=True,
        extendpols=True,
        correlation="",
        flagdimension="freqtime",
        growtime=95.0,
        growfreq=95.0,
        timecutoff=4.0,
        timefit="line",
        freqfit="poly",
        maxnpieces=5,
        combinescans=False,
        ntime="scan",
        extendflags=False,
    )
    print("Extending flags to all correlations")
    flagdata(
        vis=visname,
        mode="extend",
        action="apply",
        display="report",
        flagbackup=False,
        extendpols=True,
        correlation="",
        growtime=95.0,
        growfreq=95.0,
        growaround=True,
        flagneartime=False,
        flagnearfreq=False,
        combinescans=False,
        ntime="scan",
    )
    return


def split_ms(src_dir, img_dir, visname, msname, epoch, ATCA_band, pri, sec, tar, n_spw):
    os.system(f"rm -r {msname}")
    os.system(f"rm -r {msname}.flagversions")
    os.system("rm -r *.last")
    # have removed n_spw for mstransform and included it in the split just before imaging
    mstransform(
        vis=visname,
        outputvis=msname,
        datacolumn="data",
        field=f"{pri},{sec},{tar}",
        # nspw=n_spw,
        regridms=True,
        # field=f"{sec},{tar}",
        # scan="3,>90"
    )
    listobs(
        vis=msname,
        listfile=f"{src_dir}/listobs_{epoch}_{ATCA_band}_{tar}.dat",
        overwrite=True,
    )
    flagmanager(vis=msname, mode="save", versionname="after_transform")
    return


def calibrate_ms(src_dir, msname, epoch, ATCA_band, ref, pri, sec, tar):
    os.system(f"rm -r {src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.F0")
    setjy(
        vis=msname,
        field=pri,
        scalebychan=True,
        standard="Perley-Butler 2010",
        usescratch=True,
    )
    print(f"Performing gain calibration on {pri}")
    gaincal(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G0",
        field=pri,
        refant=ref,
        gaintype="G",
        calmode="p",
        parang=True,
        # minblperant=3,
        solint="60s",
    )
    print(f"Performing bandpass calibration on {pri}")
    bandpass(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B0",
        field=pri,
        refant=ref,
        solnorm=True,
        solint="120s",
        bandtype="B",
        gaintable=[f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G0"],
        parang=True,
    )
    print(f"Determining gains on {sec}")
    gaincal(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G1",
        field=pri + "," + sec,
        refant=ref,
        gaintype="G",
        calmode="ap",
        parang=True,
        solint="120s",
        gaintable=[f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B0"],
    )
    bandpass(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1",
        field=pri,
        refant=ref,
        solnorm=True,
        solint="120s",
        bandtype="B",
        gaintable=[f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G1"],
        parang=True,
    )
    print(f"Deriving gain calibration using {pri}")
    gaincal(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G2",
        field=pri,
        refant=ref,
        gaintype="G",
        calmode="ap",
        parang=True,
        solint="60s",
        # minblperant=3,
        gaintable=[f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1"],
    )
    print(f"Deriving gain calibration using {sec}")
    gaincal(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G2",
        field=sec,
        refant=ref,
        gaintype="G",
        calmode="ap",
        parang=True,
        solint="60s",
        # minblperant=3,
        gaintable=[f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1"],
        append=True,
    )
    print(
        "Correcting the flux scale using comparison between the primary and secondary calibrator."
    )
    fluxscale(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G2",
        fluxtable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.F0",
        reference=pri,
    )
    flagmanager(vis=msname, mode="save", versionname="before_applycal")
    return


def applycal_ms(src_dir, msname, epoch, ATCA_band, pri, sec, tar):
    applycal(
        vis=msname,
        gaintable=[
            f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1",
            f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.F0",
        ],
        gainfield=[pri, pri, pri],
        field=f"{pri}",
        parang=True,
        flagbackup=False,
    )
    applycal(
        vis=msname,
        gaintable=[
            f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1",
            f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.F0",
        ],
        gainfield=[pri, pri, sec],
        field=f"{sec},{tar}",
        parang=True,
        flagbackup=False,
    )
    return


def flagcal_ms(img_dir, msname, epoch, ATCA_band, pri, sec):
    flagmanager(vis=msname, mode="save", versionname="before_rflag")
    flagdata(
        vis=msname,
        mode="rflag",
        field=pri,
        datacolumn="corrected",
        action="apply",
        display="report",
        correlation="ABS_ALL",
        timedevscale=3.0,
        freqdevscale=3.0,
        winsize=3,
        combinescans=True,
        ntime="9999999min",
        extendflags=False,
        flagbackup=False,
    )
    flagdata(
        vis=msname,
        mode="rflag",
        field=sec,
        datacolumn="corrected",
        action="apply",
        display="report",
        correlation="ABS_ALL",
        timedevscale=3.0,
        freqdevscale=3.0,
        winsize=3,
        combinescans=True,
        ntime="9999999min",
        extendflags=False,
        flagbackup=False,
    )
    flagdata(
        vis=msname,
        mode="extend",
        field=pri + "," + sec,
        action="apply",
        display="report",
        flagbackup=False,
        extendpols=True,
        correlation="",
        growtime=95.0,
        growfreq=95.0,
        growaround=True,
        flagneartime=False,
        flagnearfreq=False,
        combinescans=True,
        ntime="9999999min",
    )

    return


def flagcaltar_ms(src_dir, img_dir, msname, epoch, ATCA_band, pri, sec, tar):
    applycal(
        vis=msname,
        gaintable=[
            f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1",
            f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.F0",
        ],
        gainfield=[pri, pri, sec],
        field=tar,
        parang=True,
        flagbackup=False,
    )
    flagdata(
        vis=msname,
        mode="rflag",
        field=tar,
        datacolumn="corrected",
        action="apply",
        display="report",
        correlation="ABS_ALL",
        timedevscale=3.0,
        freqdevscale=3.0,
        winsize=3,
        combinescans=True,
        ntime="9999999min",
        extendflags=False,
        flagbackup=False,
    )
    return


def imgmfs_ms(src_dir, msname, targetms, epoch, ATCA_band, n_spw, tar):
    mode = "mfs"
    nterms = 1
    niter = 3000
    antenna = "0~6,0~6"
    if (
        (tar in ["J001513", "J224408", "J223933"])
        and (epoch == "epoch2")
        and (ATCA_band == "X")
    ):
        n_spw = n_spw - 1
    else:
        n_spw = n_spw
    if ATCA_band == "L":
        imsize = 2250
        cell = "1arcsec"
        threshold = "2e-2Jy"
    if ATCA_band == "C":
        imsize = 2250
        cell = "0.5arcsec"
        threshold = "2e-2Jy"
    if ATCA_band == "X":
        imsize = 1152
        cell = "0.2arcsec"
        threshold = "2e-2Jy"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = True
    gain = 0.01
    if epoch == "2021-10-15":
        uvrange = "<1000"
        cell = "1arcsec"
        imsize = 1152
    else:
        uvrange = ""
    print(uvrange)
    os.system(f"rm -r {targetms}*")
    # split(vis=msname, datacolumn='corrected', field=tar, outputvis=targetms)
    mstransform(
        vis=msname,
        outputvis=targetms,
        datacolumn="corrected",
        nspw=n_spw,
        regridms=True,
        field=tar,
    )
    listobs(
        vis=targetms,
        listfile=f"{src_dir}/listobs_{epoch}_{ATCA_band}_{tar}_preimage.dat",
        overwrite=True,
    )
    os.system(f"rm -r {src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs*")
    imagename = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs"
    flagmanager(vis=targetms, mode="save", versionname="before_selfcal")
    print("Initiating interactive cleaning on {0}".format(imagename))

    tclean(
        vis=targetms,
        imagename=imagename,
        gain=gain,
        specmode=mode,
        nterms=nterms,
        niter=niter,
        threshold=threshold,
        imsize=imsize,
        cell=cell,
        stokes=stokes,
        weighting=weighting,
        robust=robust,
        antenna=antenna,
        interactive=interactive,
        savemodel="modelcolumn",
        pbcor=False,
        uvrange=uvrange,
    )
    tclean(
        vis=targetms,
        imagename=imagename,
        gain=gain,
        specmode=mode,
        nterms=nterms,
        niter=0,
        threshold=threshold,
        imsize=imsize,
        cell=cell,
        stokes=stokes,
        weighting=weighting,
        robust=robust,
        antenna=antenna,
        interactive=interactive,
        savemodel="modelcolumn",
        pbcor=False,
        calcres=False,
        calcpsf=False,
        uvrange=uvrange,
    )
    return


def img_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar):
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Preself Imaging  +\n+ + + + + + + + + + + + + + + + +"
    )
    antenna = "0~6,0~6"
    mode = "mfs"
    nterms = 1
    niter = 3000
    if (
        (tar in ["J001513", "J224408", "J223933"])
        and (epoch == "epoch2")
        and (ATCA_band == "X")
    ):
        n_spw = n_spw - 1
    else:
        n_spw = n_spw
    if ATCA_band == "L":
        imsize = 2250
        cell = "1arcsec"
    if ATCA_band == "C":
        imsize = 2250
        cell = "0.5arcsec"
    if ATCA_band == "X":
        imsize = 1152
        cell = "0.2arcsec"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = False
    gain = 0.01
    threshold = "5e-3Jy"
    if epoch == "2021-10-15":
        uvrange = "<1000"
        cell = "1arcsec"
        imsize = 1152
    else:
        uvrange = ""
    flagmanager(vis=targetms, mode="save", versionname="preself")
    for i in range(0, n_spw):
        spw = str(i)
        print("Cleaning on band: " + str(spw))
        os.system(
            f"rm -r {src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_preself*"
        )
        imagename = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_preself",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=niter,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
            uvrange=uvrange,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_preself",
            selectdata=True,
            mask="",
            spw=spw,
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=0,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
            uvrange=uvrange,
        )
        model_im = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_preself.model"
        plt.subplots(1, 1, figsize=(18, 12))
        pix, w = buildImage(model_im)
        ax = plt.subplot(1, 1, 1, projection=w)
        p1 = int(pix.shape[0] * 0.25)
        p2 = int(pix.shape[0] * 0.75)

        im = ax.imshow(
            pix[p1:p2, p1:p2].transpose(), origin="lower", cmap=plt.cm.plasma
        )
        plt.colorbar(im, ax=ax)
        ax.set_xlabel("Right Ascension", fontsize=30)
        ax.set_ylabel("Declination", fontsize=30)
        plt.title(f"{tar} {epoch} {ATCA_band} mask", fontsize=30)
        plt.savefig(
            f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_preself_model.png"
        )
        plt.close()
    return


def slefcal_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar):
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal Round 1  +\n+ + + + + + + + + + + + + + + + +"
    )
    mode = "mfs"
    nterms = 2
    niter = 3000
    antenna = "0~6,0~6"
    if (
        (tar in ["J001513", "J224408", "J223933"])
        and (epoch == "epoch2")
        and (ATCA_band == "X")
    ):
        n_spw = n_spw - 1
    else:
        n_spw = n_spw
    if ATCA_band == "L":
        imsize = 2250
        solint = "60s"
        minsnr = 3.0
        cell = "1arcsec"
        minblperant = 3
    if ATCA_band == "C":
        imsize = 2250
        solint = "60s"
        minsnr = 3.0
        minblperant = 3
        cell = "0.5arcsec"
    if ATCA_band == "X":
        imsize = 1152
        solint = "60s"
        minsnr = 3.0
        minblperant = 3
        cell = "0.2arcsec"
    if epoch == "2021-10-15":
        print(epoch)
        uvrange = "<1000"
        cell = "1arcsec"
        imsize = 1152
    else:
        uvrange = ""
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = False
    gain = 0.01
    threshold = "5e-4Jy"

    rmtables(f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}")
    gaincal(
        vis=targetms,
        caltable=f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}",
        combine="scan,spw",
        spwmap=[0] * n_spw,
        gaintype="G",
        calmode="p",
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant,
    )
    applycal(
        vis=targetms,
        gaintable=f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}",
        spwmap=[0] * n_spw,
        parang=True,
        applymode="calonly",
        flagbackup=False,
    )
    flagmanager(vis=targetms, mode="save", versionname="post self1")

    for i in range(0, n_spw):
        spw = str(i)
        os.system(f"rm -r {src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_self1*")
        print("Cleaning on band: " + str(spw))
        imagename = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_self1",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            startmodel=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_preself.model",
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=niter,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
            uvrange=uvrange,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_self1",
            selectdata=True,
            mask="",
            spw=spw,
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=0,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
            uvrange=uvrange,
        )

    threshold = "5e-5Jy"
    rmtables(f"{src_dir}/cal_tables/pcal2_{tar}_{epoch}_{ATCA_band}")
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal Round 2  +\n+ + + + + + + + + + + + + + + + +"
    )

    gaincal(
        vis=targetms,
        caltable=f"{src_dir}/cal_tables/pcal2_{tar}_{epoch}_{ATCA_band}",
        gaintable=f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}",
        combine="scan,spw",
        spwmap=[0] * n_spw,
        gaintype="G",
        calmode="p",
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant,
    )
    applycal(
        vis=targetms,
        gaintable=[
            f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}",
            f"{src_dir}/cal_tables/pcal2_{tar}_{epoch}_{ATCA_band}",
        ],
        spwmap=[[0] * n_spw, [0] * n_spw],
        parang=True,
        applymode="calonly",
        flagbackup=False,
    )
    flagmanager(vis=targetms, mode="save", versionname="post self2")
    for i in range(0, n_spw):
        spw = str(i)
        os.system(f"rm -r {src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_self2*")
        print("Cleaning on band: " + str(spw))
        imagename = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_self2",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            startmodel=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_self1.model",
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=niter,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
            uvrange=uvrange,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_self2",
            selectdata=True,
            mask="",
            spw=spw,
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=0,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
            uvrange=uvrange,
        )

    threshold = "5e-6Jy"
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal Round 3  +\n+ + + + + + + + + + + + + + + + +"
    )

    rmtables(f"{src_dir}/cal_tables/pcal3_{tar}_{epoch}_{ATCA_band}")
    gaincal(
        vis=targetms,
        caltable=f"{src_dir}/cal_tables/pcal3_{tar}_{epoch}_{ATCA_band}",
        gaintable=[
            f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}",
            f"{src_dir}/cal_tables/pcal2_{tar}_{epoch}_{ATCA_band}",
        ],
        combine="scan,spw",
        spwmap=[[0] * n_spw, [0] * n_spw],
        gaintype="G",
        calmode="p",
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant,
    )
    applycal(
        vis=targetms,
        gaintable=[
            f"{src_dir}/cal_tables/pcal1_{tar}_{epoch}_{ATCA_band}",
            f"{src_dir}/cal_tables/pcal2_{tar}_{epoch}_{ATCA_band}",
            f"{src_dir}/cal_tables/pcal3_{tar}_{epoch}_{ATCA_band}",
        ],
        spwmap=[[0] * n_spw, [0] * n_spw, [0] * n_spw],
        parang=True,
        applymode="calonly",
        flagbackup=False,
    )
    flagmanager(vis=targetms, mode="save", versionname="post self3")
    for i in range(0, n_spw):
        spw = str(i)
        os.system(f"rm -r {src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_self3*")
        print("Cleaning on band: " + str(spw))
        imagename = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_self3",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs.mask",
            startmodel=f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}_self1.model",
            spw=spw,
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=niter,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
            uvrange=uvrange,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_self3",
            selectdata=True,
            mask="",
            spw=spw,
            gain=gain,
            specmode=mode,
            nterms=nterms,
            niter=0,
            threshold=threshold,
            imsize=imsize,
            cell=cell,
            stokes=stokes,
            weighting=weighting,
            robust=robust,
            antenna=antenna,
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
            uvrange=uvrange,
        )
    return


def pbcor_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar):
    if (
        (tar in ["J001513", "J224408", "J223933"])
        and (epoch == "epoch2")
        and (ATCA_band == "X")
    ):
        n_spw = n_spw - 1
    else:
        n_spw = n_spw
    for i in range(0, n_spw):
        spw = str(i)
        imagename = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_{spw}"
        os.system("rm -r " + imagename + "_self3_pbcor")
        impbcor(
            imagename=f"{imagename}_self3.image",
            pbimage=f"{imagename}_self3.pb",
            outfile=f"{imagename}_self3_pbcor",
            cutoff=0.1,
            overwrite=True,
        )
    return


def measureflux_ms(
    src_dir, targetms, tar_ms, epoch, ATCA_band, sourcepar, n_spw, tar, calibrator, timerange="",
):
    os.system(f"rm -r {tar_ms}")
    if calibrator == "TARGET":
        split(
            vis=targetms, datacolumn="corrected", field=tar, outputvis=tar_ms
        )  # , nspw=n_spw)
    else:
        mstransform(
            vis=targetms,
            outputvis=tar_ms,
            datacolumn="corrected",
            nspw=n_spw,
            regridms=True,
            field=tar,
        )

    int_flux_c = []
    if tar == "J215436":
        uvrange = "<1000"
    else:
        uvrange = ""
    # if (
    #     (tar in ["J001513", "J224408", "J223933", "1934_cal_l"])
    #     and (epoch == "epoch1")
    #     and (ATCA_band == "L")
    # ):
    #     n_spw = n_spw
    # else:
    #     n_spw = n_spw
    for i in range(n_spw):
        spw = str(i)
        # If things look like theyre not working, then check the source position! Chances are it can't find the source too far away from the phase centre
        outfile = f"{src_dir}/casa_files/{tar}_{ATCA_band}_{epoch}_{spw}{timerange}.cl"
        os.system(f"rm -r {outfile}")
        uvmodelfit(
            vis=tar_ms,
            niter=25,
            comptype="P",
            spw=spw,
            sourcepar=sourcepar,
            outfile=outfile,
            uvrange=uvrange,
            field="0",
            selectdata=True,
            timerange=timerange,
        )
        tbl = table(outfile)
        flux = tbl.getcell("Flux", 0)[0].astype("float64")
        int_flux_c.append(flux)
        print(flux)
    if ATCA_band == "C":
        np.savetxt(
            f"{src_dir}/{tar}_{epoch}_{ATCA_band}{timerange}.csv",
            int_flux_c,
            delimiter=",",
            header="S_Cband",
        )
        print(int_flux_c)
    elif ATCA_band == "X":
        np.savetxt(
            f"{src_dir}/{tar}_{epoch}_{ATCA_band}{timerange}.csv",
            int_flux_c,
            delimiter=",",
            header="S_Xband",
        )
        print(int_flux_c)
    elif ATCA_band == "L":
        # int_flux_l = np.array(int_flux_c[::-1])
        int_flux_l = int_flux_c
        np.savetxt(
            f"{src_dir}/{tar}_{epoch}_{ATCA_band}{timerange}.csv",
            int_flux_l,
            header="S_Lband",
            delimiter=",",
        )
        print(int_flux_l)
    return


def inspection_plots(
    src_dir, img_dir, visname, msname, targetms, epoch, ATCA_band, pri, sec, tar
):
    # plotms(
    #     vis=visname,
    #     field=pri,
    #     xaxis="channel",
    #     yaxis="amp",
    #     correlation="xy,yx",
    #     ydatacolumn="data",
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvschan_postRFI_flag.png",
    #     showgui=False,
    #     overwrite=True,
    # )
    # plotms(
    #     vis=visname,
    #     field=sec,
    #     xaxis="channel",
    #     yaxis="amp",
    #     correlation="xy,yx",
    #     ydatacolumn="data",
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{sec}_ampvschan_postRFI_flag.png",
    #     showgui=False,
    #     overwrite=True,
    # )
    # plotms(
    #     vis=targetms,
    #     field=tar,
    #     xaxis="channel",
    #     yaxis="amp",
    #     correlation="xy,yx",
    #     ydatacolumn="data",
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_ampvschan_postRFI_flag.png",
    #     showgui=False,
    #     overwrite=True,
    # )
    plotms(
        vis=msname,
        field=pri,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="data",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvsfreq_pre_cal.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=msname,
        field=sec,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="data",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{sec}_ampvsfreq_pre_cal.png",
        showgui=False,
        overwrite=True,
    )
    # plotms(
    #     vis=msname,
    #     field=pri,
    #     xaxis="frequency",
    #     yaxis="amp",
    #     correlation="xx,yy",
    #     ydatacolumn="corrected",
    #     coloraxis="spw",
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvsfreq_post_cal.png",
    #     showgui=False,
    #     overwrite=True,
    # )
    # plotms(
    #     vis=msname,
    #     field=sec,
    #     xaxis="frequency",
    #     yaxis="amp",
    #     correlation="xx,yy",
    #     ydatacolumn="corrected",
    #     coloraxis="spw",
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{sec}_ampvsfreq_post_cal.png",
    #     showgui=False,
    #     overwrite=True,
    # )
    # plotms(
    #     vis=msname,
    #     field=pri,
    #     xaxis="frequency",
    #     yaxis="amp",
    #     correlation="xx,yy",
    #     ydatacolumn="corrected",
    #     coloraxis="baseline",
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvsfreqbaseline_post_cal.png",
    #     showgui=False,
    #     overwrite=True,
    # )
    plotms(
        vis=msname,
        field=pri,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvsfreq_post_RFIflag.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=msname,
        field=sec,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{sec}_ampvsfreq_post_RFIflag.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=targetms,
        field=tar,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_ampvsfreq_post_cal.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=targetms,
        field=tar,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="data",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_ampvsfreq_pre_cal.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=targetms,
        field=tar,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_ampvsfreq_post_RFIflag.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=msname,
        field=tar,
        xaxis="u",
        yaxis="v",
        plotfile=f"{img_dir}/{epoch}_{tar}_uv.png",
        overwrite=True,
        showgui=False,
    )
    # plotms(
    #     vis=f"{src_dir}/pcal1",
    #     field=tar,
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_pcal1.png",
    #     overwrite=True,
    #     showgui=False,
    # )
    # plotms(
    #     vis=f"{src_dir}/pcal2",
    #     field=tar,
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_pcal2.png",
    #     overwrite=True,
    #     showgui=False,
    # )
    # plotms(
    #     vis=f"{src_dir}/pcal3",
    #     field=tar,
    #     plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar}_pcal3.png",
    #     overwrite=True,
    #     showgui=False,
    # )
    return


def export_fitspng(src_dir, n_spw, epoch, ATCA_band, tar):
    for i in range(0, n_spw):
        spw = str(i)
        imagename = f"{tar}_{epoch}_{ATCA_band}_{spw}"
        exportfits(
            imagename=f"{src_dir}/casa_files/{imagename}_preself.image",
            fitsimage=f"{src_dir}/images/{imagename}_preself.fits",
            overwrite=True,
        )
        exportfits(
            imagename=f"{src_dir}/casa_files/{imagename}_self1.image",
            fitsimage=f"{src_dir}/images/{imagename}_self1.fits",
            overwrite=True,
        )
        exportfits(
            imagename=f"{src_dir}/casa_files/{imagename}_self2.image",
            fitsimage=f"{src_dir}/images/{imagename}_self2.fits",
            overwrite=True,
        )
        exportfits(
            imagename=f"{src_dir}/casa_files/{imagename}_self3.image",
            fitsimage=f"{src_dir}/images/{imagename}_self3.fits",
            overwrite=True,
        )
        exportfits(
            imagename=f"{src_dir}/casa_files/{imagename}_self3_pbcor",
            fitsimage=f"{src_dir}/images/{imagename}_self3_pbcor.fits",
            overwrite=True,
        )
        extensions = ["preself", "self1", "self2", "self3"]
        for ext in extensions:
            imname = f"{src_dir}/casa_files/{imagename}_{ext}.image"
            plt.subplots(1, 1, figsize=(18, 12))
            pix, w = buildImage(imname)
            ax = plt.subplot(1, 1, 1, projection=w)
            p1 = int(pix.shape[0] * 0.25)
            p2 = int(pix.shape[0] * 0.75)

            norm = simple_norm(pix[p1:p2, p1:p2].transpose(), "sqrt")
            im = ax.imshow(
                pix[p1:p2, p1:p2].transpose(),
                origin="lower",
                cmap=plt.cm.plasma,
                norm=norm,
            )
            plt.colorbar(im, ax=ax)
            ax.set_xlabel("Right Ascension", fontsize=30)
            ax.set_ylabel("Declination", fontsize=30)
            plt.title(f"{imagename} {ext}", fontsize=30)
            plt.savefig(f"{src_dir}/images/{imagename}_{ext}.png")
            plt.close()
    mask_im = f"{src_dir}/casa_files/{tar}_{epoch}_{ATCA_band}_mfs.mask"
    plt.subplots(1, 1, figsize=(18, 12))
    pix, w = buildImage(mask_im)
    ax = plt.subplot(1, 1, 1, projection=w)
    p1 = int(pix.shape[0] * 0.25)
    p2 = int(pix.shape[0] * 0.75)

    im = ax.imshow(pix[p1:p2, p1:p2].transpose(), origin="lower", cmap=plt.cm.plasma)
    plt.colorbar(im, ax=ax)
    ax.set_xlabel("Right Ascension", fontsize=30)
    ax.set_ylabel("Declination", fontsize=30)
    plt.title(f"{tar} {epoch} {ATCA_band} mask", fontsize=30)
    plt.savefig(f"{src_dir}/images/{tar}_{epoch}_{ATCA_band}_mask.png")
    plt.close()
    return
