#!/data/bin/casa-6.1.2-7-pipeline-2020.1.0.36/bin/python3
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
)
from casaplotms import plotms
import numpy as np


def flag_ms(img_dir, visname, epoch, ATCA_band, pri, sec, tar, tar_nm):
    flagmanager(vis=visname, mode="save", versionname="before_online_flagging")
    print("Flagging antennae affected by shadowing...")
    flagdata(vis=visname, mode="shadow", tolerance=0.0, flagbackup=False)
    print("Flagging visibilities with zero amplitudes...")
    flagdata(vis=visname, mode="clip", clipzeros=True, flagbackup=False)
    print("Quacking visibilities ...")
    flagdata(
        vis=visname, mode="quack", quackinterval=5.0, quackmode="beg", flagbackup=False
    )
    flagmanager(vis=visname, mode="save", versionname="after_online_flagging")
    print(
        "Inspecting {0} amplitude as a function of channel to identify RFI...".format(
            pri
        )
    )
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
    plotms(
        vis=visname,
        field=pri,
        xaxis="channel",
        yaxis="amp",
        correlation="xy,yx",
        ydatacolumn="data",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvschan_postRFI_flag.png",
        showgui=True,
        overwrite=True,
    )
    plotms(
        vis=visname,
        field=sec,
        xaxis="channel",
        yaxis="amp",
        correlation="xy,yx",
        ydatacolumn="data",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{sec}_ampvschan_postRFI_flag.png",
        showgui=True,
        overwrite=True,
    )
    plotms(
        vis=visname,
        field="J001513",
        xaxis="channel",
        yaxis="amp",
        correlation="xy,yx",
        ydatacolumn="data",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar_nm}_ampvschan_postRFI_flag.png",
        showgui=False,
        overwrite=True,
    )
    return


def split_ms(
    src_dir,
    img_dir,
    visname,
    msname,
    epoch,
    ATCA_band,
    pri,
    sec,
    tar,
    tar_nm,
):
    os.system(f"rm -r {msname}")
    os.system("rm -r {0}.flagversions".format(msname))
    os.system("rm -r *.last")
    # have removed n_spw for mstransform and included it in the split just before imaging
    mstransform(
        vis=visname,
        outputvis=msname,
        regridms=True,
        datacolumn="data",
        mode="channel",
        field=f"{pri},{sec},{tar}",
    )
    listobs(
        vis=msname,
        listfile=f"{src_dir}/listobs_{epoch}_{ATCA_band}_{tar_nm}.dat",
        overwrite=True,
    )
    flagmanager(vis=msname, mode="save", versionname="after_transform")
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
    return


def calibrate_ms(src_dir, msname, epoch, ATCA_band, ref, pri, sec, tar, tar_nm):
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
        minblperant=3,
        solint="120s",
    )
    print(f"Performing bandpass calibration on {pri}")
    bandpass(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B0",
        field=pri,
        combine="spw,scan",
        refant=ref,
        solnorm=True,
        solint="inf",
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
        combine="spw,scan",
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
        combine="spw,scan",
        refant=ref,
        solnorm=True,
        solint="inf",
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
        combine="spw,scan",
        gaintype="G",
        calmode="ap",
        parang=True,
        solint="120s",
        minblperant=3,
        gaintable=[f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.B1"],
    )
    print(f"Deriving gain calibration using {sec}")
    gaincal(
        vis=msname,
        caltable=f"{src_dir}/cal_tables/cal_{pri}_{epoch}_{ATCA_band}.G2",
        field=sec,
        refant=ref,
        combine="spw,scan",
        gaintype="G",
        calmode="ap",
        parang=True,
        solint="120s",
        minblperant=3,
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


def inspectpostcal_ms(img_dir, msname, epoch, ATCA_band, pri, sec, tar):
    plotms(
        vis=msname,
        field=pri,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvsfreq_post_cal.png",
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
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{sec}_ampvsfreq_post_cal.png",
        showgui=False,
        overwrite=True,
    )
    plotms(
        vis=msname,
        field=pri,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="baseline",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{pri}_ampvsfreqbaseline_post_cal.png",
        showgui=False,
        overwrite=True,
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
    return


def flagcaltar_ms(src_dir, img_dir, msname, epoch, ATCA_band, pri, sec, tar, tar_nm):
    plotms(
        vis=msname,
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
    print(f"Inspecting {tar} amp vs freq Before RFI flagging")
    plotms(
        vis=msname,
        field=tar,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar_nm}_ampvsfreq_post_cal.png".format(
            epoch, ATCA_band, tar_nm
        ),
        showgui=False,
        overwrite=True,
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
    print(f"Inspecting {tar} amp vs freq AFTER RFI flagging")
    plotms(
        vis=msname,
        field=tar,
        xaxis="frequency",
        yaxis="amp",
        correlation="xx,yy",
        ydatacolumn="corrected",
        coloraxis="spw",
        plotfile=f"{img_dir}/{epoch}_{ATCA_band}_{tar_nm}_ampvsfreq_post_RFIflag.png",
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
    return


def imgmfs_ms(src_dir, msname, targetms, epoch, ATCA_band, n_spw, tar, tar_nm):
    mode = "mfs"
    nterms = 2
    niter = 3000
    if ATCA_band == "L":
        imsize = 2240
    if ATCA_band == "C":
        imsize = 1120
    if ATCA_band == "X":
        imsize = 960
    cell = "0.1arcsec"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = True
    gain = 0.01
    threshold = "2e-2Jy"

    os.system(f"rm -r {targetms}*")
    # split(vis=msname, datacolumn='corrected', field=tar, outputvis=targetms)
    mstransform(
        vis=msname,
        outputvis=targetms,
        regridms=True,
        datacolumn="corrected",
        mode="channel",
        nspw=n_spw,
        field=tar,
    )
    listobs(
        vis=targetms,
        listfile=f"{src_dir}/listobs_{epoch}_{ATCA_band}_{tar_nm}_preimage.dat",
        overwrite=True,
    )
    os.system(f"rm -r {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_mfs*")
    imagename = f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_mfs"
    flagmanager(vis=targetms, mode="save", versionname="before_selfcal")
    print("Initiating interactive cleaning on {0}".format(imagename))

    tclean(
        vis=targetms,
        imagename=imagename,
        selectdata=True,
        deconvolver="mtmfs",
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
        interactive=interactive,
        savemodel="modelcolumn",
        pbcor=False,
    )
    tclean(
        vis=targetms,
        imagename=imagename,
        selectdata=True,
        deconvolver="mtmfs",
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
        interactive=interactive,
        savemodel="modelcolumn",
        pbcor=False,
        calcres=False,
        calcpsf=False,
    )
    return


def img_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar, tar_nm):
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Preself Imaging  +\n+ + + + + + + + + + + + + + + + +"
    )
    mode = "mfs"
    nterms = 2
    niter = 3000
    if ATCA_band == "L":
        imsize = 2240
    if ATCA_band == "C":
        imsize = 1120
    if ATCA_band == "X":
        imsize = 960
    cell = "0.1arcsec"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = False
    gain = 0.01
    threshold = "5e-3Jy"
    flagmanager(vis=targetms, mode="save", versionname="preself")
    for i in range(0, n_spw):
        spw = str(i)
        print("Cleaning on band: " + str(spw))
        os.system(
            f"rm -r {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_preself*"
        )
        imagename = f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_preself",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_preself",
            selectdata=True,
            mask="",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
        )
    return


def slefcal_ms(src_dir, process_dir, targetms, epoch, ATCA_band, n_spw, tar, tar_nm):
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal Round 1  +\n+ + + + + + + + + + + + + + + + +"
    )
    mode = "mfs"
    nterms = 2
    niter = 3000
    if ATCA_band == "L":
        imsize = 2240
        solint = "240s"
        minsnr = 3.0
        minblperant = 3
    if ATCA_band == "C":
        imsize = 1120
        solint = "60s"
        minsnr = 3.0
        minblperant = 3
    if ATCA_band == "X":
        imsize = 960
        solint = "60s"
        minsnr = 3.0
        minblperant = 3
    cell = "0.1arcsec"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = False
    gain = 0.01
    threshold = "5e-4Jy"

    rmtables(f"{process_dir}/pcal1")
    gaincal(
        vis=targetms,
        caltable=f"{process_dir}/pcal1",
        combine="scan,spw",
        spwmap=[0] * n_spw,
        gaintype="GSPLINE",
        calmode="p",
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant,
    )
    applycal(
        vis=targetms,
        gaintable=f"{process_dir}/pcal1",
        spwmap=[0] * n_spw,
        parang=True,
        applymode="calonly",
        flagbackup=False,
    )
    flagmanager(vis=targetms, mode="save", versionname="post self1")

    for i in range(0, n_spw):
        spw = str(i)
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_preself.psf.tt0 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1.psf.tt0 "
        # )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_preself.psf.tt1 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1.psf.tt1 "
        # )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_preself.psf.tt2 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1.psf.tt2 "
        # )
        os.system(
            f"rm -r {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1*"
        )
        print("Cleaning on band: " + str(spw))
        imagename = f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_self1",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_self1",
            selectdata=True,
            mask="",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
        )

    threshold = "5e-5Jy"
    rmtables(f"{process_dir}/pcal2")
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal Round 2  +\n+ + + + + + + + + + + + + + + + +"
    )

    gaincal(
        vis=targetms,
        caltable=f"{process_dir}/pcal2",
        gaintable=f"{process_dir}/pcal1",
        combine="scan,spw",
        spwmap=[0] * n_spw,
        gaintype="GSPLINE",
        calmode="p",
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant,
    )
    applycal(
        vis=targetms,
        gaintable=[f"{process_dir}/pcal1", f"{process_dir}/pcal2"],
        spwmap=[0] * n_spw,
        parang=True,
        applymode="calonly",
        flagbackup=False,
    )
    flagmanager(vis=targetms, mode="save", versionname="post self2")
    for i in range(0, n_spw):
        spw = str(i)
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1.psf.tt0 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2.psf.tt0 "
        # )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1.psf.tt1 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2.psf.tt1 "
        # )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self1.psf.tt2 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2.psf.tt2 "
        # )
        os.system(
            f"rm -r {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2*"
        )
        print("Cleaning on band: " + str(spw))
        imagename = f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_self2",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_self2",
            selectdata=True,
            mask="",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
        )

    threshold = "5e-6Jy"
    print(
        "+ + + + + + + + + + + + + + + + +\n+  Self Cal Round 3  +\n+ + + + + + + + + + + + + + + + +"
    )

    rmtables(f"{process_dir}/pcal3")
    gaincal(
        vis=targetms,
        caltable=f"{process_dir}/pcal3",
        gaintable=[f"{process_dir}/pcal1", f"{process_dir}/pcal2"],
        combine="scan,spw",
        spwmap=[0] * n_spw,
        gaintype="GSPLINE",
        calmode="p",
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant,
    )
    applycal(
        vis=targetms,
        gaintable=[
            f"{process_dir}/pcal1",
            f"{process_dir}/pcal2",
            f"{process_dir}/pcal3",
        ],
        spwmap=[0] * n_spw,
        parang=True,
        applymode="calonly",
        flagbackup=False,
    )
    flagmanager(vis=targetms, mode="save", versionname="post self3")
    for i in range(0, n_spw):
        spw = str(i)
        os.system(
            f"rm -r {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self3*"
        )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2.psf.tt0 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self3.psf.tt0 "
        # )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2.psf.tt1 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self3.psf.tt1 "
        # )
        # os.system(
        #     f"mv {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self2.psf.tt2 {src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}_self3.psf.tt2 "
        # )
        print("Cleaning on band: " + str(spw))
        imagename = f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}"
        tclean(
            vis=targetms,
            imagename=imagename + "_self3",
            selectdata=True,
            mask=f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_mfs.mask",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            savemodel="modelcolumn",
            pbcor=False,
        )
        tclean(
            vis=targetms,
            imagename=imagename + "_self3",
            selectdata=True,
            mask="",
            spw=spw,
            deconvolver="mtmfs",
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
            interactive=interactive,
            pbcor=False,
            savemodel="modelcolumn",
            calcres=False,
            calcpsf=False,
        )
    return


def pbcor_ms(src_dir, targetms, epoch, ATCA_band, n_spw, tar, tar_nm):
    for i in range(0, n_spw):
        spw = str(i)
        imagename = f"{src_dir}/casa_files/{tar_nm}_{epoch}_{ATCA_band}_{spw}"
        os.system("rm -r " + imagename + "_pbcor")
        impbcor(
            imagename=f"{imagename}_self3.image.tt0",
            pbimage=f"{imagename}_self3.pb.tt0",
            outfile=f"{imagename}_self3_pbcor",
            cutoff=0.1,
            overwrite=True,
        )
    return


def measureflux_ms(
    src_dir, targetms, tar_ms, epoch, ATCA_band, sourcepar, n_spw, tar, tar_nm
):
    os.system(f"rm -r {tar_ms}")
    split(vis=targetms, datacolumn="corrected", outputvis=tar_ms)
    int_flux_c = []
    for i in range(n_spw):
        spw = str(i)
        # If things look like theyre not working, then check the source position! Chances are it can't find the source too far away from the phase centre
        outfile = f"{src_dir}/casa_files/{tar_nm}_{ATCA_band}_{epoch}_{spw}.cl"
        os.system(f"rm -r {outfile}")
        uvmodelfit(
            vis=tar_ms,
            niter=15,
            comptype="P",
            spw=spw,
            sourcepar=sourcepar,
            outfile=outfile,
            field="0",
        )
        tbl = table(outfile)
        flux = tbl.getcell("Flux", 0)[0].astype("float64")
        int_flux_c.append(flux)
        print(flux)
    if ATCA_band == "C":
        np.savetxt(
            f"{src_dir}/{tar_nm}_{epoch}_{ATCA_band}.csv",
            int_flux_c,
            delimiter=",",
            header="S_Cband",
        )
        print(int_flux_c)
    elif ATCA_band == "X":
        np.savetxt(
            f"{src_dir}/{tar_nm}_{epoch}_{ATCA_band}.csv",
            int_flux_c,
            delimiter=",",
            header="S_Xband",
        )
        print(int_flux_c)
    elif ATCA_band == "L":
        int_flux_l = np.array(int_flux_c[::-1])
        np.savetxt(
            f"{src_dir}/{tar_nm}_{epoch}_{ATCA_band}.csv",
            int_flux_l,
            header="S_Lband",
            delimiter=",",
        )
        print(int_flux_l)
    return
