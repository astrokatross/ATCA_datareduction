#!/usr/bin/python
# This is a script to analyse the ATCA data, it is the python script with the details of reduction
# TODO: create an executable scriptfor terminal that pipes all the needed variables into this python script and executes in casa
# Updated from B. Quici script By K.Ross 19/5/21

# Importing relevant python packages
from recipes.atcapolhelpers import qufromgain
from pathlib2 import Path
from astropy.io import fits, votable
import os
import pyfits

# Reading in variables from the bash script (TODO: setup actual bash script lol)
data_dir = os.environ['PROJECT']
src = os.environ['SRC']
epoch = os.environ['EPOCH']
ATCA_band = os.environ["BAND"]
pri = os.environ['PRIMARY_CALIBRATOR']
sec = os.environ['SECONDARY_CALIBRATOR']
tar = os.environ['TARGET']
tar_nm = os.environ['TARGET_NAME']
sourcepar = [
    0.1, 12.6, -0.9
]  # I think the best way for this is to have a dictionary of source and their position which I call based on tar_nm. This will be time consuming but probs worth and I have all the values so I just gotta transfer

# Defining constant variables for all sources
src_dir = data_dir + src
process_dir = data_dir + "processing/"
img_dir = data_dir + src + "images/"
visname = data_dir + "{0}_{1}.ms".format(epoch, ATCA_band)
msname = "{0}_{1}_{2}.ms".format(epoch, ATCA_band, tar_nm)
targetms = "{0}_{1}_{2}_img.ms".format(epoch, ATCA_band, tar_nm)
tar_ms = "{0}_{1}_{2}_selfcal.ms".format(epoch, ATCA_band, tar_nm)
ref = "CA04"
if ATCA_band == "L":
    n_spw = 8
    if_centre = 0
elif ATCA_band == "C":
    n_spw = 5
    if_centre = 0
elif ATCA_band == "X":
    n_spw = 4
    if_centre = 1

# TODO: fix this bad boy
pri_cal_no = 0
sec_cal_no = 3
tar_cal_no = 5

pri_cal_id = 0
tar_cal_id = 2
sec_cal_id = 1

pipeline_options = raw_input("Flag {0} for RFI ? ".format(visname))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    print("+ + + + + + + + + + + + +")
    print("+ Preliminary Flagging  +")
    print("+ + + + + + + + + + + + +")
    flagmanager(vis=visname, mode="save", versionname="before_online_flagging")
    print("Flagging antennae affected by shadowing...")
    flagdata(vis=visname, mode="shadow", tolerance=0.0, flagbackup=False)
    print("Flagging visibilities with zero amplitudes...")
    flagdata(vis=visname, mode="clip", clipzeros=True, flagbackup=False)
    print("Quacking visibilities ...")
    flagdata(vis=visname,
             mode="quack",
             quackinterval=5.0,
             quackmode="beg",
             flagbackup=False)
    flagmanager(vis=visname, mode="save", versionname="after_online_flagging")
    print(
        "Inspecting {0} amplitude as a function of channel to identify RFI...".
        format(pri))
    plotms(vis=visname,
           field=pri,
           xaxis="channel",
           yaxis="amp",
           correlation="xy,yx",
           ydatacolumn="data",
           plotfile="{0}_{1}_{2}_ampvschan_pre_RFIflag_.png".format(
               epoch, ATCA_band, pri),
           showgui=False,
           overwrite=True
           )  # change all the plotms outputs to encode the IF band as well
    print("+ + + + + + + + +")
    print("+ RFI Flagging  +")
    print("+ + + + + + + + +")
    print(
        "Using 'mode=tfcrop' to calculate flags based on time & frequency window."
    )
    flagdata(vis=visname,
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
             timecutoff=4.,
             timefit="line",
             freqfit="poly",
             maxnpieces=5,
             combinescans=False,
             ntime="scan",
             extendflags=False)
    print("Extending flags to all correlations")
    flagdata(vis=visname,
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
             ntime="scan")
    print(
        "Inspecting {0} amplitude as a function of channel to inspect RFI before_online_flagging..."
        .format(pri))
    plotms(vis=visname,
           field=pri,
           xaxis="channel",
           yaxis="amp",
           correlation="xy,yx",
           ydatacolumn="data",
           plotfile="{0}_{1}_{2}_ampvschan_postRFI_flag.png".format(
               epoch, ATCA_band, pri),
           showgui=False,
           overwrite=True)
    plotms(vis=visname,
           field=sec,
           xaxis="channel",
           yaxis="amp",
           correlation="xy,yx",
           ydatacolumn="data",
           plotfile="{0}_{1}_{2}_ampvschan_postRFI_flag.png".format(
               epoch, ATCA_band, sec),
           showgui=False,
           overwrite=True)
    plotms(vis=visname,
           field=tar,
           xaxis="channel",
           yaxis="amp",
           correlation="xy,yx",
           ydatacolumn="data",
           plotfile="{0}_{1}_{2}_ampvschan_postRFI_flag.png".format(
               epoch, ATCA_band, tar_nm),
           showgui=False,
           overwrite=True)
    print("+ + + + + + + + + + + + + + + +")
    print("+ Automated Flagging complete +")
    print("+ + + + + + + + + + + + + + + +")
else:
    print("Skipping...")

pipeline_options = raw_input("Split the target to make its own ms? ")
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    print("+ + + + + + + + + + + + + + + + +")
    print("+ Transforming measurement set  +")
    print("+ + + + + + + + + + + + + + + + +")
    os.system("rm -r {0}".format(msname))
    os.system("rm -r {0}.flagversions".format(msname))
    os.system("rm -r *.last")
    mstransform(vis=visname,
                outputvis=msname,
                regridms=True,
                datacolumn="data",
                mode="channel",
                spw=str(if_centre),
                nspw=n_spw,
                field="{0},{1},{2}".format(str(pri_cal_no), str(sec_cal_no),
                                           str(tar_cal_no)))
    listobs(vis=msname,
            listfile="listobs_{0}_{1}_{2}.dat".format(epoch, ATCA_band,
                                                      tar_nm),
            overwrite=True)
    flagmanager(vis=msname, mode="save", versionname="after_transform")
    plotms(vis=msname,
           field=pri,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="data",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_pre_cal.png".format(
               epoch, ATCA_band, pri),
           showgui=False,
           overwrite=True)
    plotms(vis=msname,
           field=sec,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="data",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_pre_cal.png".format(
               epoch, ATCA_band, sec),
           showgui=False,
           overwrite=True)

pipeline_options = raw_input("Begin calibration ? ")
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    print("+ + + + + + + + + + + + + + + + + +")
    print("+ Performing initial calibration  +")
    print("+ + + + + + + + + + + + + + + + + +")
    setjy(vis=msname,
          field=pri,
          scalebychan=True,
          standard="Perley-Butler 2010",
          usescratch=True)
    print("Performing gain calibration on {0}".format(pri))
    gaincal(vis=msname,
            caltable="cal_{0}_{1}.G0".format(pri, ATCA_band),
            field=pri,
            refant=ref,
            gaintype="G",
            calmode="p",
            parang=True,
            solint="60s")
    print("Performing bandpass calibration on {0}".format(pri))
    bandpass(vis=msname,
             caltable="cal_{0}_{1}.B0".format(pri, ATCA_band),
             field=pri,
             spw="",
             refant=ref,
             solnorm=True,
             solint="inf",
             bandtype="B",
             gaintable=["cal_{0}_{1}.G0".format(pri, ATCA_band)],
             parang=True)
    print("Determining gains on {0}".format(sec))
    gaincal(vis=msname,
            caltable="cal_{0}_{1}.G1".format(pri, ATCA_band),
            field=pri + "," + sec,
            refant=ref,
            spw="*",
            gaintype="G",
            calmode="ap",
            parang=True,
            solint="60s",
            gaintable=["cal_{0}_{1}.B0".format(pri, ATCA_band)])
    print("Solving for polarisation leakages.")
    print("+ + + + + + + + + + + + + + + + +")
    print("+ Final calibration & applying  +")
    print("+ + + + + + + + + + + + + + + + +")
    qu = qufromgain("cal_{0}_{1}.G1".format(pri, ATCA_band),
                    fieldids=[sec_cal_id
                              ])  # MAKE SURE THIS IS SET TO THE CORRECT FIELD
    smodel = [1, qu[sec_cal_id][0], qu[sec_cal_id][1], 0]
    print("smodel parameters = {0}".format(smodel))
    print(
        "Repeating entire calibration using best estimates for gains, leakages and secondary polarisation"
    )
    bandpass(vis=msname,
             caltable="cal_{0}_{1}.B1".format(pri, ATCA_band),
             field=pri,
             spw="",
             refant=ref,
             solnorm=True,
             solint="inf",
             bandtype="B",
             gaintable=["cal_{0}_{1}.G1".format(pri, ATCA_band)],
             parang=True)
    print("Deriving gain calibration using {0}".format(pri))
    gaincal(vis=msname,
            caltable="cal_{0}_{1}.G2".format(pri, ATCA_band),
            field=pri,
            refant=ref,
            spw="*",
            gaintype="G",
            calmode="ap",
            parang=True,
            solint="60s",
            gaintable=["cal_{0}_{1}.B1".format(pri, ATCA_band)])
    print("Deriving gain calibration using {0}".format(sec))
    gaincal(vis=msname,
            caltable="cal_{0}_{1}.G2".format(pri, ATCA_band),
            field=sec,
            refant=ref,
            spw="*",
            gaintype="G",
            calmode="ap",
            parang=True,
            solint="60s",
            gaintable=["cal_{0}_{1}.B1".format(pri, ATCA_band)],
            smodel=smodel,
            append=True)
    print(
        "Correcting the flux scale using comparison between the primary and secondary calibrator."
    )
    flux_transfer = fluxscale(vis=msname,
                              caltable="cal_{0}_{1}.G2".format(pri, ATCA_band),
                              fluxtable="cal_{0}_{1}.F0".format(
                                  pri, ATCA_band),
                              reference=pri)
    flagmanager(vis=msname, mode="save", versionname="before_applycal")

else:
    print("Skipping...")

pipeline_options = raw_input("Apply calibration?")
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    print("Applying calibration tables to {0}".format(pri))
    applycal(vis=msname,
             gaintable=[
                 "cal_{0}_{1}.B1".format(pri, ATCA_band),
                 "cal_{0}_{1}.F0".format(pri, ATCA_band)
             ],
             gainfield=[str(pri_cal_id),
                        str(pri_cal_id),
                        str(pri_cal_id)],
             field="{0}".format(pri_cal_id),
             parang=True,
             flagbackup=False)
    applycal(vis=msname,
             gaintable=[
                 "cal_{0}_{1}.B1".format(pri, ATCA_band),
                 "cal_{0}_{1}.F0".format(pri, ATCA_band)
             ],
             gainfield=[str(pri_cal_id),
                        str(pri_cal_id),
                        str(sec_cal_id)],
             field="{0},{1}".format(sec_cal_id, tar_cal_id),
             parang=True,
             flagbackup=False)
    print("Applying calibration tables to {0}".format(sec))
else:
    print("Skipping...")

pipeline_options = raw_input("Inspect calibrated {0} & {1} data ? ".format(
    pri, sec))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    plotms(vis=msname,
           field=pri,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_post_cal.png".format(
               epoch, ATCA_band, pri),
           showgui=False,
           overwrite=True)
    plotms(vis=msname,
           field=sec,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_post_cal.png".format(
               epoch, ATCA_band, sec),
           showgui=False,
           overwrite=True)
    plotms(vis=msname,
           field=pri,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="baseline",
           plotfile="{0}_{1}_{2}_ampvsfreqbaseline_post_cal.png".format(
               epoch, ATCA_band, pri),
           showgui=False,
           overwrite=True)
else:
    print("Skipping...")

pipeline_options = raw_input(
    "Flag calibrated {0} & {1} data for RFI ? ".format(pri, sec))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    # add print statements here.
    flagmanager(vis=msname, mode="save", versionname="before_rflag")
    spw = "0~{0}".format(str(n_spw - 1))
    flagdata(vis=msname,
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
             flagbackup=False)
    flagdata(vis=msname,
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
             flagbackup=False)
    flagdata(vis=msname,
             mode="extend",
             field=pri + "," + sec,
             spw="4:0~338",
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
             ntime="9999999min")
    plotms(vis=msname,
           field=pri,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_post_RFIflag.png".format(
               epoch, ATCA_band, pri),
           showgui=False,
           overwrite=True)
    plotms(vis=msname,
           field=sec,
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_post_RFIflag.png".format(
               epoch, ATCA_band, sec),
           showgui=False,
           overwrite=True)
else:
    print("Skipping...")

pipeline_options = raw_input("RFI flag and inspect {0} data ?".format(tar))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    print("Applying calibration tables to {0}".format(tar))
    plotms(vis=msname,
           field=str(tar_cal_id),
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="data",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_pre_cal.png".format(
               epoch, ATCA_band, tar),
           showgui=False,
           overwrite=True)
    applycal(vis=msname,
             gaintable=[
                 "cal_{0}_{1}.B1".format(pri, ATCA_band),
                 "cal_{0}_{1}.F0".format(pri, ATCA_band)
             ],
             gainfield=[pri, pri, sec],
             field=str(tar_cal_id),
             parang=True,
             flagbackup=False)
    print("Inspecting {0} amp vs freq Before RFI flagging".format(tar))
    plotms(vis=msname,
           field=str(tar_cal_id),
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="spw",
           plotfile="{0}_{1}_{2}_ampvsfreq_post_cal.png".format(
               epoch, ATCA_band, tar),
           showgui=False,
           overwrite=True)
    flagdata(vis=msname,
             mode="rflag",
             field=str(tar_cal_id),
             spw=spw,
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
             flagbackup=False)
    print("Inspecting {0} amp vs freq AFTER RFI flagging".format(tar))
    plotms(vis=msname,
           field=str(tar_cal_id),
           xaxis="frequency",
           yaxis="amp",
           correlation="xx,yy",
           ydatacolumn="corrected",
           coloraxis="spw",
           plotfile="{0}_{1}_ampvsfreq_post_RFIflag.png".format(
               epoch, ATCA_band, tar),
           showgui=False,
           overwrite=True)
    plotms(vis=msname,
           field=str(tar_cal_id),
           xaxis="u",
           yaxis="v",
           plotfile="{0}_{1}_uv.png".format(epoch, tar),
           overwrite=True,
           showgui=False)
else:
    print("Skipping...")

pipeline_options = raw_input("Begin automated imaging of {0} ?".format(tar))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    print("+ + + + + +")
    print("+ Imaging +")
    print("+ + + + + +")
    print(
        "This round is just to get hte mask, its using the mfs so you have high S/N and then it doesnt need the rest, use the mask generated here for all future rounds"
    )
    # set up parameters for imaging.
    mode = "mfs"
    nterms = 2
    niter = 3000
    if ATCA_band == "L":
        imsize = 2240
    if ATCA_band == "C":
        imsize = 1120
    if ATCA_band == "X":
        imsize == 940
    cell = "0.1arcsec"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = True
    gain = 0.01
    threshold = "2e-2Jy"

    os.system('rm -r {0}*'.format(targetms))
    split(vis=msname, datacolumn='corrected', field=tar, outputvis=targetms)
    listobs(vis=targetms,
            listfile="listobs_{0}_{1}_{2}_{3}.dat".format(
                epoch, ATCA_band, tar_nm, "preimage"),
            overwrite=True)
    os.system("rm -r {0}_{1}_{2}_mfs*".format(tar_nm, epoch, ATCA_band))
    imagename = "{0}_{1}_{2}_mfs".format(tar_nm, epoch, ATCA_band)
    spw = "0~{0}".format(str(n_spw - 1))
    flagmanager(vis=targetms, mode="save", versionname="before_selfcal")
    print("Initiating interactive cleaning on {0}".format(imagename))

    tclean(vis=targetms,
           imagename=imagename,
           selectdata=True,
           deconvolver="mtmfs",
           gain=gain,
           spw=spw,
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
           pbcor=False)
    tclean(vis=targetms,
           imagename=imagename,
           selectdata=True,
           deconvolver="mtmfs",
           gain=gain,
           spw=spw,
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
           calcpsf=False)

    pipeline_options = raw_input(
        "Begin first round of imaging for {0} ?".format(tar))
    if pipeline_options in ["Y", "y", "Yes", "yes"]:
        os.system("rm -r {0}_{1}_{2}_preself*".format(tar_nm, epoch, ATCA_band,
                                                      str(i)))
        mode = "mfs"
        nterms = 2
        niter = 2000
        threshold = "5e-3Jy"
        cell = "0.1arcsec"
        stokes = "I"
        weighting = "briggs"
        robust = 0.5
        interactive = False
        gain = 0.01
        flagmanager(vis=msname, mode="save", versionname="preself")

        for i in range(0, n_spw):
            spw = str(i)
            print("Cleaning on band: " + str(spw))
            imagename = "{0}_{1}_{2}_{3}".format(tar_nm, epoch, ATCA_band,
                                                 str(i))
            tclean(vis=targetms,
                   imagename=imagename + "_preself",
                   selectdata=True,
                   mask="{0}_{1}_{2}_mfs.mask".format(tar_nm, epoch,
                                                      ATCA_band),
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
                   pbcor=False)
            tclean(vis=targetms,
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
                   calcpsf=False)

    pipeline_options = raw_input(
        "Begin first round of selfcal for {0} ?".format(tar))
    if pipeline_options in ["Y", "y", "Yes", "yes"]:
        os.system("rm -r {0}_{1}_{2}_self1*".format(tar_nm, epoch, ATCA_band,
                                                    str(i)))
        mode = "mfs"
        nterms = 2
        niter = 2000
        threshold = "5e-4Jy"
        cell = "0.1arcsec"
        stokes = "I"
        weighting = "briggs"
        robust = 0.5
        interactive = False
        gain = 0.01

        # for i in range(0,n_spw):
        rmtables("pcal1")
        gaincal(vis=targetms,
                caltable="pcal1",
                gaintype="G",
                calmode="p",
                solint="60s",
                minsnr=3.0)
        applycal(vis=targetms,
                 gaintable="pcal1",
                 parang=True,
                 flagbackup=False)
        flagmanager(vis=targetms, mode="save", versionname="post self1")

        for i in range(0, n_spw):
            spw = str(i)
            print("Cleaning on band: " + str(spw))
            tclean(vis=targetms,
                   imagename=imagename + "_self1",
                   selectdata=True,
                   mask="{0}_{1}_{2}_mfs.mask".format(tar_nm, epoch,
                                                      ATCA_band),
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
                   pbcor=False)
            tclean(vis=targetms,
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
                   calcpsf=False)

    pipeline_options = raw_input(
        "Begin second round of selfcal for {0} ?".format(tar))
    if pipeline_options in ["Y", "y", "Yes", "yes"]:
        os.system("rm -r {0}_{1}_{2}_self2*".format(tar_nm, epoch, ATCA_band,
                                                    str(i)))
        mode = "mfs"
        nterms = 2
        niter = 2000
        threshold = "5e-5Jy"
        cell = "0.1arcsec"
        stokes = "I"
        weighting = "briggs"
        robust = 0.5
        interactive = False
        gain = 0.01
        rmtables("pcal2")
        gaincal(vis=targetms,
                caltable="pcal2",
                gaintable="pcal1",
                gaintype="G",
                calmode="p",
                solint="60s",
                minsnr=3.0)
        applycal(vis=targetms,
                 gaintable=["pcal1", "pcal2"],
                 parang=True,
                 flagbackup=False)
        flagmanager(vis=targetms, mode="save", versionname="post self2")
        for i in range(0, n_spw):
            spw = str(i)
            print("Cleaning on band: " + str(spw))
            tclean(vis=targetms,
                   imagename=imagename + "_self2",
                   selectdata=True,
                   mask="{0}_{1}_{2}_mfs.mask".format(tar_nm, epoch,
                                                      ATCA_band),
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
                   pbcor=False)
            tclean(vis=targetms,
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
                   calcpsf=False)

    pipeline_options = raw_input(
        "Begin third round of selfcal for {0} ?".format(tar))
    if pipeline_options in ["Y", "y", "Yes", "yes"]:

        os.system("rm -r {0}_{1}_{2}_self3*".format(tar_nm, epoch, ATCA_band,
                                                    str(i)))
        mode = "mfs"
        nterms = 2
        niter = 2000
        threshold = "5e-6Jy"
        cell = "0.1arcsec"
        stokes = "I"
        weighting = "briggs"
        robust = 0.5
        interactive = False
        rmtables("pcal3")
        gain = 0.01
        gaincal(vis=targetms,
                caltable="pcal3",
                gaintable=["pcal1", "pcal2"],
                gaintype="G",
                calmode="p",
                solint="60s",
                minsnr=3.0)
        applycal(vis=targetms,
                 gaintable=["pcal1", "pcal2", "pcal3"],
                 parang=True,
                 flagbackup=False)
        flagmanager(vis=targetms, mode="save", versionname="post self3")
        for i in range(0, n_spw):
            spw = str(i)
            print("Cleaning on band: " + str(spw))
            tclean(vis=targetms,
                   imagename=imagename + "_self3",
                   selectdata=True,
                   mask="{0}_{1}_{2}_mfs.mask".format(tar_nm, epoch,
                                                      ATCA_band),
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
                   pbcor=False)
            tclean(vis=targetms,
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
                   calcpsf=False)

pipeline_options = raw_input("PBcorr of {0}?".format(tar))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    mode = "mfs"
    nterms = 2
    niter = 2000
    threshold = "5e-5Jy"
    cell = "0.1arcsec"
    stokes = "I"
    weighting = "briggs"
    robust = 0.5
    interactive = False
    gain = 0.01
    for i in range(0, n_spw):
        os.system("rm -r " + imagename + "_pbcor")
        impbcor(imagename=imagename + "_self3.image.tt0",
                pbimage=imagename + "_self3.pb.tt0",
                outfile=imagename + "_self3_pbcor",
                cutoff=0.1,
                overwrite=True)

pipeline_options = raw_input(
    "Measure flux from visibilities of {0}?".format(tar))
if pipeline_options in ["Y", "y", "Yes", "yes"]:
    split(vis=targetms, datacolumn='corrected', outputvis=tar_ms)
    int_flux_c = []
    err_int_flux_c = []
    for i in range(n_spw):
        spw = str(i)
        # If things look like theyre not working, then check the source position! Chances are it can't find the source too far away from the phase centre
        outfile = '{0}_{1}_{2}_{3}.cl'.format(tar_nm, ATCA_band, epoch, spw)
        uvmodelfit(vis=tar_ms,
                   niter=10,
                   comptype='P',
                   spw=spw,
                   sourcepar=sourcepar,
                   outfile=outfile)
        cl.open(outfile)
        flux = cl.getfluxvalue(0)[0]
        flx_err = cl.getfluxerror(0)[0]
        int_flux_c.append(flux)
        err_int_flux_c.append(flx_err)
    if ATCA_band == 'C':
        t = Table()
        t['S_Cband'] = int_flux_c
        t['err_S_Cband'] = err_int_flux_c
        votable.writeto(t,
                        '{0}_{1}_{2}.votable'.format(tar_nm, epoch, ATCA_band))
        print(int_flux_c)
    elif ATCA_band == 'X':
        t = Table()
        t['S_Xband'] = int_flux_c
        t['err_S_Xband'] = err_int_flux_c
        votable.writeto(t,
                        '{0}_{1}_{2}.votable'.format(tar_nm, epoch, ATCA_band))
        print(int_flux_c)

    elif ATCA_band == 'L':
        t = Table()
        int_flux_l = int_flux_c[::-1]
        err_int_flux_l = err_int_flux_c[::-1]
        t['S_Lband'] = int_flux_l
        t['err_S_Lband'] = err_int_flux_c
        votable.writeto(t,
                        '{0}_{1}_{2}.votable'.format(tar_nm, epoch, ATCA_band))
        print(int_flux_l)

        os.system("mv *.png {1}".format(process_dir))
        os.system("mv *.fits {1}".format(img_dir))
        os.system("mv *.votable {1}".format(src_dir))
