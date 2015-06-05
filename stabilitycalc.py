#!/usr/bin/env python

import os
from os.path import join as pjoin
import time
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from htmltagutils import *
import shutil
import nibabel as nib

from mako.template import Template
from mako.lookup import TemplateLookup
makolookup = TemplateLookup(directories=['./tpl'])

import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
import stabilityfuncs as sf
import studyinfo


def stabilitycalc(dirname, filename, starttime, initxcenter=None, initycenter=None, initzcenter=None):
    """create the stability report for a scan"""

    isindividualcoil = False

    if initxcenter is not None:
        initxcenter, initycenter, initzcenter = [float(x) for x in (initxcenter, initycenter, initzcenter)]
        isindividualcoil = True

    nim = nib.load(pjoin(dirname, filename))
    nim_hdr = nim.get_header()
    xdim, ydim, slicethickness, tr = nim_hdr['pixdim'][1:5]
    xsize, ysize, numslices = nim_hdr['dim'][1:4]

    dims = nim.get_data().shape
    if len(dims) == 4:
        # normal case
        selecteddata = nim.get_data()[:, :, :, starttime:].transpose(3, 2, 1, 0)
    elif len(dims) == 3:
        # single slice
        selecteddata = nim.get_data()[:, :, starttime:]
        # nifti loses the z dimension; we must reconstitute it
        sdshape = selecteddata.shape
        selecteddata = selecteddata.reshape(sdshape[0], sdshape[1], 1, sdshape[2]).transpose(3, 2, 1, 0)

    numtimepoints = selecteddata.shape[0]

    info = studyinfo.extract_studyinfo(pjoin(dirname, filename))

    if info['ElementName'] == '':
        info['ElementName'] = 'unknown'
    if info['Coil'] != '':
        # TODO use a real datetime
        year = info['StudyDate'][0:4]
        month = info['StudyDate'][4:6]
        day = info['StudyDate'][6:8]
        hour = info['StudyTime'][0:2]
        minute = info['StudyTime'][2:4]
        second = info['StudyTime'][4:6]
        datetime = info['StudyDate'] + "T" + hour + minute + second
        formatteddate = month + "/" + day + "/" + year
        formattedtime = hour + ":" + minute + ":" + second

    # make empty results directory
    logging.debug("initializing output directory...")
    shutil.rmtree(pjoin(dirname, 'procresults'), ignore_errors=True)
    os.mkdir(pjoin(dirname, 'procresults'))

    # initialize the output files
    outputfile = pjoin(dirname, "procresults/output.html")
    thisdate = time.strftime("%m/%d/%Y %H:%M:%S", time.localtime())
    outfp = open(outputfile, "w")

    outfp.writelines(
        "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n\n")

    outfp.writelines("<head>\n")
    outfp.writelines("<title>Stability report for " + info['Coil'] + " generated on " + thisdate + "</title>\n")
    outfp.writelines("<style type=\"text/css\">\n")
    outfp.writelines("h1 {font-family:courier new;text-decoration:underline;}\n")
    outfp.writelines("h2 {font-family:courier new;color: teal; text-decoration:underline;}\n")
    outfp.writelines("h3 {font-family:courier new;color: maroon; text-decoration:none;}\n")
    outfp.writelines("h4 {font-family:courier new;text-decoration:none;}\n")
    outfp.writelines("p {font-family:courier new;color:black; font-size:16px;text-decoration:none;}\n")
    outfp.writelines("td {font-family:courier new;color:black; font-size:12px;text-decoration:none;}\n")
    outfp.writelines("</style>\n")
    outfp.writelines("</head>\n\n")

    outfp.writelines("<body>\n")
    outfp.writelines("<h2>Stability analysis of " + filename + "</h2>\n")
    if info['Coil'] != '':
        row1str = tablerowtag(bigtableentrytag("Coil:") + bigtableentrytag(info['Coil']))
        row2str = tablerowtag(bigtableentrytag("Element:") + bigtableentrytag(info['ElementName']))
        row3str = tablerowtag(bigtableentrytag("Date:") + bigtableentrytag(formatteddate))
        row4str = tablerowtag(bigtableentrytag("Time:") + bigtableentrytag(formattedtime))
        outfp.writelines(tablepropstag(row1str + row2str + row3str + row4str, 700, "left"))

    #############################
    #
    #  Calculate various statistical images
    #

    # calculate the mean, stddev, variance and ptp images
    logging.debug("calculating mean, stddev, and variance...")
    meanslice = np.mean(selecteddata, 0)
    stddevslice = np.std(selecteddata, 0)
    varslice = np.var(selecteddata, 0)
    ppslice = np.ptp(selecteddata, 0)

    # calculate a mask from meanimage and find its center
    (threshguess, threshfrac) = sf.findsepval(meanslice)
    initmask = sf.makemask(meanslice, threshfrac, 0)
    threshmean = sf.getnzfracval(initmask * meanslice, 0.02)
    if not isindividualcoil:
        objectmask = sf.makemask(meanslice, threshmean, 1)
    else:
        objectmask = sf.makemask(meanslice, 0.01 * threshmean, 1)

    logging.debug("calculating normalized standard deviation and sfnr...")
    with np.errstate(invalid='ignore'):
        normstdslice = objectmask * np.nan_to_num(100.0 * stddevslice / meanslice)
        minstddev = sf.nzrobust(objectmask * meanslice)[1] / 5000.0
        sfnrslice = np.where(objectmask * stddevslice > minstddev, meanslice / stddevslice, 0.0)

    # Now determine where the object is and how big it is
    slicecenter = sf.findCOM(objectmask)
    zcenterf = slicecenter[2]
    zcenter = int(round(zcenterf))
    slicecenter = sf.findCOM(objectmask[zcenter, :, :])
    xcenterf = slicecenter[0]
    ycenterf = slicecenter[1]
    xcenter, ycenter, zcenter = [int(round(x)) for x in xcenterf, ycenterf, zcenterf]
    xvec = objectmask[zcenter, ycenter, :]
    yvec = objectmask[zcenter, :, xcenter]
    # noinspection PyUnresolvedReferences
    xmin, xmax = np.nonzero(xvec)[0][[0, -1]]
    xcenterf = (xmin + xmax) / 2
    objectradiusx = (xmax - xmin + 1.0) / 2.0
    objectradiusx_mm = xdim * objectradiusx
    # noinspection PyUnresolvedReferences
    ymin, ymax = np.nonzero(yvec)[0][[0, -1]]
    ycenterf = (ymin + ymax) / 2
    objectradiusy = (ymax - ymin + 1.0) / 2.0
    objectradiusy_mm = ydim * objectradiusy
    xcenter, ycenter = [int(round(x)) for x in xcenterf, ycenterf]
    origslicecenter = (xcenterf, ycenterf, zcenterf)

    if isindividualcoil:
        # reset everything to assumed values
        # TODO get this from config
        objectmask[:, :, :] = 0.0
        objectradiusx_mm = 85.0
        objectradiusy_mm = 85.0
        objectradiusx, objectradiusy = objectradiusx_mm / xdim, objectradiusy_mm / ydim
        xcenterf, ycenterf, zcenterf = initxcenter, initycenter, initzcenter
        for i in range(xsize):
            ival = (float(i) - xcenterf) * xdim
            isq = ival * ival
            for j in range(ysize):
                jval = (float(j) - xcenterf) * ydim
                jsq = jval * jval
                for k in range(numslices):
                    kval = (float(k) - zcenterf) * slicethickness
                    ksq = kval * kval
                    if np.sqrt(isq + jsq + ksq) <= objectradiusx_mm:
                        objectmask[k, j, i] = 1.0
        xcenter, ycenter, zcenter = [int(round(x)) for x in xcenterf, ycenterf, zcenterf]

    # define the canonical limits
    limits = sf.getlimits(info['Coil'])

    # Try to figure out what we're looking at
    unknown = 0
    birn_phantom = 1
    head = 3
    objecttype = unknown
    objectname = "Unknown"

    object_radius_mm = np.sqrt(objectradiusx_mm * objectradiusy_mm)
    object_shape = objectradiusy / objectradiusx

    birn_phantom_radiuscheck = sf.limitcheck(object_radius_mm, limits['BIRNphantom_rad'])
    birn_phantom_shapecheck = sf.limitcheck(object_shape, limits['BIRNphantom_shape'])
    if (birn_phantom_radiuscheck < 2) and (birn_phantom_shapecheck < 2):
        objecttype = birn_phantom
        objectname = "BIRN_phantom"
        logging.debug("setting objecttype to birn_phantom")

    head_radiuscheck = sf.limitcheck(object_radius_mm, limits['head_rad'])
    head_shapecheck = sf.limitcheck(object_shape, limits['head_shape'])
    if (head_radiuscheck < 2) or (head_shapecheck < 2):
        objecttype = head
        objectname = "Head"
        logging.debug("setting objecttype to head")

    is_birn_sequence = True
    is_birn_protocol = False
    if (xsize != 64) or (ysize != 64) or (numslices != 28) or (tr != 2.0):
        is_birn_sequence = False
    if is_birn_sequence and (objecttype == birn_phantom):
        logging.debug("Assuming this is a BIRN protocol")
        is_birn_protocol = True
        protocolname = "fBIRN"
    else:
        logging.debug("Assuming this is NOT a BIRN protocol")
        protocolname = "Unknown"

    #############################
    #
    #       Odd-even SNR - Modified to match BIRN
    #
    logging.debug("calculating even/odd snr...")
    evenims = selecteddata[0:-1:2]
    oddims = selecteddata[1:-1:2]
    evenlength = evenims.shape[0]
    oddlength = oddims.shape[0]
    if oddlength < evenlength:
        evenims = evenims[:evenlength - 1]
    eodiffimage = np.sum(oddims, 0) - np.sum(evenims, 0)
    with np.errstate(invalid='ignore', over='ignore', divide='ignore'):
        eodiffpcimage = 100.0 * np.nan_to_num(eodiffimage / (objectmask * meanslice))

    #############################
    #
    #       Weisskoff analysis - Modified to match BIRN
    #
    numrois = 21
    weissstddevs = np.zeros(numrois)
    weisscvs = np.zeros(numrois)
    roisizes = range(1, numrois + 1)
    roiareas = np.zeros(numrois)
    projstddevs = np.zeros(numrois)
    projcvs = np.zeros(numrois)
    timepoints = np.arange(0.0, tr * numtimepoints, tr)
    for i in roisizes:
        roi = sf.setroilims(xcenter, ycenter, i)
        timecourse = sf.getroimeantc(selecteddata, roi, zcenter)
        weisskoffcoffs = np.polyfit(timepoints, timecourse, 2)
        fittc = sf.trendgen(timepoints, weisskoffcoffs)
        detrendedweisskofftc = timecourse - fittc
        dtweissmean = np.mean(detrendedweisskofftc)
        roiareas[i - 1] = np.square(roisizes[i - 1])
        weissstddevs[i - 1] = np.nan_to_num(np.std(detrendedweisskofftc))
        weisscvs[i - 1] = weissstddevs[i - 1] / dtweissmean
        projstddevs[i - 1] = weissstddevs[0] / i
        projcvs[i - 1] = weisscvs[0] / i
    weissrdc = weisscvs[0] / weisscvs[numrois - 1]

    #############################
    #
    #       Image analysis
    #

    with np.errstate(invalid='ignore'):
        meanstats = sf.nzstats(meanslice * objectmask)
        stddevstats = sf.nzstats(stddevslice * objectmask)
        varstats = sf.nzstats(varslice * objectmask)
        sfnrstats = sf.nzstats(np.nan_to_num(sfnrslice * objectmask))
        normstdstats = sf.nzstats(normstdslice * objectmask)
        eodiffstats = sf.nzstats(eodiffimage * objectmask)
        eodiffpcstats = sf.nzstats(np.nan_to_num(eodiffpcimage * objectmask))
        ppstats = sf.nzstats(ppslice * objectmask)

        objectmax = np.max(objectmask)
        objectmin = np.min(objectmask)
        rawmeanmax = sf.completerobust(meanslice)[1]
        [meanmin, meanmax] = sf.nzrobust(meanslice * objectmask)
        [stddevmin, stddevmax] = sf.nzrobust(stddevslice * objectmask)
        [varmin, varmax] = sf.nzrobust(varslice * objectmask)
        [sfnrmin, sfnrmax] = sf.nzrobust(np.nan_to_num(sfnrslice * objectmask))
        [normstdmin, normstdmax] = sf.nzrobust(normstdslice * objectmask)
        [eodiffmin, eodiffmax] = sf.nzrobust(eodiffimage * objectmask)
        [eodiffpcmin, eodiffpcmax] = sf.nzrobust(np.nan_to_num(eodiffpcimage * objectmask))
        [ppmin, ppmax] = sf.nzrobust(ppslice * objectmask)

    #############################
    #
    #       Corner (noise region) analysis
    #
    roislice = 0.5 * meanslice
    cornerroisize = 5
    cornerxpos = int(int(cornerroisize / 2.0) + 1)
    cornerypos = int(int(cornerroisize / 2.0) + 1)
    cornerroi = sf.setroilims(cornerxpos, cornerypos, cornerroisize)
    if not isindividualcoil:
        sf.markroi(cornerroi, zcenter, roislice, 0.91 * rawmeanmax)
    cornertc = sf.getroistdtc(selecteddata, cornerroi, zcenter)

    #############################
    #
    #       Central ROI analysis
    #
    logging.debug("Analyzing central ROI...")
    centralroisize = 10
    centralroi = sf.setroilims(xcenter, ycenter, centralroisize)
    if not isindividualcoil:
        sf.markroi(centralroi, zcenter, roislice, 0.92 * rawmeanmax)
    centtc = sf.getroimeantc(selecteddata, centralroi, zcenter)
    centsnrvec = centtc / cornertc
    centsnr = np.mean(centsnrvec)
    centsfnr = sf.getroival(sfnrslice, centralroi, zcenter)
    timepoints = np.arange(0.0, tr * numtimepoints, tr)
    centfitcoffs = np.polyfit(timepoints, timecourse, 2)
    fittc = sf.trendgen(timepoints, centfitcoffs)
    detrendedcenttc = centtc - fittc

    centmean = np.mean(centtc)
    centdrift = 100.0 * (np.max(fittc) - np.min(fittc)) / centmean
    centstddev = np.std(centtc)
    centmin = np.min(centtc)
    centmax = np.max(centtc)
    centpp = np.ptp(centtc)
    centmeanstr = "mean=%4.4f" % float(centmean)
    centdriftstr = "drift=%4.4f" % float(centdrift)

    stddevquality = sf.limitcheck(centstddev / centmean * 100.0, limits['central_roi_raw_std%'])
    centstddev_qualitytag = qualitytag("(%4.4f%%)", stddevquality) % (centstddev / centmean * 100.0)
    centstddevstr = "stddev=%4.4f " % centstddev + boldtag(centstddev_qualitytag)
    ppquality = sf.limitcheck(centpp / centmean * 100.0, limits['central_roi_raw_p-p%'])
    centpp_qualitytag = qualitytag("(%4.4f%%)", ppquality) % (centpp / centmean * 100.0)
    centppstr = "p-p=%4.4f " % centpp + boldtag(centpp_qualitytag)
    centtc_summary = centmeanstr + breaktag(centstddevstr) + breaktag(centppstr) + breaktag(centdriftstr)

    centmean_dt = np.mean(detrendedcenttc)
    centstddev_dt = np.std(detrendedcenttc)
    centmin_dt = np.min(detrendedcenttc)
    centmax_dt = np.max(detrendedcenttc)
    centpp_dt = np.ptp(detrendedcenttc)

    centstddev_dt_qualitytag = qualitytag("(%4.4f%%)", sf.limitcheck(centstddev_dt / centmean_dt * 100.0, limits['central_roi_detrended_std%']))
    centpp_dt_qualitytag = qualitytag("(%4.4f%%)", sf.limitcheck(centpp_dt / centmean_dt * 100.0, limits['central_roi_detrended_p-p%']))
    tcstats_format = "mean=%4.4f" + breaktag("stddev=%4.4f " + centstddev_dt_qualitytag) + breaktag("p-p=%4.4f " + centpp_dt_qualitytag)
    centtc_dt_summary = tcstats_format % (centmean_dt, centstddev_dt, centstddev_dt / centmean_dt * 100.0, centpp_dt, centpp_dt / centmean_dt * 100.0)

    centsnr_summary = "mean=%4.4f" % centsnr

    #############################
    #
    #       Maximum value ROI analysis
    #
    logging.debug("Finding and analyzing maximum signal ROI...")
    # TODO get from config
    maxlocroi_rawpplimits = ((0, 0.5), (0, 0.6))
    maxlocroi_dtpplimits = ((0, 0.5), (0, 0.6))
    maxlocroi_rawstddevlimits = ((0, 0.125), (0, 0.15))
    maxlocroi_dtstddevlimits = ((0, 0.125), (0, 0.15))
    maxlocroisize = 5
    maxlocradfrac = 0.7

    # find the maximum region
    if isindividualcoil:
        elementmask = sf.makemask(meanslice, threshmean, 1)
        elementcenter = sf.findCOM(elementmask)
        elementdirvec = (elementcenter[0] - origslicecenter[0],
                         elementcenter[1] - origslicecenter[1],
                         elementcenter[2] - origslicecenter[2])
        elementdirnormfac = sf.vecnorm(elementdirvec)
        maxlocoffsetscl = maxlocradfrac * objectradiusx / elementdirnormfac
        elementmaxpos = (origslicecenter[0] + maxlocoffsetscl * elementdirvec[0],
                         origslicecenter[1] + maxlocoffsetscl * elementdirvec[1],
                         origslicecenter[2] + maxlocoffsetscl * elementdirvec[2])
        if elementmaxpos[2] > numslices - 1:
            newmaxlocoffsetscl = (float(numslices - 1) - origslicecenter[2]) / elementdirvec[2]
            elementmaxpos = (origslicecenter[0] + newmaxlocoffsetscl * elementdirvec[0],
                             origslicecenter[1] + newmaxlocoffsetscl * elementdirvec[1],
                             origslicecenter[2] + newmaxlocoffsetscl * elementdirvec[2])
            logging.debug("maxpos adjusted to fall within valid image region")

        if elementmaxpos[2] < 0:
            newmaxlocoffsetscl = -origslicecenter[2] / elementdirvec[2]
            elementmaxpos = (origslicecenter[0] + newmaxlocoffsetscl * elementdirvec[0],
                             origslicecenter[1] + newmaxlocoffsetscl * elementdirvec[1],
                             origslicecenter[2] + newmaxlocoffsetscl * elementdirvec[2])
            logging.debug("maxpos adjusted to fall within valid image region")
        maxloccenterx, maxloccentery, maxloccenterz = [int(round(x)) for x in elementmaxpos[:3]]
        maxlocroi = sf.setroilims(maxloccenterx, maxloccentery, maxlocroisize)
        sf.markroi(maxlocroi, maxloccenterz, roislice, 0.92 * rawmeanmax)
        maxloctc = sf.getroimeantc(selecteddata, maxlocroi, zcenter)
        maxlocsnrvec = maxloctc / cornertc
        maxlocsnr = np.mean(maxlocsnrvec)
        maxlocsfnr = sf.getroival(sfnrslice, maxlocroi, zcenter)
        timepoints = np.arange(0.0, tr * numtimepoints, tr)
        maxlocfitcoffs = np.polyfit(timepoints, timecourse, 2)
        fittc = sf.trendgen(timepoints, maxlocfitcoffs)
        detrendedmaxloctc = maxloctc - fittc

        maxlocmean = np.mean(maxloctc)
        maxlocstddev = np.std(maxloctc)
        maxlocpp = np.ptp(maxloctc)
        maxlocmeanstr = "mean=%4.4f" % maxlocmean
        stddevquality = sf.limitcheck(maxlocstddev / maxlocmean * 100.0, maxlocroi_rawstddevlimits)
        maxlocstddevstr = "stddev=%4.4f " % maxlocstddev + boldtag(
            qualitytag("(%4.4f%%)", stddevquality) % (maxlocstddev / maxlocmean * 100.0))
        ppquality = sf.limitcheck(maxlocpp / maxlocmean * 100.0, maxlocroi_rawpplimits)
        maxlocppstr = "p-p=%4.4f " % maxlocpp + boldtag(
            qualitytag("(%4.4f%%)", ppquality) % (maxlocpp / maxlocmean * 100.0))
        maxloctc_summary = maxlocmeanstr + breaktag(maxlocstddevstr) + breaktag(maxlocppstr)

        maxlocmean_dt = np.mean(detrendedmaxloctc)
        maxlocstddev_dt = np.std(detrendedmaxloctc)
        maxlocmin_dt = np.min(detrendedmaxloctc)
        maxlocmax_dt = np.max(detrendedmaxloctc)
        maxlocpp_dt = np.ptp(detrendedmaxloctc)
        stddevformat = qualitytag("(%4.4f%%)",
                                  sf.limitcheck(maxlocstddev_dt / maxlocmean_dt * 100.0, maxlocroi_dtstddevlimits))
        ppformat = qualitytag("(%4.4f%%)", sf.limitcheck(maxlocpp_dt / maxlocmean_dt * 100.0, maxlocroi_dtpplimits))
        tcstats_format = "mean=%4.4f" + breaktag("stddev=%4.4f " + stddevformat) + breaktag("p-p=%4.4f " + ppformat)
        maxloctc_dt_summary = tcstats_format % (
            maxlocmean_dt, maxlocstddev_dt, maxlocstddev_dt / maxlocmean_dt * 100.0, maxlocpp_dt,
            maxlocpp_dt / maxlocmean_dt * 100.0)

        maxlocsnr_summary = "mean=%4.4f" % maxlocsnr

    #############################
    #
    # Individual coil assessment ROIs
    #
    logging.debug("Analyzing phased array ROIs...")

    isphasedarray = False

    coildata = sf.getphasedarraydata(info['Coil'])
    if coildata and numslices > 1:
        isphasedarray = True
        numphasedarray = len(coildata)

    if isphasedarray:
        paindices = np.arange(0.0, numphasedarray, 1.0)
        # TODO move these into a config file
        phasedarrayroi_rawpplimits = {'good_max': 0.5, 'good_min': 0, 'warn_max': 0.6, 'warn_min': 0}
        phasedarrayroi_dtpplimits = {'good_max': 0.5, 'good_min': 0, 'warn_max': 0.6, 'warn_min': 0}
        phasedarrayroi_rawstddevlimits = {'good_max': 0.125, 'good_min': 0, 'warn_max': 0.15, 'warn_min': 0}
        phasedarrayroi_dtstddevlimits = {'good_max': 0.125, 'good_min': 0, 'warn_max': 0.15, 'warn_min': 0}
        phasedarraysize = 5

        phasedarrayroimeans = np.zeros(numphasedarray)
        phasedarrayroistddevs = np.zeros(numphasedarray)
        phasedarrayroimins = np.zeros(numphasedarray)
        phasedarrayroimaxs = np.zeros(numphasedarray)
        phasedarrayroipps = np.zeros(numphasedarray)
        phasedarrayroimeans_dt = np.zeros(numphasedarray)
        phasedarrayroistddevs_dt = np.zeros(numphasedarray)
        phasedarrayroimins_dt = np.zeros(numphasedarray)
        phasedarrayroimaxs_dt = np.zeros(numphasedarray)
        phasedarrayroipps_dt = np.zeros(numphasedarray)
        phasedarrayroisfnrs = np.zeros(numphasedarray)
        phasedarrayroisnrs = np.zeros(numphasedarray)
        phasedarraytcs = np.zeros((numphasedarray, len(timepoints)), dtype=float)
        phasedarrayfittcs = np.zeros((numphasedarray, len(timepoints)), dtype=float)
        phasedarraydttcs = np.zeros((numphasedarray, len(timepoints)), dtype=float)
        phasedarraydttcs_demeaned = np.zeros((numphasedarray, len(timepoints)), dtype=float)
        phasedarraytc_summary = []
        phasedarraytc_dt_summary = []
        for i, ele in enumerate(coildata):
            roi = sf.setroilims(round(coildata[ele]['xloc']), round(coildata[ele]['yloc']), phasedarraysize)
            if not isindividualcoil:
                sf.markroi(roi, round(coildata[ele]['zloc']), roislice, 0.95 * rawmeanmax)
            timecourse = sf.getroimeantc(selecteddata, roi, zcenter)
            phasedarraytcs[i, :] = timecourse[:]
            snrvec = sf.getroisnr(selecteddata, roi, round(coildata[ele]['zloc']))
            phasedarrayroisnrs[i] = np.mean(snrvec)
            phasedarrayroisfnrs[i] = sf.getroival(sfnrslice, roi, round(coildata[ele]['zloc']))
            phasedarrayroimeans[i] = np.mean(timecourse)
            phasedarrayroistddevs[i] = np.std(timecourse)
            phasedarrayroimins[i] = np.min(timecourse)
            phasedarrayroimaxs[i] = np.max(timecourse)
            phasedarrayroipps[i] = np.ptp(timecourse)
            phasedarrayfitcoffs = np.polyfit(timepoints, timecourse, 2)
            phasedarrayfittcs[i, :] = sf.trendgen(timepoints, phasedarrayfitcoffs)
            phasedarraydttcs[i, :] = phasedarraytcs[i, :] - phasedarrayfittcs[i, :]
            phasedarrayroimeans_dt[i] = np.mean(phasedarraydttcs[i, :])
            phasedarraydttcs_demeaned[i, :] = phasedarraydttcs[i, :] - phasedarrayroimeans_dt[i]
            phasedarrayroistddevs_dt[i] = np.std(phasedarraydttcs[i, :])
            phasedarrayroimins_dt[i] = np.min(phasedarraydttcs[i, :])
            phasedarrayroimaxs_dt[i] = np.max(phasedarraydttcs[i, :])
            phasedarrayroipps_dt[i] = np.ptp(phasedarraydttcs[i, :])

            # do average timecourse calculations
            phasedarraymeanstr = "mean=%4.4f" % phasedarrayroimeans[i]
            stddevquality = sf.limitcheck(phasedarrayroistddevs[i] / phasedarrayroimeans[i] * 100.0,
                                          phasedarrayroi_rawstddevlimits)
            phasedarraystddevstr = "stddev=%4.4f " % (phasedarrayroistddevs[i]) + boldtag(qualitytag("(%4.4f%%)",
                                                                                                     stddevquality) % (
                phasedarrayroistddevs[i] /
                phasedarrayroimeans[
                    i] * 100.0))
            ppquality = sf.limitcheck(phasedarrayroipps[i] / phasedarrayroimeans[i] * 100.0, phasedarrayroi_rawpplimits)
            phasedarrayppstr = "p-p=%4.4f " % (phasedarrayroipps[i]) + boldtag(
                qualitytag("(%4.4f%%)", ppquality) % (phasedarrayroipps[i] / phasedarrayroimeans[i] * 100.0))
            phasedarraytc_summary.append(
                phasedarraymeanstr + breaktag(phasedarraystddevstr) + breaktag(phasedarrayppstr))

            stddevformat = qualitytag("(%4.4f%%)",
                                      sf.limitcheck(phasedarrayroistddevs_dt[i] / phasedarrayroimeans_dt[i] * 100.0,
                                                    phasedarrayroi_dtstddevlimits))
            ppformat = qualitytag("(%4.4f%%)",
                                  sf.limitcheck(phasedarrayroipps_dt[i] / phasedarrayroimeans_dt[i] * 100.0,
                                                phasedarrayroi_dtpplimits))
            tcstats_format = "mean=%4.4f" + breaktag("stddev=%4.4f " + stddevformat) + breaktag("p-p=%4.4f " + ppformat)
            phasedarraytc_dt_summary.append(tcstats_format % (phasedarrayroimeans_dt[i], phasedarrayroistddevs_dt[i],
                                                              phasedarrayroistddevs_dt[i] / phasedarrayroimeans_dt[
                                                                  i] * 100.0, phasedarrayroipps_dt[i],
                                                              phasedarrayroipps_dt[i] / phasedarrayroimeans_dt[
                                                                  i] * 100.0))

        phasedarrayroipps_percent = 100.0 * phasedarrayroipps / phasedarrayroimeans
        phasedarrayroipps_dt_percent = 100.0 * phasedarrayroipps_dt / phasedarrayroimeans_dt

        # finally calculate the correlation between the timeseries
        coilccmatrix = np.corrcoef(phasedarraydttcs_demeaned)

    #############################
    #
    # Peripheral ROIs
    #
    # TODO get these from config file
    peripheralroisize = 3
    peripheralradfrac = 0.8
    logging.debug("Analyzing peripheral ROIs...")
    voxperroi = peripheralroisize * peripheralroisize
    peripheralradiusx = peripheralradfrac * objectradiusx
    peripheralradiusy = peripheralradfrac * objectradiusy
    numperiph = 32

    periphindex = np.arange(numperiph)
    periphangles = (np.pi * 2 * periphindex) / numperiph
    periphanglesd = (360.0 * periphindex) / numperiph
    xlocs = xcenterf + peripheralradiusx * np.sin(periphangles)
    ylocs = ycenterf + peripheralradiusy * np.cos(periphangles)
    periphangmeans = np.zeros(numperiph)
    periphangsfnrs = np.zeros(numperiph)
    periphangsnrs = np.zeros(numperiph)
    avgperiphtc = centtc * 0.0
    avgperiphsnrvec = centtc * 0.0
    periphvoxels = np.zeros((selecteddata.shape[0], voxperroi * numperiph))
    for i in periphindex:
        roi = sf.setroilims(round(xlocs[i]), round(ylocs[i]), peripheralroisize)
        if not isindividualcoil:
            sf.markroi(roi, zcenter, roislice, 0.95 * rawmeanmax)
        timecourse = sf.getroimeantc(selecteddata, roi, zcenter)
        newvoxels = sf.getroivoxels(selecteddata, roi, zcenter)
        for j in range(voxperroi):
            periphvoxels[:, i * voxperroi + j] = newvoxels[:, j]
        snrvec = sf.getroisnr(selecteddata, roi, zcenter)
        snr = np.mean(snrvec)
        avgperiphtc += timecourse / (1.0 * numperiph)
        avgperiphsnrvec += snrvec / (1.0 * numperiph)
        sfnrval = sf.getroival(sfnrslice, roi, zcenter)
        periphangmeans[i] = 100.0 * np.mean(timecourse) / centmean
        periphangsfnrs[i] = sfnrval
        periphangsnrs[i] = snr

    avgperiphtc2 = np.mean(periphvoxels, 1)
    avgperiphsnrvec = avgperiphtc2 / cornertc

    # do average timecourse calculations
    periphfitcoffs = np.polyfit(timepoints, avgperiphtc, 2)
    periphfittc = sf.trendgen(timepoints, periphfitcoffs)
    detrendedperiphtc = avgperiphtc - periphfittc

    periphmean = np.mean(avgperiphtc)
    periphdrift = 100.0 * (np.max(periphfittc) - np.min(periphfittc)) / periphmean
    periphstddev = np.std(avgperiphtc)
    periphmin = np.min(avgperiphtc)
    periphmax = np.max(avgperiphtc)
    periphpp = np.ptp(avgperiphtc)
    periphmeanstr = "mean=%4.4f" % float(periphmean)
    periphdriftstr = "drift=%4.4f" % float(periphdrift)
    stddevquality = sf.limitcheck(periphstddev / periphmean * 100.0, limits['peripheral_roi_raw_std%'])
    periph_qualitytag = qualitytag("(%4.4f%%)", stddevquality) % (periphstddev / periphmean * 100.0)
    periphstddevstr = "stddev=%4.4f " % periphstddev + boldtag(periph_qualitytag)
    periphppquality = sf.limitcheck(periphpp / periphmean * 100.0, limits['peripheral_roi_raw_p-p%'])
    periphpp_qualitytag = qualitytag("(%4.4f%%)", periphppquality) % (periphpp / periphmean * 100.0)
    periphppstr = "p-p=%4.4f " % periphpp + boldtag(periphpp_qualitytag)
    periphtc_summary = periphmeanstr + breaktag(periphstddevstr) + breaktag(periphppstr) + breaktag(periphdriftstr)

    periphmean_dt = np.mean(detrendedperiphtc)
    periphstddev_dt = np.std(detrendedperiphtc)
    periphmin_dt = np.min(detrendedperiphtc)
    periphmax_dt = np.max(detrendedperiphtc)
    periphpp_dt = np.ptp(detrendedperiphtc)
    periph_dt_qualitytag = qualitytag("(%4.4f%%)", sf.limitcheck(periphstddev_dt / periphmean_dt * 100.0, limits['peripheral_roi_detrended_std%']))
    periphpp_dt_qualitytag = qualitytag("(%4.4f%%)", sf.limitcheck(periphpp_dt / periphmean_dt * 100.0, limits['peripheral_roi_detrended_p-p%']))
    tcstats_format = "mean=%4.4f" + breaktag("stddev=%4.4f " + periph_dt_qualitytag) + breaktag("p-p=%4.4f " + periphpp_dt_qualitytag)
    periphtc_dt_summary = tcstats_format % (
        periphmean_dt, periphstddev_dt, periphstddev_dt / periphmean_dt * 100.0, periphpp_dt,
        periphpp_dt / periphmean_dt * 100.0)

    # do calculations regarding angular dependance
    meanangperiphval = np.mean(periphangmeans)
    meanangperiphsfnr = np.mean(periphangsfnrs)
    meanangperiphsnr = np.mean(avgperiphsnrvec)

    ptpangperiphval = np.ptp(periphangmeans)
    periphangintensitymeanstr = "mean=%4.4f" % meanangperiphval
    periphangintensityppquality = sf.limitcheck(ptpangperiphval / meanangperiphval * 100.0, limits['peripheral_angle_p-p%'])
    periphangintensity_qualitytag = qualitytag("(%4.4f%%)", periphangintensityppquality) % (ptpangperiphval / meanangperiphval * 100.0)
    ptpangperiphvalstr = "p-p=%4.4f " % ptpangperiphval + boldtag(periphangintensity_qualitytag)
    periphangintensity_summary = periphangintensitymeanstr + breaktag(ptpangperiphvalstr)

    periphang_sfnr_ptp = np.ptp(periphangsfnrs)
    periphang_sfnr_meanstr = "mean=%4.4f" % meanangperiphsfnr
    periphang_sfnr_ppquality = sf.limitcheck((periphang_sfnr_ptp / meanangperiphsfnr) * 100.0,
                                             limits['peripheral_angle_SFNR_p-p%'])
    periphang_sfnr_pp_qualitytag = qualitytag("(%4.4f%%)", periphang_sfnr_ppquality) % ((periphang_sfnr_ptp / meanangperiphsfnr) * 100.0)
    periphang_sfnr_ptpstr = "p-p=%4.4f " % periphang_sfnr_ptp + boldtag(periphang_sfnr_pp_qualitytag)
    periphang_sfnr_summary = periphang_sfnr_meanstr + breaktag(periphang_sfnr_ptpstr)

    periphsnr_summary = "mean=%4.4f" % meanangperiphsnr

    #############################
    #
    #       Ghost ROI analysis
    #
    ghostroisize = 4
    ghostevenxpos = int(int(xsize / 2) + 1)
    ghostoddxpos = int(xcenterf - objectradiusx + ghostroisize - 1)
    ghostypos = int(1 + ghostroisize / 2.0)
    evenghostroi = sf.setroilims(ghostevenxpos, ghostypos, ghostroisize)
    oddghostroi = sf.setroilims(ghostoddxpos, ghostypos, ghostroisize)
    if not isindividualcoil:
        sf.markroi(evenghostroi, zcenter, roislice, 0.97 * rawmeanmax)
        sf.markroi(oddghostroi, zcenter, roislice, 0.97 * rawmeanmax)

    evenghosttc = sf.getroimeantc(selecteddata, evenghostroi, zcenter)
    oddghosttc = sf.getroimeantc(selecteddata, oddghostroi, zcenter)

    relevenghosttc = 100.0 * evenghosttc / centtc
    reloddghosttc = 100.0 * oddghosttc / centtc

    oddghostmean = np.mean(reloddghosttc)
    oddghoststddev = np.std(reloddghosttc)
    oddghostmin = np.min(reloddghosttc)
    oddghostmax = np.max(reloddghosttc)
    oddghostpp = np.ptp(reloddghosttc)
    evenghostmean = np.mean(relevenghosttc)
    evenghoststddev = np.std(relevenghosttc)
    evenghostmin = np.min(relevenghosttc)
    evenghostmax = np.max(relevenghosttc)
    evenghostpp = np.ptp(relevenghosttc)
    tcstats_format2 = "mean=%4.4f" + breaktag("stddev=%4.4f") + breaktag("min=%4.4f") + breaktag(
        "max=%4.4f") + breaktag(
        "p-p=%4.4f ")
    oddghosttc_summary = tcstats_format2 % (oddghostmean, oddghoststddev, oddghostmin, oddghostmax, oddghostpp)
    evenghosttc_summary = tcstats_format2 % (evenghostmean, evenghoststddev, evenghostmin, evenghostmax, evenghostpp)

    #############################
    #
    # Output
    #

    # sample type and protocol type description
    if objecttype == unknown:
        objecttypestr = "unknown"
    elif objecttype == head:
        objecttypestr = "head"
    elif objecttype == birn_phantom:
        objecttypestr = "BIRN phantom"

    row1str = tablerowtag(bigtableentrytag("Object center of mass:") + bigtableentrytag(
        str(xcenterf) + "," + str(ycenterf) + "," + str(zcenterf)))
    row2str = tablerowtag(bigtableentrytag("Object type:") + bigtableentrytag(objecttypestr))
    row3str = tablerowtag(bigtableentrytag("Object mean radius:") + bigtableentrytag(str(object_radius_mm)))
    row4str = tablerowtag(bigtableentrytag("Object shape factor:") + bigtableentrytag(str(object_shape)))
    outfp.writelines(tablepropstag(row1str + row2str + row3str + row4str, 700, "left"))
    if is_birn_protocol:
        outfp.writelines("<h3>Imaging protocol: BIRN stability</h3>\n")

    # statistical images section

    # noinspection PyShadowingNames
    def slicepic(inputslice, caption, minp, maxp, dirname, outputname, colormap):
        plt.figure(figsize=plt.figaspect(1.0))
        sf.showslice2(inputslice, caption, minp, maxp, colormap)
        plt.savefig(pjoin(dirname, 'procresults', outputname + '.png'), format='png')
        plt.close()

    roimin = np.min(roislice)
    roimax = np.max(roislice)

    with np.errstate(invalid='ignore'):
        slicepic(roislice, "ROI locations", roimin, roimax, dirname, 'roiimage', 1)
        slicepic(normstdslice, "Normalized stddev % image", normstdmin, normstdmax, dirname, 'normstdimage', 0)
        slicepic(objectmask, "Object mask", objectmin, objectmax, dirname, 'objectmaskimage', 0)
        slicepic(varslice, "Variance image", varmin, varmax, dirname, 'varimage', 0)
        slicepic(stddevslice, "Stddev image", stddevmin, stddevmax, dirname, 'stddevimage', 0)
        slicepic(meanslice, "Mean image", meanmin, meanmax, dirname, 'meanimage', 0)
        slicepic(sfnrslice, "SFNR image", sfnrmin, sfnrmax, dirname, 'sfnrimage', 0)
        slicepic(eodiffimage, "Even odd diff image", eodiffmin, eodiffmax, dirname, 'eodiffimage', 0)
        slicepic(np.nan_to_num(objectmask * eodiffpcimage), "Even odd diff percent image", eodiffpcmin, eodiffpcmax, dirname, 'eodiffpcimage', 0)
        slicepic(ppslice, "Peak to peak image", ppmin, ppmax, dirname, 'ppimage', 0)

    def makefig(figk):
        """make a figure using the figure dictionary"""
        plt.figure()
        # noinspection PyCallingNonCallable
        figs[figk]['fn'](*figs[figk]['args'])
        plt.savefig(pjoin(dirname, 'procresults', figk + '.png'), format='png')
        plt.close()

    # @formatter:off
    figs = {
        'weisskoffplot': {'fn': sf.showweisskoff, 'args': (roiareas, weissstddevs, projstddevs, "Weisskoff plot")},
        'oddghostroiplot': {'fn': sf.showtc2, 'args': (timepoints, reloddghosttc, 0.0 * timepoints + oddghostmean, "Relative odd ghost ROI amplitude plot (%)")},
        'evenghostroiplot': {'fn': sf.showtc2, 'args': (timepoints, relevenghosttc, 0.0 * timepoints + evenghostmean, "Relative even ghost ROI amplitude plot (%)")}}

    if isphasedarray:
        figs.update({
            'coilccmatrix': {'fn': sf.showslice3, 'args': (coilccmatrix, "Phased array element correlation matrix", 0.0, 1.0, 0)},
            'phasedarrayroisnrplot': {'fn': sf.showtc, 'args': (paindices, phasedarrayroisnrs, "Phased array SNR by element")},
            'phasedarrayroisfnrplot': {'fn': sf.showtc, 'args': (paindices, phasedarrayroisfnrs, "Phased array SFNR by element")},
            'phasedarrayroippplot': {'fn': sf.showtc, 'args': (paindices, phasedarrayroipps_percent, "Phased array p-p% variation by element")},
            'phasedarrayroippdtplot': {'fn': sf.showtc, 'args': (paindices, phasedarrayroipps_dt_percent, "Phased array p-p% variation by element (after detrending)")}})

    if isindividualcoil:
        figs.update({
            'maxlocroiplot': {'fn': sf.showtc2, 'args': (timepoints, maxloctc, maxlocmean_dt + fittc, "Max sensitivity ROI plot (%)")},
            'maxlocroisnrplot': {'fn': sf.showtc2, 'args': (timepoints, maxlocsnrvec, 0.0 * maxlocsnrvec + maxlocsnr, "Max sensitivity ROI SNR over time")},
            'maxlocroidtplot': {'fn': sf.showtc, 'args': (timepoints, maxloctc - fittc, "Detrended max sensitivity ROI plot (%)")}})
    else:
        figs.update({
            'centroiplot': {'fn': sf.showtc2, 'args': (timepoints, centtc, centmean_dt + fittc, "Central ROI plot (%)")},
            'centroisnrplot': {'fn': sf.showtc2, 'args': (timepoints, centsnrvec, 0.0 * centsnrvec + centsnr, "Central ROI SNR over time")},
            'centroidtplot': {'fn': sf.showtc, 'args': (timepoints, centtc - fittc, "Detrended central ROI plot (%)")},
            'periphroiplot': {'fn': sf.showtc2, 'args': (periphanglesd, periphangmeans, 0.0 * periphanglesd + meanangperiphval, "Relative peripheral image intensity (%)")},
            'periphroisfnrplot': {'fn': sf.showtc2, 'args': (periphanglesd, periphangsfnrs, 0.0 * periphanglesd + meanangperiphsfnr, "Absolute peripheral SFNR")},
            'periphroitcplot': {'fn': sf.showtc2, 'args': (timepoints, avgperiphtc, periphmean_dt + periphfittc, "Peripheral ROI plot (%)")},
            'periphroisnrplot': {'fn': sf.showtc2, 'args': (timepoints, avgperiphsnrvec, 0.0 * timepoints + meanangperiphsnr, "Peripheral ROI SNR over time")},
            'periphroidttcplot': {'fn': sf.showtc, 'args': (timepoints, avgperiphtc - periphfittc, "Detrended peripheral ROI plot (%)")}})

    # @formatter:on

    for k in figs:
        makefig(k)

    ########################################################
    #
    # Compose the image table
    #
    calcimagehdrstr = hruletag() + headertag("Calculated images")
    myimwidth = 500  # TODO get from config

    meanimagestr = sf.makecaptionedimage("Mean over time:", meanstats, "meanimage.png", myimwidth)
    stddevimagestr = sf.makecaptionedimage("Standard deviation over time:", stddevstats, "stddevimage.png", myimwidth)
    varianceimagestr = sf.makecaptionedimage("Variance over time:", varstats, "varimage.png", myimwidth)
    normstdimagestr = sf.makecaptionedimage("Normalized % stddev over time:", normstdstats, "normstdimage.png",
                                            myimwidth)
    sfnrimagestr = sf.makecaptionedimage("SFNR:", sfnrstats, "sfnrimage.png", myimwidth)
    eoimagestr = sf.makecaptionedimage("Even-odd difference:", eodiffstats, "eodiffimage.png", myimwidth)
    eopcimagestr = sf.makecaptionedimage("Even-odd difference percent:", eodiffpcstats, "eodiffpcimage.png", myimwidth)
    ppimagestr = sf.makecaptionedimage("Peak to peak:", ppstats, "ppimage.png", myimwidth)
    objectstr = sf.makecaptionedimage("Object mask", [], "objectmaskimage.png", myimwidth)
    roistr = sf.makecaptionedimage("ROI locations:", [], "roiimage.png", myimwidth)

    row1str = tablerowtag(tableentrytag(meanimagestr) + tableentrytag(stddevimagestr))
    row2str = tablerowtag(tableentrytag(varianceimagestr) + tableentrytag(sfnrimagestr))
    row3str = tablerowtag(tableentrytag(eoimagestr) + tableentrytag(eopcimagestr))
    row4str = tablerowtag(tableentrytag(ppimagestr) + tableentrytag(normstdimagestr))
    row5str = tablerowtag(tableentrytag(roistr) + tableentrytag(objectstr))
    outfp.writelines(tablepropstag(calcimagehdrstr + row1str + row2str + row3str + row4str + row5str, 500, "left"))

    ########################################################
    #
    # Central ROI output
    #
    centralroihdrstr = hruletag() + headertag("Central ROI Analysis")
    croiparastr = paratag(
        "This is an analysis of the temporal fluctuation in a central voxels. A " + str(centralroisize) + "x" + str(
            centralroisize) +
        "x1 voxel is automatically positioned at the center of gravity, and the average voxel value is plotted a function of time. Linear and quadratic terms are fit to the drift. The voxel statistics are reported with and without the drift removed.\n")

    rawroicellstr = tableentrytag(
        paratag(boldtag("Raw central ROI plot")) +
        paratag(centtc_summary) +
        imagetag("centroiplot.png", myimwidth))
    snrroicellstr = tableentrytag(
        paratag(boldtag("Central ROI over time")) +
        paratag(centsnr_summary) +
        imagetag("centroisnrplot.png", myimwidth))
    detrendedroicellstr = tableentrytag(
        paratag(boldtag("Detrended central ROI plot")) +
        paratag(centtc_dt_summary) +
        imagetag("centroidtplot.png", myimwidth))
    row1str = tablerowtag(rawroicellstr + detrendedroicellstr)
    row2str = tablerowtag(snrroicellstr)
    if not isindividualcoil:
        outfp.writelines(tablepropstag(centralroihdrstr + croiparastr + row1str + row2str, 500, "left"))

    ########################################################
    #
    # Peripheral ROI output
    #
    peripheralroihdrstr = hruletag() + headertag("Peripheral ROI Analysis")
    periphroiparastr = paratag(
        "This is an analysis of the image intensity and SFNR variation in a set of " + str(
            peripheralroisize) + "x" + str(
            peripheralroisize) + "x1 voxels at a fixed radius from the center of the object in the central axial slice (" +
        str(
            peripheralradfrac) + " of the distance from the center to the edge of the phantom). For phased array coils this is likely to have better signal to noise than an ROI at the center of the coil.  The variation in SNR as a function of angle is also displayed.\n")

    rawperiphroicellstr = tableentrytag(
        paratag(boldtag("Raw peripheral ROI plot")) +
        paratag(periphtc_summary) +
        imagetag("periphroitcplot.png", myimwidth))
    detrendedperiphroicellstr = tableentrytag(
        paratag(boldtag("Detrended peripheral ROI plot")) +
        paratag(periphtc_dt_summary) +
        imagetag("periphroidttcplot.png", myimwidth))
    row1str = tablerowtag(rawperiphroicellstr + detrendedperiphroicellstr)
    periphroicellstr = tableentrytag(
        paratag(boldtag("Peripheral ROI intensity plot")) +
        paratag(periphangintensity_summary) +
        imagetag("periphroiplot.png", myimwidth))
    periphsfnrcellstr = tableentrytag(
        paratag(boldtag("Peripheral ROI SFNR plot")) +
        paratag(periphang_sfnr_summary) +
        imagetag("periphroisfnrplot.png", myimwidth))
    row2str = tablerowtag(periphroicellstr + periphsfnrcellstr)
    snrperiphroicellstr = tableentrytag(
        paratag(boldtag("Peripheral ROI SNR plot")) +
        paratag(periphsnr_summary) +
        imagetag("periphroisnrplot.png", myimwidth))
    row3str = tablerowtag(snrperiphroicellstr)
    if not isindividualcoil:
        outfp.writelines(
            tablepropstag(peripheralroihdrstr + periphroiparastr + row1str + row2str + row3str, 500, "left"))

    ########################################################
    #
    # Phased array ROI output
    #
    if isphasedarray:
        phasedarrayroihdrstr = hruletag() + headertag("Phased array element maximum sensitivity region ROI Analysis")
        phasedarrayroiparastr = paratag(
            "This is an analysis of the variation in SNR, SFNR, and stability parameters across the individual coil elements in a phased array. A set of " + str(
                maxlocroisize) + "x" + str(maxlocroisize) + "x1 voxels are positioned " + str(maxlocradfrac) +
            " of the distance from the phantom center to the edge along the direction of maximum sensitivity for the coil element. The average voxel values are plotted a function of coil element. The timecourses from each location are then crosscorrelated to assess common mode noise. Non-zero off-diagonal elements indicate correlation between channels, either due to geometric overlap or common mode system noise.\n")

        snrroicellstr = tableentrytag(
            paratag(boldtag("SNR at region of maximum sensitivity for each phased array element")) +
            imagetag("phasedarrayroisnrplot.png", myimwidth))
        sfnrroicellstr = tableentrytag(
            paratag(boldtag("SFNR at region of maximum sensitivity for each phased array element")) +
            imagetag("phasedarrayroisfnrplot.png", myimwidth))
        pproicellstr = tableentrytag(
            paratag(boldtag("p-p% variation at region of maximum sensitivity for each phased array element")) +
            imagetag("phasedarrayroippplot.png", myimwidth))
        ppdtroicellstr = tableentrytag(
            paratag(
                boldtag(
                    "p-p% variation after detrending at region of maximum sensitivity for each phased array element")) +
            imagetag("phasedarrayroippdtplot.png", myimwidth))
        crosscorrcellstr = tableentrytag(
            paratag(boldtag("Cross correlation of time data at each coil element's region of max sensitivity")) +
            imagetag("coilccmatrix.png", myimwidth))
        row1str = tablerowtag(snrroicellstr + sfnrroicellstr)
        row2str = tablerowtag(pproicellstr + ppdtroicellstr)
        row3str = tablerowtag(crosscorrcellstr)
        outfp.writelines(
            tablepropstag(phasedarrayroihdrstr + phasedarrayroiparastr + row1str + row2str + row3str, 500, "left"))

    ########################################################
    #
    # Max sensitivity ROI output
    #
    if isindividualcoil:
        maxlocroihdrstr = hruletag() + headertag("Maximum sensitivity region ROI Analysis")
        maxlocroiparastr = paratag(
            "This is an analysis of the temporal fluctuation in the ROI of maximum sensitivity for the individual coil element. A 3 x 3 x 1 voxel is automatically positioned at the center of gravity, and the average voxel value is plotted a function of time. Linear and quadratic terms are fit to the drift. The voxel statistics are reported with and without the drift removed.\n")

        rawroicellstr = tableentrytag(
            paratag(boldtag("Raw max sensitivity ROI plot")) +
            paratag(maxloctc_summary) +
            imagetag("maxlocroiplot.png", myimwidth))
        snrroicellstr = tableentrytag(
            paratag(boldtag("Central ROI over time")) +
            paratag(maxlocsnr_summary) +
            imagetag("maxlocroisnrplot.png", myimwidth))
        detrendedroicellstr = tableentrytag(
            paratag(boldtag("Detrended max sensitivity ROI plot")) +
            paratag(maxloctc_dt_summary) +
            imagetag("maxlocroidtplot.png", myimwidth))
        row1str = tablerowtag(rawroicellstr + detrendedroicellstr)
        row2str = tablerowtag(snrroicellstr)
        outfp.writelines(tablepropstag(maxlocroihdrstr + maxlocroiparastr + row1str + row2str, 500, "left"))

    ########################################################
    #
    # Ghost ROI output
    #
    ghostroihdrstr = hruletag() + headertag("Ghost ROI Analysis")
    ghostroiparastr = paratag(
        "This is an analysis of the amplitude and time variation of image ghosts.  Odd and even ghosts are assessed by calculating the ratio of the average signal in a " + str(
            ghostroisize) + "x" + str(ghostroisize) +
        "x1 ghost roi to the average amplitude in the center of the phantom. The ghost roi is placed at the edge of the field of view in the phase encode direction, outside the phantom, and in the center (even ghost) or at the edge (odd ghost) of the phantom position in the readout direction.\n")
    oddghostroicellstr = tableentrytag(
        paratag(boldtag("Odd ghost ROI plot")) +
        paratag(oddghosttc_summary) +
        imagetag("oddghostroiplot.png", myimwidth))
    evenghostroicellstr = tableentrytag(
        paratag(boldtag("Even ghost ROI plot")) +
        paratag(evenghosttc_summary) +
        imagetag("evenghostroiplot.png", myimwidth))
    row1str = tablerowtag(oddghostroicellstr + evenghostroicellstr)
    outfp.writelines(tablepropstag(ghostroihdrstr + ghostroiparastr + row1str, 500, "left"))

    ########################################################
    #
    # Weisskoff output
    #
    weisshdrstr = hruletag() + headertag("Weisskoff analysis")
    weissimagestr = imagetag("weisskoffplot.png", myimwidth)
    weisstablestr = sf.weisstable(roiareas, weisscvs, projcvs)
    row1str = tablerowtag(tableentrytag(weissimagestr) + tableentrytag(weisstablestr))
    row2str = tablerowtag(tableentrytag(bigtag(boldtag("RDC=" + str(weissrdc)))))
    outfp.writelines(tablepropstag(weisshdrstr + row1str + row2str, 500, "left"))

    datadict = {'Coil': info['Coil'],
                'Date': formatteddate,
                'Time': formattedtime,
                'DateTime': datetime,
                'Object': objectname,
                'Protocol': protocolname,
                'Element': info['ElementName'],
                'processed_as_individual': isindividualcoil,
                'object_radius_mm': object_radius_mm,
                'object_shape': object_shape,
                'center_of_mass_x': xcenterf,
                'center_of_mass_y': ycenterf,
                'center_of_mass_z': zcenterf}

    if not isindividualcoil:
        datadict['central_roi_raw_mean'] = centmean
        datadict['central_roi_raw_std'] = centstddev
        datadict['central_roi_raw_std%'] = 100.0 * centstddev / centmean
        datadict['central_roi_raw_min'] = centmin
        datadict['central_roi_raw_max'] = centmax
        datadict['central_roi_raw_p-p'] = centpp
        datadict['central_roi_raw_p-p%'] = 100.0 * centpp / centmean
        datadict['central_roi_detrended_mean'] = centmean_dt
        datadict['central_roi_detrended_std'] = centstddev_dt
        datadict['central_roi_detrended_std%'] = 100.0 * centstddev_dt / centmean_dt
        datadict['central_roi_detrended_min'] = centmin_dt
        datadict['central_roi_detrended_max'] = centmax_dt
        datadict['central_roi_detrended_p-p'] = centpp_dt
        datadict['central_roi_detrended_p-p%'] = 100.0 * centpp_dt / centmean_dt
        datadict['central_roi_SNR'] = centsnr
        datadict['central_roi_SFNR'] = centsfnr
        datadict['central_roi_polyfit_lin'] = 100.0 * centfitcoffs[1] / centfitcoffs[2]
        datadict['central_roi_polyfit_quad'] = 100.0 * centfitcoffs[0] / centfitcoffs[2]
        datadict['peripheral_roi_raw_mean'] = periphmean
        datadict['peripheral_roi_raw_std'] = periphstddev
        datadict['peripheral_roi_raw_std%'] = 100.0 * periphstddev / periphmean
        datadict['peripheral_roi_raw_min'] = periphmin
        datadict['peripheral_roi_raw_max'] = periphmax
        datadict['peripheral_roi_raw_p-p'] = periphpp
        datadict['peripheral_roi_raw_p-p%'] = 100.0 * periphpp / periphmean
        datadict['peripheral_roi_detrended_mean'] = periphmean_dt
        datadict['peripheral_roi_detrended_std'] = periphstddev_dt
        datadict['peripheral_roi_detrended_std%'] = 100.0 * periphstddev_dt / periphmean_dt
        datadict['peripheral_roi_detrended_min'] = periphmin_dt
        datadict['peripheral_roi_detrended_max'] = periphmax_dt
        datadict['peripheral_roi_detrended_p-p'] = periphpp_dt
        datadict['peripheral_roi_detrended_p-p%'] = 100.0 * periphpp_dt / periphmean_dt
        datadict['peripheral_roi_SNR'] = meanangperiphsnr
        datadict['peripheral_roi_SFNR'] = meanangperiphsfnr
        datadict['peripheral_roi_polyfit_lin'] = 100.0 * periphfitcoffs[1] / periphfitcoffs[2]
        datadict['peripheral_roi_polyfit_quad'] = 100.0 * periphfitcoffs[0] / periphfitcoffs[2]
    else:
        datadict['maxloc_roi_x'] = elementmaxpos[0]
        datadict['maxloc_roi_y'] = elementmaxpos[1]
        datadict['maxloc_roi_z'] = elementmaxpos[2]
        datadict['maxloc_roi_dirvec_x'] = elementdirvec[0] / elementdirnormfac
        datadict['maxloc_roi_dirvec_y'] = elementdirvec[1] / elementdirnormfac
        datadict['maxloc_roi_dirvec_z'] = elementdirvec[2] / elementdirnormfac
        datadict['maxloc_roi_mean'] = maxlocmean_dt
        datadict['maxloc_roi_std'] = maxlocstddev_dt
        datadict['maxloc_roi_std%'] = 100.0 * maxlocstddev_dt / maxlocmean_dt
        datadict['maxloc_roi_min'] = maxlocmin_dt
        datadict['maxloc_roi_max'] = maxlocmax_dt
        datadict['maxloc_roi_p-p'] = maxlocpp_dt
        datadict['maxloc_roi_p-p%'] = 100.0 * maxlocpp_dt / maxlocmean_dt
        datadict['maxloc_roi_SNR'] = maxlocsnr
        datadict['maxloc_roi_SFNR'] = maxlocsfnr
        datadict['maxloc_roi_polyfit_lin'] = 100.0 * maxlocfitcoffs[1] / maxlocfitcoffs[2]
        datadict['maxloc_roi_polyfit_quad'] = 100.0 * maxlocfitcoffs[0] / maxlocfitcoffs[2]
    datadict['odd_ghost_mean'] = oddghostmean
    datadict['odd_ghost_std'] = oddghoststddev
    datadict['odd_ghost_min'] = oddghostmin
    datadict['odd_ghost_max'] = oddghostmax
    datadict['odd_ghost_p-p'] = oddghostpp
    datadict['odd_ghost_p-p%'] = 100.0 * oddghostpp / oddghostmean
    datadict['even_ghost_mean'] = evenghostmean
    datadict['even_ghost_std'] = evenghoststddev
    datadict['even_ghost_min'] = evenghostmin
    datadict['even_ghost_max'] = evenghostmax
    datadict['even_ghost_p-p'] = evenghostpp
    datadict['even_ghost_p-p%'] = 100.0 * evenghostpp / evenghostmean
    datadict['weissrdc'] = weissrdc
    datadict['central_roi_drift%'] = centdrift
    datadict['peripheral_roi_drift%'] = periphdrift

    afp = open(pjoin(dirname, "procresults/dataquality.txt"), "w")
    for k in datadict:
        try:
            entrydesc = sf.formatlimits(limits[k])
            entryval = datadict[k]
            entryquality = sf.limitcheck(entryval, limits[k])
            afp.writelines(','.join((entrydesc, str(entryval), {0: 'Pass', 1: 'Warn', 2: 'Fail'}[entryquality])) + "\n")
        except KeyError:
            pass

    ########################################################
    #
    # Write summary text file
    summaryfile = pjoin(dirname, "procresults/analysissummary.txt")

    sumfp = open(summaryfile, "w")
    sumfp.writelines("Filename	" + summaryfile + "\n")
    sumfp.writelines("Coil	" + datadict['Coil'] + "\n")
    sumfp.writelines("Date	" + datadict['Date'] + "\n")
    sumfp.writelines("Time	" + datadict['Time'] + "\n")
    sumfp.writelines("DateTime	" + datadict['DateTime'] + "\n")
    sumfp.writelines("Object	" + datadict['Object'] + "\n")
    sumfp.writelines("Protocol	" + datadict['Protocol'] + "\n")
    sumfp.writelines("Element	" + datadict['Element'] + "\n")
    sumfp.writelines("processed_as_individual	" + str(datadict['processed_as_individual']) + "\n")
    sumfp.writelines("object_radius_mm	" + str(datadict['object_radius_mm']) + "\n")
    sumfp.writelines("object_shape	" + str(datadict['object_shape']) + "\n")
    sumfp.writelines("center_of_mass_x	" + str(datadict['center_of_mass_x']) + "\n")
    sumfp.writelines("center_of_mass_y	" + str(datadict['center_of_mass_y']) + "\n")
    sumfp.writelines("center_of_mass_z	" + str(datadict['center_of_mass_z']) + "\n")
    if not isindividualcoil:
        sumfp.writelines("central_roi_raw_mean	" + str(datadict['central_roi_raw_mean']) + "\n")
        sumfp.writelines("central_roi_raw_std	" + str(datadict['central_roi_raw_std']) + "\n")
        sumfp.writelines("central_roi_raw_std%	" + str(datadict['central_roi_raw_std%']) + "\n")
        sumfp.writelines("central_roi_raw_min	" + str(datadict['central_roi_raw_min']) + "\n")
        sumfp.writelines("central_roi_raw_max	" + str(datadict['central_roi_raw_max']) + "\n")
        sumfp.writelines("central_roi_raw_p-p	" + str(datadict['central_roi_raw_p-p']) + "\n")
        sumfp.writelines("central_roi_raw_p-p%	" + str(datadict['central_roi_raw_p-p%']) + "\n")
        sumfp.writelines("central_roi_detrended_mean	" + str(datadict['central_roi_detrended_mean']) + "\n")
        sumfp.writelines("central_roi_detrended_std	" + str(datadict['central_roi_detrended_std']) + "\n")
        sumfp.writelines("central_roi_detrended_std%	" + str(datadict['central_roi_detrended_std%']) + "\n")
        sumfp.writelines("central_roi_detrended_min	" + str(datadict['central_roi_detrended_min']) + "\n")
        sumfp.writelines("central_roi_detrended_max	" + str(datadict['central_roi_detrended_max']) + "\n")
        sumfp.writelines("central_roi_detrended_p-p	" + str(datadict['central_roi_detrended_p-p']) + "\n")
        sumfp.writelines("central_roi_detrended_p-p%	" + str(datadict['central_roi_detrended_p-p%']) + "\n")
        sumfp.writelines("central_roi_SNR	" + str(datadict['central_roi_SNR']) + "\n")
        sumfp.writelines("central_roi_SFNR	" + str(datadict['central_roi_SFNR']) + "\n")
        sumfp.writelines("central_roi_polyfit_lin	" + str(datadict['central_roi_polyfit_lin']) + "\n")
        sumfp.writelines("central_roi_polyfit_quad	" + str(datadict['central_roi_polyfit_quad']) + "\n")
        sumfp.writelines("peripheral_roi_raw_mean	" + str(datadict['peripheral_roi_raw_mean']) + "\n")
        sumfp.writelines("peripheral_roi_raw_std	" + str(datadict['peripheral_roi_raw_std']) + "\n")
        sumfp.writelines("peripheral_roi_raw_std%	" + str(datadict['peripheral_roi_raw_std%']) + "\n")
        sumfp.writelines("peripheral_roi_raw_min	" + str(datadict['peripheral_roi_raw_min']) + "\n")
        sumfp.writelines("peripheral_roi_raw_max	" + str(datadict['peripheral_roi_raw_max']) + "\n")
        sumfp.writelines("peripheral_roi_raw_p-p	" + str(datadict['peripheral_roi_raw_p-p']) + "\n")
        sumfp.writelines("peripheral_roi_raw_p-p%	" + str(datadict['peripheral_roi_raw_p-p%']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_mean	" + str(datadict['peripheral_roi_detrended_mean']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_std	" + str(datadict['peripheral_roi_detrended_std']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_std%	" + str(datadict['peripheral_roi_detrended_std%']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_min	" + str(datadict['peripheral_roi_detrended_min']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_max	" + str(datadict['peripheral_roi_detrended_max']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_p-p	" + str(datadict['peripheral_roi_detrended_p-p']) + "\n")
        sumfp.writelines("peripheral_roi_detrended_p-p%	" + str(datadict['peripheral_roi_detrended_p-p%']) + "\n")
        sumfp.writelines("peripheral_roi_SNR	" + str(datadict['peripheral_roi_SNR']) + "\n")
        sumfp.writelines("peripheral_roi_SFNR	" + str(datadict['peripheral_roi_SFNR']) + "\n")
        sumfp.writelines("peripheral_roi_polyfit_lin	" + str(datadict['peripheral_roi_polyfit_lin']) + "\n")
        sumfp.writelines("peripheral_roi_polyfit_quad	" + str(datadict['peripheral_roi_polyfit_quad']) + "\n")
    else:
        sumfp.writelines("maxloc_roi_x	" + str(datadict['maxloc_roi_x']) + "\n")
        sumfp.writelines("maxloc_roi_y	" + str(datadict['maxloc_roi_y']) + "\n")
        sumfp.writelines("maxloc_roi_z	" + str(datadict['maxloc_roi_z']) + "\n")
        sumfp.writelines("maxloc_roi_dirvec_x	" + str(datadict['maxloc_roi_dirvec_x']) + "\n")
        sumfp.writelines("maxloc_roi_dirvec_y	" + str(datadict['maxloc_roi_dirvec_y']) + "\n")
        sumfp.writelines("maxloc_roi_dirvec_z	" + str(datadict['maxloc_roi_dirvec_z']) + "\n")
        sumfp.writelines("maxloc_roi_mean	" + str(datadict['maxloc_roi_mean']) + "\n")
        sumfp.writelines("maxloc_roi_std	" + str(datadict['maxloc_roi_std']) + "\n")
        sumfp.writelines("maxloc_roi_std%	" + str(datadict['maxloc_roi_std%']) + "\n")
        sumfp.writelines("maxloc_roi_min	" + str(datadict['maxloc_roi_min']) + "\n")
        sumfp.writelines("maxloc_roi_max	" + str(datadict['maxloc_roi_max']) + "\n")
        sumfp.writelines("maxloc_roi_p-p	" + str(datadict['maxloc_roi_p-p']) + "\n")
        sumfp.writelines("maxloc_roi_p-p%	" + str(datadict['maxloc_roi_p-p%']) + "\n")
        sumfp.writelines("maxloc_roi_SNR	" + str(datadict['maxloc_roi_SNR']) + "\n")
        sumfp.writelines("maxloc_roi_SFNR	" + str(datadict['maxloc_roi_SFNR']) + "\n")
        sumfp.writelines("maxloc_roi_polyfit_lin	" + str(datadict['maxloc_roi_polyfit_lin']) + "\n")
        sumfp.writelines("maxloc_roi_polyfit_quad	" + str(datadict['maxloc_roi_polyfit_quad']) + "\n")
    sumfp.writelines("odd_ghost_mean	" + str(datadict['odd_ghost_mean']) + "\n")
    sumfp.writelines("odd_ghost_std	" + str(datadict['odd_ghost_std']) + "\n")
    sumfp.writelines("odd_ghost_min	" + str(datadict['odd_ghost_min']) + "\n")
    sumfp.writelines("odd_ghost_max	" + str(datadict['odd_ghost_max']) + "\n")
    sumfp.writelines("odd_ghost_p-p	" + str(datadict['odd_ghost_p-p']) + "\n")
    sumfp.writelines("odd_ghost_p-p%	" + str(datadict['odd_ghost_p-p%']) + "\n")
    sumfp.writelines("even_ghost_mean	" + str(datadict['even_ghost_mean']) + "\n")
    sumfp.writelines("even_ghost_std	" + str(datadict['even_ghost_std']) + "\n")
    sumfp.writelines("even_ghost_min	" + str(datadict['even_ghost_min']) + "\n")
    sumfp.writelines("even_ghost_max	" + str(datadict['even_ghost_max']) + "\n")
    sumfp.writelines("even_ghost_p-p	" + str(datadict['even_ghost_p-p']) + "\n")
    sumfp.writelines("even_ghost_p-p%	" + str(datadict['even_ghost_p-p%']) + "\n")
    sumfp.writelines("weissrdc	" + str(datadict['weissrdc']) + "\n")
    sumfp.writelines("central_roi_drift%	" + str(datadict['central_roi_drift%']) + "\n")
    sumfp.writelines("peripheral_roi_drift%	" + str(datadict['peripheral_roi_drift%']) + "\n")
    outfp.writelines("</body>\n")

    tpl = makolookup.get_template('stability.html')
    with open(pjoin(dirname, 'procresults/mako.html'), 'w') as fp:
        fp.write(tpl.render(**locals()))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Calculate stability values and create the output web page.')
    parser.add_argument('dirname', help='the directory where the 4D NIFTI file is located')
    parser.add_argument('filename', help='the name of the 4D NIFTI file')
    parser.add_argument('starttime', help='the number of tr periods to skip at the beginning of the file')

    com = parser.add_argument_group('use this location as the center of mass of the phantom')
    com.add_argument('initxcenter', nargs='?', metavar='xcenter')
    com.add_argument('initycenter', nargs='?', metavar='ycenter')
    com.add_argument('initzcenter', nargs='?', metavar='zcenter')

    args = parser.parse_args()
    if None in (args.initxcenter, args.initycenter, args.initzcenter) and not args.initxcenter == args.initycenter == args.initzcenter:
        parser.error('If you set one center of mass parameter, you must set all three.')

    stabilitycalc(args.dirname, args.filename, int(args.starttime), args.initxcenter, args.initycenter, args.initzcenter)
