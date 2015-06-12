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

from mako.lookup import TemplateLookup

makolookup = TemplateLookup(directories=['./tpl'])

import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
import stabilityfuncs as sf
import studyinfo


def stabilitycalc(dirname, filename, starttime, initxcenter=None, initycenter=None, initzcenter=None):
    """create the stability report for a scan"""

    logging.debug('dirname: {}\n'
                  'filename: {}\n'
                  'starttime: {}\n'
                  'centers: {}, {}, {}'.format(dirname, filename, starttime, initxcenter, initycenter, initzcenter))

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

    thisdate = time.strftime("%m/%d/%Y %H:%M:%S", time.localtime())

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
    objectname = "Unknown"

    object_radius_mm = np.sqrt(objectradiusx_mm * objectradiusy_mm)
    object_shape = objectradiusy / objectradiusx

    birn_phantom_radiuscheck = sf.limitcheck(object_radius_mm, limits['BIRNphantom_rad'])
    birn_phantom_shapecheck = sf.limitcheck(object_shape, limits['BIRNphantom_shape'])
    if (birn_phantom_radiuscheck < 2) and (birn_phantom_shapecheck < 2):
        objectname = "BIRN phantom"
        logging.debug("setting objectname to BIRN phantom")

    head_radiuscheck = sf.limitcheck(object_radius_mm, limits['head_rad'])
    head_shapecheck = sf.limitcheck(object_shape, limits['head_shape'])
    if (head_radiuscheck < 2) or (head_shapecheck < 2):
        objectname = "Head"
        logging.debug("setting objectname to Head")

    is_birn_sequence = True
    is_birn_protocol = False

    if (xsize != 64) or (ysize != 64) or (numslices != 28) or (tr != 2.0):
        is_birn_sequence = False
    if is_birn_sequence and (objectname == 'BIRN phantom'):
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

    def qualitypercent(n, basis, lim):
        percent = n / basis * 100.0
        quality = sf.limitcheck(percent, limits[lim])
        return qualitytag('(%4.4f%%)', quality) % percent

    centstddev_qualitytag = qualitypercent(centstddev, centmean, 'central_roi_raw_std%')
    centpp_qualitytag = qualitypercent(centpp, centmean, 'central_roi_raw_p-p%')

    centmean_dt = np.mean(detrendedcenttc)
    centstddev_dt = np.std(detrendedcenttc)
    centmin_dt = np.min(detrendedcenttc)
    centmax_dt = np.max(detrendedcenttc)
    centpp_dt = np.ptp(detrendedcenttc)

    centstddev_dt_qualitytag = qualitypercent(centstddev_dt, centmean_dt, 'central_roi_detrended_std%')
    centpp_dt_qualitytag = qualitypercent(centpp_dt, centmean_dt, 'central_roi_detrended_p-p%')

    #############################
    #
    #       Maximum value ROI analysis
    #
    logging.debug("Finding and analyzing maximum signal ROI...")

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

        maxloc_qualitytag = qualitypercent(maxlocstddev, maxlocmean, 'maxlocroi_rawstddev')
        maxloc_pp_qualitytag = qualitypercent(maxlocpp, maxlocmean, 'maxlocroi_rawpp')

        maxlocmean_dt = np.mean(detrendedmaxloctc)
        maxlocstddev_dt = np.std(detrendedmaxloctc)
        maxlocmin_dt = np.min(detrendedmaxloctc)
        maxlocmax_dt = np.max(detrendedmaxloctc)
        maxlocpp_dt = np.ptp(detrendedmaxloctc)

        maxloc_dt_qualitytag = qualitypercent(maxlocstddev_dt, maxlocmean_dt, 'maxlocroi_dtstddev')
        maxlocpp_dt_qualitytag = qualitypercent(maxlocpp_dt, maxlocmean_dt, 'maxlocroi_dtpp')

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
            if selecteddata.shape[1] == 1 and coildata[ele]['zloc'] != 0:
                # single slice
                coildata[ele]['zloc'] = 0
                logging.debug('single slice data, changing zloc from {} to 0'.format(coildata[ele]['zloc']))
            roi = sf.setroilims(round(coildata[ele]['xloc']), round(coildata[ele]['yloc']), phasedarraysize)
            if not isindividualcoil:
                sf.markroi(roi, round(coildata[ele]['zloc']), roislice, 0.95 * rawmeanmax)
            timecourse = sf.getroimeantc(selecteddata, roi, zcenter)
            phasedarraytcs[i, :] = timecourse[:]
            snrvec = sf.getroisnr(selecteddata, roi, round(coildata[ele]['zloc']))
            phasedarrayroisnrs[i] = np.mean(snrvec)
            phasedarrayroisfnrs[i] = sf.getroival(sfnrslice, roi, round(coildata[ele]['zloc']))
            phasedarrayroimeans[i] = np.mean(timecourse)
            phasedarrayroipps[i] = np.ptp(timecourse)
            phasedarrayfitcoffs = np.polyfit(timepoints, timecourse, 2)
            phasedarrayfittcs[i, :] = sf.trendgen(timepoints, phasedarrayfitcoffs)
            phasedarraydttcs[i, :] = phasedarraytcs[i, :] - phasedarrayfittcs[i, :]
            phasedarrayroimeans_dt[i] = np.mean(phasedarraydttcs[i, :])
            phasedarraydttcs_demeaned[i, :] = phasedarraydttcs[i, :] - phasedarrayroimeans_dt[i]
            phasedarrayroipps_dt[i] = np.ptp(phasedarraydttcs[i, :])

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
    periph_qualitytag = qualitypercent(periphstddev, periphmean, 'peripheral_roi_raw_std%')
    periphpp_qualitytag = qualitypercent(periphpp, periphmean, 'peripheral_roi_raw_p-p%')

    periphmean_dt = np.mean(detrendedperiphtc)
    periphstddev_dt = np.std(detrendedperiphtc)
    periphmin_dt = np.min(detrendedperiphtc)
    periphmax_dt = np.max(detrendedperiphtc)
    periphpp_dt = np.ptp(detrendedperiphtc)
    periph_dt_qualitytag = qualitypercent(periphstddev_dt, periphmean_dt, 'peripheral_roi_detrended_std%')
    periphpp_dt_qualitytag = qualitypercent(periphpp_dt, periphmean_dt, 'peripheral_roi_detrended_p-p%')

    # do calculations regarding angular dependance
    meanangperiphval = np.mean(periphangmeans)
    meanangperiphsfnr = np.mean(periphangsfnrs)
    meanangperiphsnr = np.mean(avgperiphsnrvec)

    ptpangperiphval = np.ptp(periphangmeans)
    periphangintensity_qualitytag = qualitypercent(ptpangperiphval, meanangperiphval, 'peripheral_angle_p-p%')

    periphang_sfnr_ptp = np.ptp(periphangsfnrs)
    periphang_sfnr_pp_qualitytag = qualitypercent(periphang_sfnr_ptp, meanangperiphsfnr, 'peripheral_angle_SFNR_p-p%')

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

    #############################
    #
    # Output
    #

    # statistical images section

    # noinspection PyShadowingNames
    def slicepic(inputslice, caption, minp, maxp, dirname, outputname, colormap):
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
        slicepic(np.nan_to_num(objectmask * eodiffpcimage), "Even odd diff percent image", eodiffpcmin, eodiffpcmax,
                 dirname, 'eodiffpcimage', 0)
        slicepic(ppslice, "Peak to peak image", ppmin, ppmax, dirname, 'ppimage', 0)

    def makefig(figk):
        """make a figure using the figure dictionary"""
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

    # data quality report

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
        datadict.update({'central_roi_raw_mean': centmean,
                         'central_roi_raw_std': centstddev,
                         'central_roi_raw_std%': 100.0 * centstddev / centmean,
                         'central_roi_raw_min': centmin,
                         'central_roi_raw_max': centmax,
                         'central_roi_raw_p-p': centpp,
                         'central_roi_raw_p-p%': 100.0 * centpp / centmean,
                         'central_roi_detrended_mean': centmean_dt,
                         'central_roi_detrended_std': centstddev_dt,
                         'central_roi_detrended_std%': 100.0 * centstddev_dt / centmean_dt,
                         'central_roi_detrended_min': centmin_dt,
                         'central_roi_detrended_max': centmax_dt,
                         'central_roi_detrended_p-p': centpp_dt,
                         'central_roi_detrended_p-p%': 100.0 * centpp_dt / centmean_dt,
                         'central_roi_SNR': centsnr,
                         'central_roi_SFNR': centsfnr,
                         'central_roi_polyfit_lin': 100.0 * centfitcoffs[1] / centfitcoffs[2],
                         'central_roi_polyfit_quad': 100.0 * centfitcoffs[0] / centfitcoffs[2],
                         'peripheral_roi_raw_mean': periphmean,
                         'peripheral_roi_raw_std': periphstddev,
                         'peripheral_roi_raw_std%': 100.0 * periphstddev / periphmean,
                         'peripheral_roi_raw_min': periphmin,
                         'peripheral_roi_raw_max': periphmax,
                         'peripheral_roi_raw_p-p': periphpp,
                         'peripheral_roi_raw_p-p%': 100.0 * periphpp / periphmean,
                         'peripheral_roi_detrended_mean': periphmean_dt,
                         'peripheral_roi_detrended_std': periphstddev_dt,
                         'peripheral_roi_detrended_std%': 100.0 * periphstddev_dt / periphmean_dt,
                         'peripheral_roi_detrended_min': periphmin_dt,
                         'peripheral_roi_detrended_max': periphmax_dt,
                         'peripheral_roi_detrended_p-p': periphpp_dt,
                         'peripheral_roi_detrended_p-p%': 100.0 * periphpp_dt / periphmean_dt,
                         'peripheral_roi_SNR': meanangperiphsnr,
                         'peripheral_roi_SFNR': meanangperiphsfnr,
                         'peripheral_roi_polyfit_lin': 100.0 * periphfitcoffs[1] / periphfitcoffs[2],
                         'peripheral_roi_polyfit_quad': 100.0 * periphfitcoffs[0] / periphfitcoffs[2]})
    else:
        datadict.update({'maxloc_roi_x': elementmaxpos[0],
                         'maxloc_roi_y': elementmaxpos[1],
                         'maxloc_roi_z': elementmaxpos[2],
                         'maxloc_roi_dirvec_x': elementdirvec[0] / elementdirnormfac,
                         'maxloc_roi_dirvec_y': elementdirvec[1] / elementdirnormfac,
                         'maxloc_roi_dirvec_z': elementdirvec[2] / elementdirnormfac,
                         'maxloc_roi_mean': maxlocmean_dt,
                         'maxloc_roi_std': maxlocstddev_dt,
                         'maxloc_roi_std%': 100.0 * maxlocstddev_dt / maxlocmean_dt,
                         'maxloc_roi_min': maxlocmin_dt,
                         'maxloc_roi_max': maxlocmax_dt,
                         'maxloc_roi_p-p': maxlocpp_dt,
                         'maxloc_roi_p-p%': 100.0 * maxlocpp_dt / maxlocmean_dt,
                         'maxloc_roi_SNR': maxlocsnr,
                         'maxloc_roi_SFNR': maxlocsfnr,
                         'maxloc_roi_polyfit_lin': 100.0 * maxlocfitcoffs[1] / maxlocfitcoffs[2],
                         'maxloc_roi_polyfit_quad': 100.0 * maxlocfitcoffs[0] / maxlocfitcoffs[2]})
    datadict.update({'odd_ghost_mean': oddghostmean,
                     'odd_ghost_std': oddghoststddev,
                     'odd_ghost_min': oddghostmin,
                     'odd_ghost_max': oddghostmax,
                     'odd_ghost_p-p': oddghostpp,
                     'odd_ghost_p-p%': 100.0 * oddghostpp / oddghostmean,
                     'even_ghost_mean': evenghostmean,
                     'even_ghost_std': evenghoststddev,
                     'even_ghost_min': evenghostmin,
                     'even_ghost_max': evenghostmax,
                     'even_ghost_p-p': evenghostpp,
                     'even_ghost_p-p%': 100.0 * evenghostpp / evenghostmean,
                     'weissrdc': weissrdc,
                     'central_roi_drift%': centdrift,
                     'peripheral_roi_drift%': periphdrift})

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
    # Write summary text file and report

    tpl = makolookup.get_template('analysissummary.txt')
    summaryfile = pjoin(dirname, "procresults/analysissummary.txt")
    with open(summaryfile, 'w') as fp:
        fp.write(tpl.render(**locals()))

    tpl = makolookup.get_template('stability.html')
    with open(pjoin(dirname, 'procresults/output.html'), 'w') as fp:
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
    if None in (args.initxcenter, args.initycenter,
                args.initzcenter) and not args.initxcenter == args.initycenter == args.initzcenter:
        parser.error('If you set one center of mass parameter, you must set all three.')

    stabilitycalc(args.dirname, args.filename, int(args.starttime), args.initxcenter, args.initycenter,
                  args.initzcenter)
