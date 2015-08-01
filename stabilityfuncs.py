#!/usr/bin/env python

"""stabilityfuncs.py: provides general helper functions, especially for stabilitycalc"""

import csv
from collections import OrderedDict
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
from collections import namedtuple
import ConfigParser


Stats = namedtuple('Stats', 'mean stddev var max min ptp')


def dict_from_tsvfile(filename):
    """open a tab separated two column file, return it as a str->str dict"""

    d = {}
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            d[row[0]] = row[1]

    return d


def stabilityparms(option, section='paths'):
    """
    retrieve a value from the configuration file.
    for legacy use, assume we're looking for a path.
    """

    config = ConfigParser.ConfigParser()

    try:
        config.read('config/stability.ini')
    except IOError:
        logging.critical('Failed to open configuration file.')
        exit(1)

    try:
        return config.get(section, option.lower())
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        return ''


def getlimits(coil):
    # noinspection PyShadowingNames
    def get_lim_csv(coil):
        d = {}
        try:
            with open('config/{}_limits.csv'.format(coil)) as csvfile:
                reader = csv.DictReader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
                for r in reader:
                    d[r['var']] = r
        except IOError:
            logging.critical("getlimits: coil not recognized!")
            exit(1)
        return d

    limitdict = get_lim_csv('coil_independent')
    limitdict.update(get_lim_csv('sample'))
    limitdict.update(get_lim_csv(coil))

    return limitdict


def getphasedarraydata(coil):
    d = OrderedDict()
    try:
        with open('config/{}_coildata.csv'.format(coil)) as csvfile:
            reader = csv.DictReader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
            for r in reader:
                d[r['element']] = r
        return d
    except IOError:
        logging.debug("getphasedarraydata: Not a phased array.")
        return False


def makemask(inputim, inputthresh, useabs):
    if useabs < 1:
        thethreshval = getfracval(inputim, inputthresh)
    else:
        thethreshval = inputthresh
    themask = sp.where(inputim > thethreshval, 1.0, 0.0)
    return themask


def vecnorm(thevec):
    return np.sqrt(np.square(thevec).sum())


def formatlimits(lim):
    return '"{}",{},{},{},{}'.format(lim['description'],
                                     lim['warn_min'],
                                     lim['good_min'],
                                     lim['good_max'],
                                     lim['warn_max'])


def limitcheck(n, lim):
    # check to see if a parameter falls within preset limits.

    # check for legacy input mode and convert if needed
    if type(lim) is tuple:
        lim = {'good_min': lim[0][0], 'good_max': lim[0][1], 'warn_min': lim[1][0], 'warn_max': lim[1][1]}

    retval = 2  # start with the assumption that the data is bad
    if (float(n) >= float(lim['warn_min'])) and (float(n) <= float(lim['warn_max'])):
        retval = 1  # number falls within the warning limits
    if (float(n) >= float(lim['good_min'])) and (float(n) <= float(lim['good_max'])):
        retval = 0  # number falls within the good limits
    return retval


def trendgen(thexvals, thefitcoffs):
    # generate the polynomial fit timecourse from the coefficients
    theshape = thefitcoffs.shape
    order = theshape[0] - 1
    thepoly = thexvals
    thefit = 0.0 * thexvals
    if order > 0:
        for i in range(1, order + 1):
            thefit = thefit + thefitcoffs[order - i] * thepoly
            thepoly = np.multiply(thepoly, thexvals)
    return thefit


def completerobust(thearray):
    # calculate the robust range of the all voxels
    themin = getfracval(thearray, 0.02)
    themax = getfracval(thearray, 0.98)
    return [themin, themax]


def nzrobust(thearray):
    # calculate the robust range of the non-zero voxels
    themin = getnzfracval(thearray, 0.02)
    themax = getnzfracval(thearray, 0.98)
    return [themin, themax]


def nzminmax(thearray):
    # calculate the min and max of the non-zero voxels
    flatarray = np.ravel(thearray)
    nzindices = np.nonzero(flatarray)
    theflatarray = flatarray[nzindices]
    themax = np.max(theflatarray)
    themin = np.min(theflatarray)
    return [themin, themax]


def nzstats(thearray):
    # calculate the stats of the non-zero voxels
    flatarray = np.ravel(thearray)
    nzindices = np.nonzero(np.ravel(thearray))
    return Stats(np.mean(flatarray[nzindices]),
                 np.std(flatarray[nzindices]),
                 np.var(flatarray[nzindices]),
                 np.max(flatarray[nzindices]),
                 np.min(flatarray[nzindices]),
                 np.ptp(flatarray[nzindices]))


def showtc(thexvals, theyvals, thelabel):
    # show an roi timecourse plot
    w, h = plt.figaspect(0.25)
    roiplot = plt.figure(figsize=(w, h))
    roisubplot = roiplot.add_subplot(111)
    roisubplot.plot(thexvals, theyvals, 'b')
    roisubplot.grid(True)
    # roisubplot.axes.Subplot.set_pad(0.1)
    for tick in roisubplot.xaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    for tick in roisubplot.yaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    roisubplot.set_title(thelabel, fontsize=30)


def showtc2(thexvals, theyvals, thefitvals, thelabel):
    # show an roi timecourse plot and a fit line
    w, h = plt.figaspect(0.25)
    roiplot = plt.figure(figsize=(w, h))
    roisubplot = roiplot.add_subplot(111)
    roisubplot.plot(thexvals, theyvals, 'b', thexvals, thefitvals, 'g')
    roisubplot.grid(True)
    # roisubplot.axes.Subplot.set_pad(0.1)
    for tick in roisubplot.xaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    for tick in roisubplot.yaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    roisubplot.set_title(thelabel, fontsize=30)


def showweisskoff(theareas, thestddevs, theprojstddevs, thelabel):
    # initialize and show a loglog Weiskoff plot
    logging.debug("Generating plot for {}".format(thelabel))
    w, h = plt.figaspect(1.0)
    roiplot = plt.figure(figsize=(w, h))
    roiplot.subplots_adjust(hspace=0.35)
    roisubplot = roiplot.add_subplot(111)
    thestddevs += 0.00000001
    roisubplot.loglog(theareas, thestddevs, 'r', theareas, theprojstddevs, 'k', basex=10)
    roisubplot.grid(True)


def showslice2(thedata, thelabel, minval, maxval, colormap):
    # initialize and show a 2D slice from a dataset in greyscale
    plt.figure(figsize=plt.figaspect(1.0))
    theshape = thedata.shape
    numslices = theshape[0]
    ysize = theshape[1]
    xsize = theshape[2]
    slicesqrt = int(np.ceil(np.sqrt(numslices)))
    theslice = np.zeros((ysize * slicesqrt, xsize * slicesqrt))
    for i in range(numslices):
        ypos = int(i / slicesqrt) * ysize
        xpos = int(i % slicesqrt) * xsize
        theslice[ypos:ypos + ysize, xpos:xpos + xsize] = thedata[i, :, :]
    if plt.isinteractive():
        plt.ioff()
    plt.axis('off')
    plt.axis('equal')
    plt.subplots_adjust(hspace=0.0)
    plt.axes([0, 0, 1, 1], frameon=False)
    if colormap == 0:
        thecmap = cm.gray
    else:
        mycmdata1 = {
            'red': ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 1., 1.)),
            'green': ((0., 0., 0.), (0.5, 1.0, 1.0), (1., 0., 0.)),
            'blue': ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 0., 0.))
        }
        thecmap = colors.LinearSegmentedColormap('mycm', mycmdata1)
    plt.imshow(theslice, vmin=minval, vmax=maxval, interpolation='nearest', label=thelabel, aspect='equal',
               cmap=thecmap)


def showslice3(thedata, thelabel, minval, maxval, colormap):
    # initialize and show a 2D slice from a dataset in greyscale
    theshape = thedata.shape
    ysize = theshape[0]
    xsize = theshape[1]
    np.zeros((ysize, xsize))
    if plt.isinteractive():
        plt.ioff()
    plt.axis('off')
    plt.axis('equal')
    plt.subplots_adjust(hspace=0.0)
    plt.axes([0, 0, 1, 1], frameon=False)
    if colormap == 0:
        thecmap = cm.gray
    else:
        mycmdata1 = {
            'red': ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 1., 1.)),
            'green': ((0., 0., 0.), (0.5, 1.0, 1.0), (1., 0., 0.)),
            'blue': ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 0., 0.))
        }
        thecmap = colors.LinearSegmentedColormap('mycm', mycmdata1)
    plt.imshow(thedata, vmin=minval, vmax=maxval, interpolation='nearest', label=thelabel, aspect='equal', cmap=thecmap)


def smooth(x, window_len=11, window='hanning'):
    # this routine comes from a scipy.org Cookbook
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also: 

    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window should be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[2 * x[0] - x[window_len:1:-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='same')
    return y[window_len - 1:-window_len + 1]


def findsepval(datamat):
    # Find the image intensity value that cleanly separates background from image
    numbins = 200
    themax = datamat.max()
    themin = datamat.min()
    (meanhist, bins) = np.histogram(datamat, bins=numbins, range=(themin, themax))
    smoothhist = smooth(meanhist)
    currentpos = int(numbins * 0.05)
    for i in range(currentpos + 1, numbins):
        if smoothhist[i] < smoothhist[currentpos]:
            currentpos = i
        if smoothhist[i] > 1.2 * smoothhist[currentpos]:
            break
    cummeanhist = np.cumsum(meanhist)
    cumfrac = (1.0 * cummeanhist[currentpos]) / (1.0 * cummeanhist[numbins - 1])
    sepval = bins[currentpos]
    return [sepval, cumfrac]


def getfracval(datamat, thefrac):
    # Find the image intensity value which thefrac of the non-zero voxels in the image exceed
    numbins = 200
    themax = datamat.max()
    themin = datamat.min()
    (meanhist, bins) = np.histogram(datamat, bins=numbins, range=(themin, themax))
    cummeanhist = np.cumsum(meanhist)
    target = cummeanhist[numbins - 1] * thefrac
    for i in range(numbins):
        if cummeanhist[i] >= target:
            return bins[i]
    return 0.0


def getnzfracval(datamat, thefrac):
    # Find the image intensity value which thefrac of the non-zero voxels in the image exceed
    numbins = 200
    (themin, themax) = nzminmax(datamat)
    (meanhist, bins) = np.histogram(datamat, bins=numbins, range=(themin, themax))
    cummeanhist = np.cumsum(meanhist)
    target = cummeanhist[numbins - 1] * thefrac
    for i in range(numbins):
        if cummeanhist[i] >= target:
            return bins[i]
    return 0.0


# noinspection PyPep8Naming,PyUnresolvedReferences
def findCOM(datamat):
    # find the center of mass of a 2D or 3D image
    Mx = 0.0
    My = 0.0
    Mz = 0.0
    mass = 0.0
    arrdims = np.shape(datamat)

    if datamat.ndim == 2:
        for i in range(arrdims[0]):
            for j in range(arrdims[1]):
                val = datamat[i, j]
                My += (i * val)
                Mx += (j * val)
                mass += val
        COM = (Mx / mass, My / mass, 0.0)
    if datamat.ndim == 3:
        for i in range(arrdims[0]):
            for j in range(arrdims[1]):
                for k in range(arrdims[2]):
                    val = datamat[i, j, k]
                    Mz += (i * val)
                    My += (j * val)
                    Mx += (k * val)
                    mass += val
        COM = (Mx / mass, My / mass, Mz / mass)
    return COM


def markroi(theinputroi, zpos, roislice, theval):
    # given an roi and a position, mark an roi
    xstart = theinputroi[0][0]
    xend = theinputroi[1][0]
    ystart = theinputroi[0][1]
    yend = theinputroi[1][1]
    roislice[zpos, ystart:yend, xstart:xend] = theval


def setroilims(xpos, ypos, size):
    # given a location and a size, define the corners of an roi
    if (size % 2) == 0:
        halfsize = size / 2
        return (((int(round(xpos - halfsize)), int(round(ypos - halfsize))),
                 (int(round(xpos + halfsize)), int(round(ypos + halfsize)))))
    else:
        halfsize = (size - 1) / 2
        return (((int(round(xpos - halfsize)), int(round(ypos - halfsize))),
                 (int(round(xpos + halfsize + 1)), int(round(ypos + halfsize + 1)))))


def getroisnr(theimage, theroi, zpos):
    # get an snr timecourse from the voxels of an roi
    xstart = theroi[0][0]
    xend = theroi[1][0]
    ystart = theroi[0][1]
    yend = theroi[1][1]
    thesubreg = theimage[:, zpos, ystart:yend, xstart:xend]
    theshape = thesubreg.shape
    numtimepoints = theshape[0]
    themeans = np.zeros(numtimepoints)
    thestddevs = np.zeros(numtimepoints)
    themax = np.zeros(numtimepoints)
    themin = np.zeros(numtimepoints)
    thesnrs = np.zeros(numtimepoints)
    timeindex = np.arange(numtimepoints)
    for i in timeindex:
        themeans[i] = np.mean(np.ravel(thesubreg[i, :, :]))
        thestddevs[i] = np.std(np.ravel(thesubreg[i, :, :]))
        themax[i] = np.max(np.ravel(thesubreg[i, :, :]))
        themin[i] = np.min(np.ravel(thesubreg[i, :, :]))
        thesnrs[i] = themeans[i] / thestddevs[i]
    return thesnrs


def getroivoxels(theimage, theroi, zpos):
    # get all the voxels from an roi and return a 2d (time by space) array
    xstart = theroi[0][0]
    xend = theroi[1][0]
    ystart = theroi[0][1]
    yend = theroi[1][1]
    thesubreg = theimage[:, zpos, ystart:yend, xstart:xend]
    theshape = thesubreg.shape
    numtimepoints = theshape[0]
    thevoxels = np.zeros((numtimepoints, theshape[1] * theshape[2]))
    timeindex = np.arange(numtimepoints)
    for i in timeindex:
        thevoxels[i, :] = np.ravel(thesubreg[i, :, :])
    return thevoxels


def getroistdtc(theimage, theroi, zpos):
    # get a standard deviation timecourse from the voxels of an roi
    xstart = theroi[0][0]
    xend = theroi[1][0]
    ystart = theroi[0][1]
    yend = theroi[1][1]
    thesubreg = theimage[:, zpos, ystart:yend, xstart:xend]
    theshape = thesubreg.shape
    numtimepoints = theshape[0]
    thestds = np.zeros(numtimepoints)
    timeindex = np.arange(numtimepoints)
    for i in timeindex:
        thestds[i] = np.std(np.ravel(thesubreg[i, :, :]))
    return thestds


def getroimeantc(theimage, theroi, zpos):
    # get an average timecourse from the voxels of an roi
    xstart = theroi[0][0]
    xend = theroi[1][0]
    ystart = theroi[0][1]
    yend = theroi[1][1]
    thesubreg = theimage[:, zpos, ystart:yend, xstart:xend]
    theshape = thesubreg.shape
    numtimepoints = theshape[0]
    themeans = np.zeros(numtimepoints)
    timeindex = np.arange(numtimepoints)
    for i in timeindex:
        themeans[i] = np.mean(np.ravel(thesubreg[i, :, :]))
    return themeans


def getroival(theimage, theroi, zpos):
    # get the average value from an roi in a 3D image
    xstart = theroi[0][0]
    xend = theroi[1][0]
    ystart = theroi[0][1]
    yend = theroi[1][1]
    theroival = np.mean(theimage[zpos, ystart:yend, xstart:xend])
    return theroival

def qualitytag(thestring, thequality):
    colors = ('00ff00', 'ffff00', 'ff0000')
    return '<FONT COLOR="{}">{}</FONT>'.format(colors[thequality], thestring)

