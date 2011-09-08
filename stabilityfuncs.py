#!/opt/local/bin/python
#
#       $Author: frederic $
#       $Date: 2011/03/25 15:06:35 $
#       $Id: stabilityfuncs.py,v 1.8 2011/03/25 15:06:35 frederic Exp $
#
import sys
import os
import time
import matplotlib
import numpy as np
import scipy as sp
import pylab as P
from nifti import *
from htmltagutils import *
from pylab import plot, legend, show, hold

########################################################################
#
#
#  Subroutine definitions
#
#
########################################################################
def setlimits(coil):
    limitdict={}

    # sample dependant quantities
    limitdict['DopedWaterPhantom_rad']=			((82.0,84.5),(81.0,86.5),(0,0),"Doped water phantom radius","")
    limitdict['DopedWaterPhantom_shape']=		((0.97,1.03),(0.95,1.05),(0,0),"Doped water phantom shape","")
    limitdict['DopedWaterPhantom_snr']=			((0.97,1.03),(0.95,1.05),(0,0),"Doped water phantom SNR","")	# do not yet have normative limits
    
    limitdict['BIRNphantom_rad']=			((80.9209448163338,86.6421984038845),(79.967402551742,87.5957406684762),(0,0),"BIRN phantom radius","")
    limitdict['BIRNphantom_shape']=			((0.941427953309648,1.05175367332324),(0.923040333307383,1.07014129332551),(0,0),"BIRN phantom shape","")
    limitdict['BIRNphantom_snr']=			((0.97,1.03),(0.95,1.05),(0,0),"BIRN phantom SNR","")	# do not yet have normative limits
    
    limitdict['head_rad']=				((60.0,70.0),(50.0,80.0),(0,0),"Head radius","")	# do not yet have normative limits
    limitdict['head_shape']=				((1.4,1.5),(1.2,1.7),(0,0),"Head shape","")	# do not yet have normative limits
    limitdict['head_snr']=				((1.4,1.5),(1.2,1.7),(0,0),"Head SNR","")	# do not yet have normative limits
    
    # coil independent quantities
    limitdict['center_of_mass_x']=			((30.0528824213927,34.071011383917),(29.383194260972,34.7406995443377),(1,0),"Center of mass x","")
    limitdict['center_of_mass_y']=			((28.0140073333838,32.9638687728109),(27.1890304268127,33.788845679382),(1,0),"Center of mass y","")
    limitdict['center_of_mass_z']=			((12.8178829000925,14.5672481911823),(12.5263220182442,14.8588090730306),(0,0),"Center of mass z","")

    # coil dependent quantities
    coilrecognized=0
    if (coil=='32Ch_Head'):
	coilrecognized=1
        limitdict['peripheral_angle_p-p%']=		((0,10.0),(0,10.0),(0,0),"Peripheral angle intensity p-p %","")
        limitdict['peripheral_angle_SFNR_p-p%']=	((0,10.0),(0,10.0),(0,0),"Peripheral angle SFNR p-p %","")

        limitdict['central_roi_mean']=			((100.0,200.0),(50.0,300.0),(0,0),"Central ROI mean","")
        limitdict['central_roi_raw_p-p%']=		((0.0,0.5),(0.0,0.6),(0,0),"Central ROI raw p-p %","")
        limitdict['central_roi_raw_std%']=		((0.0,0.15),(0.0,0.25),(0,0),"Central ROI raw stddev %","")
        limitdict['central_roi_detrended_mean']=	((600.0,900.0),(500.0,1200.0),(1,1),"Central ROI detrended mean","M")
        limitdict['central_roi_detrended_p-p%']=	((0.0,0.5),(0.0,0.6),(1,1),"Central ROI detrended p-p %","P")
        limitdict['central_roi_detrended_std%']=	((0.0,0.15),(0.0,0.25),(1,0),"Central ROI detrended stddev %","")
        limitdict['central_roi_SNR']=			((100.0,100000.0),(75.0,100000.0),(1,1),"Central ROI SNR","S")
        limitdict['central_roi_SFNR']=			((238.754898030312,365.602121158555),(217.613694175605,386.743325013262),(0,0),"Central ROI SFNR","F")
        limitdict['central_roi_polyfit_lin']=		((-0.00306061155726304,0.00338754176512665),(-0.00413530377766132,0.00446223398552493),(0,0),"Central ROI polyfit linear term","")
        limitdict['central_roi_polyfit_quad']=		((-5.21634333271468E-06,5.54843534508463E-06),(-7.0104731123479E-06,7.34256512471785E-06),(0,0),"Central ROI polyfit quadratic term","")
        limitdict['central_roi_drift%']=		((0.0,0.536),(0.0,0.688),(1,0),"Central ROI drift %","D")
    
        limitdict['peripheral_roi_mean']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI mean","")
        limitdict['peripheral_roi_raw_p-p%']=		((0.0,0.5),(0.0,0.6),(0,0),"Peripheral ROI raw p-p %","")
        limitdict['peripheral_roi_raw_std%']=		((0,0.125),(0,0.15),(0,0),"Peripheral ROI raw stddev %","")
        limitdict['peripheral_roi_detrended_p-p%']=	((0.0,0.5),(0.0,0.6),(1,1),"Peripheral ROI detrended p-p %","p")
        limitdict['peripheral_roi_detrended_std%']=	((0,0.125),(0,0.15),(0,0),"Peripheral ROI detrended stddev %","")
        limitdict['peripheral_roi_SNR']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI SNR","")
        limitdict['peripheral_roi_SFNR']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI SFNR","")
        limitdict['peripheral_roi_polyfit_lin']=	((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI polyfit linear term","")
        limitdict['peripheral_roi_polyfit_quad']=	((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI polyfit quadratic term","")
        limitdict['peripheral_roi_drift%']=		((0.0,0.397),(0.0,0.507),(1,0),"Peripheral ROI drift %","d")
    
        limitdict['odd_ghost_mean']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost mean","")
        limitdict['odd_ghost_std']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost stddev","")
        limitdict['odd_ghost_min']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost min","")
        limitdict['odd_ghost_max']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost max","")
        limitdict['odd_ghost_p-p']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost p-p","")
        limitdict['odd_ghost_p-p%']=			((0.0,5.0),(0.0,10.0),(1,0),"Odd ghost p-p %","")
        limitdict['even_ghost_mean']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost mean","")
        limitdict['even_ghost_std']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost stddev","")
        limitdict['even_ghost_min']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost min","")
        limitdict['even_ghost_max']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost max","")
        limitdict['even_ghost_p-p']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost p-p","")
        limitdict['even_ghost_p-p%']=			((0.0,5.0),(0.0,10.0),(1,0),"Even ghost p-p %","")
        limitdict['weissrdc']=				((0.0,5.0),(0.0,10.0),(1,0),"Weisskoff radius of decorrelation","W")
    if (coil=='HeadMatrix'):
	coilrecognized=1
        limitdict['peripheral_angle_p-p%']=		((0,10.0),(0,10.0),(0,0),"Peripheral angle intensity p-p %","")
        limitdict['peripheral_angle_SFNR_p-p%']=	((0,10.0),(0,10.0),(0,0),"Peripheral angle SFNR p-p %","")

        limitdict['central_roi_mean']=			((100.0,200.0),(50.0,300.0),(0,0),"Central ROI mean","")
        limitdict['central_roi_raw_p-p%']=		((0.0,0.5),(0.0,0.6),(0,0),"Central ROI raw p-p %","")
        limitdict['central_roi_raw_std%']=		((0.0,0.15),(0.0,0.25),(0,0),"Central ROI raw stddev %","")
        limitdict['central_roi_detrended_mean']=	((1100.0,1600.0),(900.0,1800.0),(1,1),"Central ROI detrended mean","M")
        limitdict['central_roi_detrended_p-p%']=	((0.0,0.5),(0.0,0.6),(1,1),"Central ROI detrended p-p %","P")
        limitdict['central_roi_detrended_std%']=	((0.0,0.15),(0.0,0.25),(1,0),"Central ROI detrended stddev %","")
        limitdict['central_roi_SNR']=			((250.0,100000.0),(200.0,100000.0),(1,1),"Central ROI SNR","S")
        limitdict['central_roi_SFNR']=			((238.754898030312,365.602121158555),(217.613694175605,386.743325013262),(0,0),"Central ROI SFNR","F")
        limitdict['central_roi_polyfit_lin']=		((-0.00306061155726304,0.00338754176512665),(-0.00413530377766132,0.00446223398552493),(0,0),"Central ROI polyfit linear term","")
        limitdict['central_roi_polyfit_quad']=		((-5.21634333271468E-06,5.54843534508463E-06),(-7.0104731123479E-06,7.34256512471785E-06),(0,0),"Central ROI polyfit quadratic term","")
        limitdict['central_roi_drift%']=		((0.0,0.536),(0.0,0.688),(1,0),"Central ROI drift %","D")
    
        limitdict['peripheral_roi_mean']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI mean","")
        limitdict['peripheral_roi_raw_p-p%']=		((0.0,0.5),(0.0,0.6),(0,0),"Peripheral ROI raw p-p %","")
        limitdict['peripheral_roi_raw_std%']=		((0,0.125),(0,0.15),(0,0),"Peripheral ROI raw stddev %","")
        limitdict['peripheral_roi_detrended_p-p%']=	((0.0,0.5),(0.0,0.6),(1,1),"Peripheral ROI detrended p-p %","p")
        limitdict['peripheral_roi_detrended_std%']=	((0,0.125),(0,0.15),(0,0),"Peripheral ROI detrended stddev %","")
        limitdict['peripheral_roi_SNR']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI SNR","")
        limitdict['peripheral_roi_SFNR']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI SFNR","")
        limitdict['peripheral_roi_polyfit_lin']=	((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI polyfit linear term","")
        limitdict['peripheral_roi_polyfit_quad']=	((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI polyfit quadratic term","")
        limitdict['peripheral_roi_drift%']=		((0.0,0.397),(0.0,0.507),(1,0),"Peripheral ROI drift %","d")
    
        limitdict['odd_ghost_mean']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost mean","")
        limitdict['odd_ghost_std']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost stddev","")
        limitdict['odd_ghost_min']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost min","")
        limitdict['odd_ghost_max']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost max","")
        limitdict['odd_ghost_p-p']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost p-p","")
        limitdict['odd_ghost_p-p%']=			((0.0,5.0),(0.0,10.0),(1,0),"Odd ghost p-p %","O")
        limitdict['even_ghost_mean']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost mean","")
        limitdict['even_ghost_std']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost stddev","")
        limitdict['even_ghost_min']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost min","")
        limitdict['even_ghost_max']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost max","")
        limitdict['even_ghost_p-p']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost p-p","")
        limitdict['even_ghost_p-p%']=			((0.0,5.0),(0.0,10.0),(1,0),"Even ghost p-p %","E")
        limitdict['weissrdc']=				((0.0,5.0),(0.0,10.0),(1,0),"Weisskoff radius of decorrelation","W")
    if (coil=='TxRx_Head'):
	coilrecognized=1
        limitdict['peripheral_angle_p-p%']=		((0,10.0),(0,10.0),(0,0),"Peripheral angle intensity p-p %","")
        limitdict['peripheral_angle_SFNR_p-p%']=	((0,10.0),(0,10.0),(0,0),"Peripheral angle SFNR p-p %","")

        limitdict['central_roi_mean']=			((100.0,200.0),(50.0,300.0),(0,0),"Central ROI mean","")
        limitdict['central_roi_raw_p-p%']=		((0.0,0.5),(0.0,0.6),(0,0),"Central ROI raw p-p %","")
        limitdict['central_roi_raw_std%']=		((0.0,0.15),(0.0,0.25),(0,0),"Central ROI raw stddev %","")
        limitdict['central_roi_detrended_mean']=	((1200.0,1750.0),(1000.0,2000.0),(1,1),"Central ROI detrended mean","M")
        limitdict['central_roi_detrended_p-p%']=	((0.0,0.5),(0.0,0.6),(1,1),"Central ROI detrended p-p %","P")
        limitdict['central_roi_detrended_std%']=	((0.0,0.15),(0.0,0.25),(1,0),"Central ROI detrended stddev %","")
        limitdict['central_roi_SNR']=			((300.0,100000.0),(250.0,100000.0),(1,1),"Central ROI SNR","S")
        limitdict['central_roi_SFNR']=			((238.754898030312,365.602121158555),(217.613694175605,386.743325013262),(0,0),"Central ROI SFNR","F")
        limitdict['central_roi_polyfit_lin']=		((-0.00306061155726304,0.00338754176512665),(-0.00413530377766132,0.00446223398552493),(0,0),"Central ROI polyfit linear term","")
        limitdict['central_roi_polyfit_quad']=		((-5.21634333271468E-06,5.54843534508463E-06),(-7.0104731123479E-06,7.34256512471785E-06),(0,0),"Central ROI polyfit quadratic term","")
        limitdict['central_roi_drift%']=		((0.0,0.536),(0.0,0.688),(1,0),"Central ROI drift %","D")
    
        limitdict['peripheral_roi_mean']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI mean","")
        limitdict['peripheral_roi_raw_p-p%']=		((0.0,0.5),(0.0,0.6),(0,0),"Peripheral ROI raw p-p %","")
        limitdict['peripheral_roi_raw_std%']=		((0,0.125),(0,0.15),(0,0),"Peripheral ROI raw stddev %","")
        limitdict['peripheral_roi_detrended_p-p%']=	((0.0,0.5),(0.0,0.6),(1,1),"Peripheral ROI detrended p-p %","p")
        limitdict['peripheral_roi_detrended_std%']=	((0,0.125),(0,0.15),(0,0),"Peripheral ROI detrended stddev %","")
        limitdict['peripheral_roi_SNR']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI SNR","s")
        limitdict['peripheral_roi_SFNR']=		((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI SFNR","f")
        limitdict['peripheral_roi_polyfit_lin']=	((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI polyfit linear term","")
        limitdict['peripheral_roi_polyfit_quad']=	((100.0,200.0),(50.0,300.0),(0,0),"Peripheral ROI polyfit quadratic term","")
        limitdict['peripheral_roi_drift%']=		((0.0,0.397),(0.0,0.507),(1,0),"Peripheral ROI drift %","d")
   	 
        limitdict['odd_ghost_mean']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost mean","")
        limitdict['odd_ghost_std']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost stddev","")
        limitdict['odd_ghost_min']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost min","")
        limitdict['odd_ghost_max']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost max","")
        limitdict['odd_ghost_p-p']=			((0.0,5.0),(0.0,10.0),(0,0),"Odd ghost p-p","")
        limitdict['odd_ghost_p-p%']=			((0.0,5.0),(0.0,10.0),(1,0),"Odd ghost p-p %","O")
        limitdict['even_ghost_mean']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost mean","")
        limitdict['even_ghost_std']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost stddev","")
        limitdict['even_ghost_min']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost min","")
        limitdict['even_ghost_max']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost max","")
        limitdict['even_ghost_p-p']=			((0.0,5.0),(0.0,10.0),(0,0),"Even ghost p-p","")
        limitdict['even_ghost_p-p%']=			((0.0,5.0),(0.0,10.0),(1,0),"Even ghost p-p %","E")
        limitdict['weissrdc']=				((0.0,5.0),(0.0,10.0),(1,0),"Weisskoff radius of decorrelation","W")
    if (coilrecognized!=1):
	print "setlimit: coil not recognized!"
	exit(1)

    return(limitdict)
    
def freqanalysis(thetimecourse):
    thefftsignal=abs(fft(thetimecourse))
    thelen=len(thefftsignal)
    thefftnoise=(thefftsignal[0:thelen-2]+thefftsignal[2:thelen])/2.0
    return()

def makemask(inputim,inputthresh,useabs):
    if (useabs < 1):
	#print "using relative threshold"
        thethreshval = getfracval(inputim,inputthresh)
        #print "%2.2f percent threshold at %2.2f" % (100.0*inputthresh,thethreshval)
    else:
	thethreshval = inputthresh
    themask = sp.where(inputim > thethreshval, 1.0, 0.0)
    return(themask)

def vecnorm(thevec):
    return(np.sqrt(np.square(thevec).sum()))

# format limits
def formatlimits(thelimits):
    limitdesc=thelimits[3]
    warnmin=str(thelimits[0][0])
    warnmax=str(thelimits[0][1])
    failmin=str(thelimits[1][0])
    failmax=str(thelimits[1][1])
    return("\""+limitdesc+"\","+failmin+","+warnmin+","+warnmax+","+failmax)

# check to see if a parameter falls within preset limits.
def limitcheck(thenumber,thelimits):
    retval=2	# start with the assumption that the data is bad
    if((float(thenumber)>=float(thelimits[1][0])) and (float(thenumber)<=float(thelimits[1][1]))):
	retval=1 # number falls within the warning limits
    if((float(thenumber)>=float(thelimits[0][0])) and (float(thenumber)<=float(thelimits[0][1]))):
	retval=0 # number falls within the good limits
    #print thelimits[1][0],thelimits[0][0],thenumber,thelimits[0][1],thelimits[1][1],"--->",retval
    return(retval)

def qualitytag(thestring,thequality):
    if(thequality == 0):
	return(greentag(thestring))
    if(thequality == 1):
	return(yellowtag(thestring))
    return(redtag(thestring))

# generate a table of weisskoff data
def weisstable(roiareas,weisscvs,projcvs):
    theshape=roiareas.shape
    numareas=theshape[0]
    
    tablestring=tablerowtag(
	tableentrytag("Region Size")+
	tableentrytag("Predicted Std Dev")+
	tableentrytag("Actual Std Dev")+
	tableentrytag("Ratio")
	)
    for i in range(0,numareas):
        tablestring=tablestring+tablerowtag(
	    tableentrytag("%d" % (roiareas[i]))+
	    tableentrytag("%.4f" %(projcvs[i]))+
	    tableentrytag("%.4f" %(weisscvs[i]))+
	    tableentrytag("%.4f" %(weisscvs[i]/projcvs[i])))
    return(smalltag(tablepropstag(tablestring,300,"center")))

# generate the polynomial fit timecourse from the coefficients
def trendgen(thexvals,thefitcoffs):
    theshape=thefitcoffs.shape
    order = theshape[0]-1
    #print "fitting to order "+str(order)
    thepoly=thexvals
    thefit=0.0*thexvals
    if order>0:
	for i in range(1,order+1):
	    #print "fitting component "+str(i)+", coff="+str(thefitcoffs[order-i])
    	    thefit = thefit + thefitcoffs[order-i]*thepoly
	    thepoly = np.multiply(thepoly, thexvals)
    return(thefit)

# calculate the robust range of the all voxels
def completerobust(thearray):
    themin=getfracval(thearray,0.02)
    themax=getfracval(thearray,0.98)
    return([themin,themax])

# calculate the robust range of the non-zero voxels
def nzrobust(thearray):
    themin=getnzfracval(thearray,0.02)
    themax=getnzfracval(thearray,0.98)
    return([themin,themax])

# calculate the min and max of the non-zero voxels
def nzminmax(thearray):
    flatarray=np.ravel(thearray)
    nzindices=np.nonzero(flatarray)
    theflatarray = flatarray[nzindices]
    themax = np.max(theflatarray)
    themin = np.min(theflatarray)
    return([themin,themax])

# calculate the stats of the non-zero voxels
def completestats(thearray):
    themean = np.mean(thearray)
    thestddev = np.std(thearray)
    thevar = np.var(thearray)
    themax = np.max(thearray)
    themin = np.min(thearray)
    theptp = np.ptp(thearray)
    return([themean,thestddev,thevar,themax,themin,theptp])

# calculate the stats of the non-zero voxels
def nzstats(thearray):
    flatarray=np.ravel(thearray)
    nzindices=np.nonzero(flatarray)
    theflatarray = flatarray[nzindices]
    themean = np.mean(theflatarray)
    thestddev = np.std(theflatarray)
    thevar = np.var(theflatarray)
    themax = np.max(theflatarray)
    themin = np.min(theflatarray)
    theptp = np.ptp(theflatarray)
    return([themean,thestddev,thevar,themax,themin,theptp])

def showstats(thestats):
    formatstring = "mean = %2.2f, stddev = %2.2f, max = %2.2f, min = %2.2f"
    interpstring = (thestats[0],thestats[1],thestats[3],thestats[4])
    return(formatstring % interpstring)
    
# calculate the mean of the non-zero voxels
def nzmean(thearray):
    flatarray=np.ravel(thearray)
    nzindices=np.nonzero(flatarray)
    return(np.mean(flatarray[nzindices]))
    
# calculate the sum of an array across space
def arrayspatialsum(thearray):
    return(np.sum(thearray))
    
# show an roi timecourse plot
def showtc(thexvals,theyvals,thelabel):
    w, h = P.figaspect(0.25)
    roiplot = P.figure(figsize=(w,h))
    roisubplot = roiplot.add_subplot(111)
    roisubplot.plot(thexvals, theyvals, 'b')
    roisubplot.grid(True)
    #roisubplot.axes.Subplot.set_pad(0.1)
    for tick in roisubplot.xaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    for tick in roisubplot.yaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    roisubplot.set_title(thelabel,fontsize=30)
    return()

# show an roi timecourse plot and a fit line
def showvals(xvecs,yvecs,legendvec,specvals,thelabel,dolegend):
    numxs=len(xvecs)
    numys=len(yvecs)
    numlegends=len(legendvec)
    numspecvals=len(specvals)
    if (numxs!=numys) or (numxs!=numlegends) or (numxs!=numspecvals):
	print "dimensions do not match"
	exit(1)
    w, h = P.figaspect(0.50)
    roiplot = P.figure(figsize=(w,h))
    roisubplot = roiplot.add_subplot(111)
    if numys==1:
        roisubplot.plot(xvecs[0], yvecs[0], specvals[0])
        hold(True)
        if dolegend:
            legend(legendvec)
        hold(False)
    if numys==2:
        roisubplot.plot(xvecs[0], yvecs[0], specvals[0], xvecs[1], yvecs[1], specvals[1])
        hold(True)
        if dolegend:
            legend(legendvec)
        hold(False)
    if numys==3:
        roisubplot.plot(xvecs[0], yvecs[0], specvals[0], xvecs[1], yvecs[1], specvals[1], xvecs[2], yvecs[2], specvals[2])
        hold(True)
        if dolegend:
            legend(legendvec)
        hold(False)
    if numys==4:
        roisubplot.plot(xvecs[0], yvecs[0], specvals[0], xvecs[1], yvecs[1], specvals[1], xvecs[2], yvecs[2], specvals[2], xvecs[3], yvecs[3], specvals[3])
        hold(True)
        if dolegend:
            legend(legendvec)
        hold(False)
    if numys==5:
        roisubplot.plot(xvecs[0], yvecs[0], specvals[0], xvecs[1], yvecs[1], specvals[1], xvecs[2], yvecs[2], specvals[2], xvecs[3], yvecs[3], specvals[3], xvecs[4], yvecs[4], specvals[4])
        hold(True)
        if dolegend:
            legend(legendvec)
        hold(False)
    if numys==6:
        roisubplot.plot(xvecs[0], yvecs[0], specvals[0], xvecs[1], yvecs[1], specvals[1], xvecs[2], yvecs[2], specvals[2], xvecs[3], yvecs[3], specvals[3], xvecs[4], yvecs[4], specvals[4], xvecs[5], yvecs[5], specvals[5])
        hold(True)
        if dolegend:
            legend(legendvec)
        hold(False)
    roisubplot.grid(True)
    for tick in roisubplot.xaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    for tick in roisubplot.yaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    roisubplot.set_title(thelabel,fontsize=30)
    return()

# show an roi timecourse plot and a fit line
def showtc2(thexvals,theyvals,thefitvals,thelabel):
    w, h = P.figaspect(0.25)
    roiplot = P.figure(figsize=(w,h))
    roisubplot = roiplot.add_subplot(111)
    roisubplot.plot(thexvals, theyvals, 'b', thexvals, thefitvals, 'g')
    roisubplot.grid(True)
    #roisubplot.axes.Subplot.set_pad(0.1)
    for tick in roisubplot.xaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    for tick in roisubplot.yaxis.get_major_ticks():
        tick.label1.set_fontsize(20)
    roisubplot.set_title(thelabel,fontsize=30)
    return()

# initialize and show a loglog Weiskoff plot
def showweisskoff(theareas,thestddevs,theprojstddevs,thelabel):
    w, h = P.figaspect(1.0)
    roiplot = P.figure(figsize=(w,h))
    roiplot.subplots_adjust(hspace=0.35)
    roisubplot = roiplot.add_subplot(111)
    thestddevs=thestddevs+0.00000001
    roisubplot.loglog(theareas, thestddevs, 'r', theareas, theprojstddevs, 'k', basex=10)
    roisubplot.grid(True)
    #roiplot.title(thelabel)
    return()

# initialize and show a 2D slice from a dataset in greyscale
def showslice2(thedata,thelabel,minval,maxval,colormap):
    theshape=thedata.shape
    numslices=theshape[0]
    ysize=theshape[1]
    xsize=theshape[2]
    slicesqrt=int(np.ceil(np.sqrt(numslices)))
    theslice=np.zeros((ysize*slicesqrt,xsize*slicesqrt))
    for i in range(0,numslices):
        ypos=int(i/slicesqrt)*ysize
        xpos=int(i%slicesqrt)*xsize
        theslice[ypos:ypos+ysize,xpos:xpos+xsize]=thedata[i,:,:]
    if P.isinteractive():
	P.ioff()
    P.axis('off')
    P.axis('equal')
    P.subplots_adjust(hspace=0.0)
    P.axes([0,0,1,1], frameon = False)
    if (colormap==0):
        thecmap=P.cm.gray
    else:
	mycmdata1 = {
	    'red'  :  ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 1., 1.)),
	    'green':  ((0., 0., 0.), (0.5, 1.0, 1.0), (1., 0., 0.)),
	    'blue' :  ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 0., 0.))
	    }
	thecmap = P.matplotlib.colors.LinearSegmentedColormap('mycm', mycmdata1)
	#thecmap=P.cm.spectral
    theimptr = P.imshow(theslice, vmin=minval, vmax=maxval, interpolation='nearest', label=thelabel, aspect='equal', cmap=thecmap)
    #P.colorbar()
    return()

# initialize and show a 2D slice from a dataset in greyscale
def showslice3(thedata,thelabel,minval,maxval,colormap):
    theshape=thedata.shape
    ysize=theshape[0]
    xsize=theshape[1]
    theslice=np.zeros((ysize,xsize))
    if P.isinteractive():
	P.ioff()
    P.axis('off')
    P.axis('equal')
    P.subplots_adjust(hspace=0.0)
    P.axes([0,0,1,1], frameon = False)
    if (colormap==0):
        thecmap=P.cm.gray
    else:
	mycmdata1 = {
	    'red'  :  ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 1., 1.)),
	    'green':  ((0., 0., 0.), (0.5, 1.0, 1.0), (1., 0., 0.)),
	    'blue' :  ((0., 0., 0.), (0.5, 1.0, 0.0), (1., 0., 0.))
	    }
	thecmap = P.matplotlib.colors.LinearSegmentedColormap('mycm', mycmdata1)
	#thecmap=P.cm.spectral
    theimptr = P.imshow(thedata, vmin=minval, vmax=maxval, interpolation='nearest', label=thelabel, aspect='equal', cmap=thecmap)
    #P.colorbar()
    return()

# show a 2D slice from a dataset in greyscale
def showslice(theslice):
    if P.isinteractive():
	P.ioff()
    P.axis('off')
    P.axis('equal')
    P.axis('tight')
    P.imshow(theslice, interpolation='nearest', aspect='equal', cmap=P.cm.gray)
    P.colorbar()
    return()

def smooth(x,window_len=11,window='hanning'):
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
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

# Find the image intensity value that cleanly separates background from image
def findsepval(datamat):
    numbins=200
    themax = datamat.max()
    themin = datamat.min()
    (meanhist,bins) = np.histogram(datamat,bins=numbins,range=(themin,themax))
    smoothhist = smooth(meanhist)
    currentpos = int(numbins*0.05)
    minval=smoothhist[currentpos]
    for i in range(currentpos+1,numbins):
	if(smoothhist[i] < smoothhist[currentpos]):
	    currentpos=i
	if(smoothhist[i] > 1.2 * smoothhist[currentpos]):
	    break
    cummeanhist = np.cumsum(meanhist)
    #print "curpos %d, cummeanhist[curpos] %2.2f, cummeanhist[numbins-1] %d" % (currentpos, cummeanhist[currentpos], cummeanhist[numbins-1])
    cummeanhist[currentpos]
    cumfrac=(1.0*cummeanhist[currentpos])/(1.0*cummeanhist[numbins-1])
    sepval=bins[currentpos]
    return([sepval,cumfrac])

# Find the image intensity value which thefrac of the non-zero voxels in the image exceed
def getfracval(datamat,thefrac):
    numbins=200
    themax = datamat.max()
    themin = datamat.min()
    (meanhist,bins) = np.histogram(datamat,bins=numbins,range=(themin,themax))
    cummeanhist = np.cumsum(meanhist)
    target = cummeanhist[numbins-1]*thefrac
    for i in range(0,numbins):
	if cummeanhist[i]>=target:
	    return(bins[i])
    return(0.0)

# Find the image intensity value which thefrac of the non-zero voxels in the image exceed
def getnzfracval(datamat,thefrac):
    numbins=200
    (themin,themax) = nzminmax(datamat)
    (meanhist,bins) = np.histogram(datamat,bins=numbins,range=(themin,themax))
    cummeanhist = np.cumsum(meanhist)
    target = cummeanhist[numbins-1]*thefrac
    for i in range(0,numbins):
	if cummeanhist[i]>=target:
	    return(bins[i])
    return(0.0)

# find the center of mass of a 2D or 3D image
def findCOM(datamat):
    Mx = 0.0
    My = 0.0
    Mz = 0.0
    mass = 0.0
    val = 0.0
    arrdims=np.shape(datamat)
 
    if datamat.ndim==2:
        for i in range(0,arrdims[0]):
            for j in range(0,arrdims[1]):
		val = datamat[i,j]
            	My += (i * val)
            	Mx += (j * val)
                mass += val
        COM = (Mx/mass , My/mass, 0.0)
    if datamat.ndim==3:
        for i in range(0,arrdims[0]):
            for j in range(0,arrdims[1]):
                for k in range(0,arrdims[2]):
		    val = datamat[i,j,k]
            	    Mz += (i * val)
            	    My += (j * val)
            	    Mx += (k * val)
                    mass += val
        COM = (Mx/mass , My/mass, Mz/mass)
    return COM

# given an roi and a position, mark an roi
def markroi(theinputroi,zpos,roislice,theval):
    xstart=theinputroi[0][0]
    xend=theinputroi[1][0]
    ystart=theinputroi[0][1]
    yend=theinputroi[1][1]
    roislice[zpos,ystart:yend,xstart:xend]=theval
    return

# given a location and a size, define the corners of an roi
def setroilims(xpos,ypos,size):
    if (size%2)==0:
	halfsize=size/2
    	return(((int(round(xpos-halfsize)),int(round(ypos-halfsize))),
	    (int(round(xpos+halfsize)),int(round(ypos+halfsize)))))
    else:
	halfsize=(size-1)/2
    	return(((int(round(xpos-halfsize)),int(round(ypos-halfsize))),
	    (int(round(xpos+halfsize+1)),int(round(ypos+halfsize+1)))))

# get an snr timecourse from the voxels of an roi
def getroisnr(theimage,theroi,zpos):
    xstart=theroi[0][0]
    xend=theroi[1][0]
    ystart=theroi[0][1]
    yend=theroi[1][1]
    thesubreg=theimage[:,zpos,ystart:yend,xstart:xend]
    theshape=thesubreg.shape
    numtimepoints=theshape[0]
    themeans=np.zeros(numtimepoints)
    thestddevs=np.zeros(numtimepoints)
    themax=np.zeros(numtimepoints)
    themin=np.zeros(numtimepoints)
    thesnrs=np.zeros(numtimepoints)
    timeindex=np.arange(0,numtimepoints)
    for i in timeindex:
        themeans[i]=np.mean(np.ravel(thesubreg[i,:,:]))
        thestddevs[i]=np.std(np.ravel(thesubreg[i,:,:]))
        themax[i]=np.max(np.ravel(thesubreg[i,:,:]))
        themin[i]=np.min(np.ravel(thesubreg[i,:,:]))
        thesnrs[i]=themeans[i]/thestddevs[i]
    return(thesnrs)

# get all the voxels from an roi and return a 2d (time by space) array
def getroivoxels(theimage,theroi,zpos):
    xstart=theroi[0][0]
    xend=theroi[1][0]
    ystart=theroi[0][1]
    yend=theroi[1][1]
    thesubreg=theimage[:,zpos,ystart:yend,xstart:xend]
    theshape=thesubreg.shape
    numtimepoints=theshape[0]
    thevoxels=np.zeros((numtimepoints,theshape[1]*theshape[2]))
    timeindex=np.arange(0,numtimepoints)
    for i in timeindex:
        thevoxels[i,:]=np.ravel(thesubreg[i,:,:])
    return(thevoxels)

# get a standard deviation timecourse from the voxels of an roi
def getroistdtc(theimage,theroi,zpos):
    xstart=theroi[0][0]
    xend=theroi[1][0]
    ystart=theroi[0][1]
    yend=theroi[1][1]
    thesubreg=theimage[:,zpos,ystart:yend,xstart:xend]
    theshape=thesubreg.shape
    numtimepoints=theshape[0]
    thestds=np.zeros(numtimepoints)
    timeindex=np.arange(0,numtimepoints)
    for i in timeindex:
        thestds[i]=np.std(np.ravel(thesubreg[i,:,:]))
    return(thestds)

# get an average timecourse from the voxels of an roi
def getroimeantc(theimage,theroi,zpos):
    xstart=theroi[0][0]
    xend=theroi[1][0]
    ystart=theroi[0][1]
    yend=theroi[1][1]
    thesubreg=theimage[:,zpos,ystart:yend,xstart:xend]
    theshape=thesubreg.shape
    numtimepoints=theshape[0]
    themeans=np.zeros(numtimepoints)
    timeindex=np.arange(0,numtimepoints)
    for i in timeindex:
        themeans[i]=np.mean(np.ravel(thesubreg[i,:,:]))
    return(themeans)

# get the average value from an roi in a 3D image
def getroival(theimage,theroi,zpos):
    xstart=theroi[0][0]
    xend=theroi[1][0]
    ystart=theroi[0][1]
    yend=theroi[1][1]
    theroival=np.mean(theimage[zpos,ystart:yend,xstart:xend])
    return(theroival)

# make a captioned image with statistics
def makecaptionedimage(imagetitle,thestats,imagename,thewidth):
    if(thestats==[]):
        imcapstring = paratag(boldtag(imagetitle))
    else:
        imcapstring = paratag(boldtag(imagetitle) + breaktag(showstats(thestats)))
    return(imcapstring + imagetag(imagename,thewidth))

# send a command to the shell
def doashellcmd(cmd):
    a = os.popen(cmd)
    while 1:
        line = a.readline()
        if not line: break
        retval = line[:-1]
        return retval
