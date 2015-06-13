#!/usr/bin/env python

import os
from os.path import join as pjoin
import shutil
import time
import logging
import subprocess

from htmltagutils import *
import stabilityfuncs as sf


def stability_eval(specs, thedictionary):
    stability_results = {}
    for theentry in thedictionary:
        # see if this spec is tracked
        try:
            entrydata = specs[theentry]
            if entrydata[2][0] == 1:
                stability_results[theentry] = {}
                stability_results[theentry]['key'] = theentry
                stability_results[theentry]['varname'] = entrydata[3]
                stability_results[theentry]['value'] = thedictionary[theentry]
                stability_results[theentry]['warnrange'] = entrydata[0]
                stability_results[theentry]['failrange'] = entrydata[1]
                stability_results[theentry]['quality'] = sf.limitcheck(thedictionary[theentry], specs[theentry])
                if entrydata[2][1] == 1:
                    stability_results[theentry]['critical'] = True
                else:
                    stability_results[theentry]['critical'] = False
                stability_results[theentry]['symbol'] = entrydata[4]
        except KeyError:
            pass

    return stability_results


def stabilitysummary(datadirectory, outputdirectory, whichscan, TargetisBIRNphantom):
    # initialize the outut directory if need be
    if not os.path.exists(pjoin(outputdirectory, whichscan)):
        os.makedirs(pjoin(outputdirectory, whichscan))

    ########################################################################
    #
    #
    # scan the data directory for stability scans
    #
    stabilitydirs = os.listdir(datadirectory)
    stabilitydirs = sorted([filename for filename in stabilitydirs if filename.startswith("stability_")])

    # print stabilitydirs

    #
    # pull all the data files into a dictionary array
    #
    datadict = {}
    filenumber_TARGET = 0
    num_cp_TARGET = 0
    num_12_TARGET = 0
    num_32_TARGET = 0
    for summaryfile in stabilitydirs:
        logging.info('Beginning processing for ' + summaryfile)
        datadict[filenumber_TARGET] = {}
        try:
            datadict[filenumber_TARGET]['datadir'] = pjoin(summaryfile, whichscan, 'procresults')
            try:
                datadict[filenumber_TARGET].update(sf.dict_from_tsvfile(pjoin(datadirectory, datadict[filenumber_TARGET]['datadir'], 'analysissummary.txt')))
                ObjectisBIRNphantom = (datadict[filenumber_TARGET]['Object'] == 'BIRN phantom')
                if ObjectisBIRNphantom == TargetisBIRNphantom:
                    if datadict[filenumber_TARGET]['Coil'] == 'TxRx_Head':
                        num_cp_TARGET += 1
                    if datadict[filenumber_TARGET]['Coil'] == '32Ch_Head':
                        num_32_TARGET += 1
                    if datadict[filenumber_TARGET]['Coil'] == 'HeadMatrix':
                        num_12_TARGET += 1
                    filenumber_TARGET += 1
            except IOError:
                pass
        except KeyError:
            pass
    logging.debug("{} CP coil runs ({} phantom)".format(num_cp_TARGET, 'BIRN' if TargetisBIRNphantom else 'NONBIRN'))
    logging.debug("{} 12 channel coil runs ({} phantom)".format(num_12_TARGET, 'BIRN' if TargetisBIRNphantom else 'NONBIRN'))
    logging.debug("{} 32 channel coil runs ({} phantom)".format(num_32_TARGET, 'BIRN' if TargetisBIRNphantom else 'NONBIRN'))

    #######################################################################################
    #
    # sort the data up by coil and write to files
    #
    mostrecenttimes = {}
    for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
        with open(pjoin(outputdirectory, whichscan, targetcoil + '_vals.txt'), 'w') as fp:
            for i in range(filenumber_TARGET):
                ObjectisBIRNphantom = (datadict[i]['Object'] == 'BIRN phantom')
                if datadict[i]['Coil'] == targetcoil and ObjectisBIRNphantom == TargetisBIRNphantom:
                    try:
                        if datadict[i]['Protocol'] != 'nothing':
                            mostrecenttimes[targetcoil] = datadict[i]['DateTime']
                            fp.write(' '.join([datadict[i][k] for k in ('Coil',
                                         'DateTime',
                                         'central_roi_detrended_p-p%',
                                         'peripheral_roi_detrended_p-p%',
                                         'central_roi_SNR',
                                         'peripheral_roi_SNR',
                                         'central_roi_SFNR',
                                         'peripheral_roi_SFNR',
                                         'odd_ghost_mean',
                                         'odd_ghost_max',
                                         'odd_ghost_min',
                                         'even_ghost_mean',
                                         'even_ghost_max',
                                         'even_ghost_min',
                                         'object_radius_mm',
                                         'object_shape',
                                         'center_of_mass_x',
                                         'center_of_mass_y',
                                         'center_of_mass_z',
                                         'central_roi_detrended_mean',
                                         'central_roi_drift%',
                                         'peripheral_roi_drift%',
                                         'weissrdc',
                                         'central_roi_detrended_mean',)]))
                            fp.write('\n')
                    except KeyError:
                        pass

                        

    #######################################################################################
    # generate plot control files to graph all interesting stability parameters
    #
    # central and peripheral SNR and SFNR
    #

    if TargetisBIRNphantom:
        wlp = 'points'
    else:
        wlp = 'linespoints'

    pointsize = 'ps 2'
    with open(outputdirectory + "/" + whichscan + "/plotcmds_snrsfnr", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [0:800]\n")
        fp.writelines("set title \"Absolute SNR and SFNR\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"p-p percent\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:5 title 'Central ROI SNR' with " + wlp + " " + pointsize + " pt 3, \\\n")
        fp.writelines(
            "\"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:6 title 'Peripheral ROI SNR' with " + wlp + " " + pointsize + " pt 4, \\\n")
        fp.writelines(
            "\"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:7 title 'Central ROI SFNR' with " + wlp + " " + pointsize + " pt 4, \\\n")
        fp.writelines(
            "\"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:8 title 'Peripheral ROI SFNR' with " + wlp + " " + pointsize + " pt 4\n")

    #
    # central and peripheral stability
    #
    with open(outputdirectory + "/" + whichscan + "/plotcmds_roistab", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [0:1.0]\n")
        fp.writelines("set title \"ROI stability: p-p variation in percent\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"p-p percent\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:3 title 'Central ROI' with " + wlp + " " + pointsize + " pt 3, \\\n")
        fp.writelines(
            "\"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:4 title 'Peripheral ROI' with " + wlp + " " + pointsize + " pt 4\n")

    #
    # central and peripheral drift
    #
    with open(outputdirectory + "/" + whichscan + "/plotcmds_roidrift", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [0:1.0]\n")
        fp.writelines("set title \"ROI linear and quadratic drift: p-p amplitude in percent\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"p-p percent\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:21 title 'Central ROI' with " + wlp + " " + pointsize + " pt 3, \\\n")
        fp.writelines(
            "\"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:22 title 'Peripheral ROI' with " + wlp + " " + pointsize + " pt 4\n")

    #
    # ghosts
    #
    with open(outputdirectory + "/" + whichscan + "/plotcmds_ghost", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [0:15.0]\n")
        fp.writelines("set title \"Ghost percentage\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"Ghost amplitude (%)\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:9:11:10 title 'Odd ghost' with yerrorbars ps 1 pt 3, \\\n")
        fp.writelines(
            "\"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:12:14:13 title 'Even ghost' with yerrorbars ps 1 pt 4\n")

    #
    # object radius
    #
    with open(outputdirectory + "/" + whichscan + "/plotcmds_objradius", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [75.0:90.0]\n")
        fp.writelines("set title \"Phantom radius\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"Phantom radius (mm)\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:15 title 'Object radius in mm' with " + wlp + " " + pointsize + " pt 3\n")

    #
    # object shape
    #
    with open(outputdirectory + "/" + whichscan + "/plotcmds_objshape", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [0.9:1.1]\n")
        fp.writelines("set title \"Phantom shape (y/x)\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"Phantom shape (y/x)\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:16 title 'Object shape (y/x)' with " + wlp + " " + pointsize + " pt 3\n")

    #
    # weisskoff
    #
    with open(outputdirectory + "/" + whichscan + "/plotcmds_weissrdc", "w") as fp:
        fp.writelines(
            "set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt \"%Y%m%dT%H:%M:%S\"\n")
        fp.writelines("set xrange [\"20091120T00:00:00\":]\n")
        fp.writelines("set yrange [0:12.0]\n")
        fp.writelines("set title \"Weisskoff Radius of decorrelation\"\n")
        fp.writelines("set xlabel \"Date\"\n")
        fp.writelines("set ylabel \"Weisskoff RDC\"\n")
        fp.writelines(
            "plot    \"" + outputdirectory + "/" + whichscan + "/graphtemp\" using 2:23 title 'Weisskoff RDC' with " + wlp + " " + pointsize + " pt 3\n")

    #
    # central mean signal intensity
    #
    with open(pjoin(outputdirectory, whichscan, 'plotcmds_centralsignal'), 'w') as fp:
        fp.writelines(
            'set terminal jpeg\n set autoscale\n unset log\n unset label\n set xdata time\n set xtics autofreq rotate\n set ytic auto\n set timefmt "%Y%m%dT%H:%M:%S"\n')
        fp.writelines('set xrange ["20091120T00:00:00":]\n')
        fp.writelines('set yrange [0:2000.0]\n')
        fp.writelines('set title "Detrended mean central signal intensity"\n')
        fp.writelines('set xlabel "Date"\n')
        fp.writelines('set ylabel "Mean signal intensity"\n')
        fp.writelines(
            'plot    "' + outputdirectory + '/' + whichscan + '/graphtemp" using 2:24 title \'Mean central signal intensity\' with ' + wlp + ' ' + pointsize + ' pt 3\n')

    for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
        outscandir = pjoin(outputdirectory, whichscan)
        shutil.copyfile(pjoin(outscandir, targetcoil + '_vals.txt'), pjoin(outscandir, 'graphtemp'))
        for plottype in ['snrsfnr', 'roistab', 'roidrift', 'ghost', 'objradius', 'objshape', 'weissrdc',
                         'centralsignal']:
            subprocess.Popen(['gnuplot', pjoin(outscandir, 'plotcmds_' + plottype)],
                             stdout=open(pjoin(outscandir, '{}_{}.jpg'.format(targetcoil, plottype)), 'wb'))

    #######################################################################################
    # generate a report file
    #
    thisdate = time.strftime("%m/%d/%Y %H:%M:%S", time.localtime())
    with open(outputdirectory + "/" + whichscan + "/" + "stabilityreport.html", "w") as fp:

        fp.writelines(
            "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n\n")

        fp.writelines("<head>\n")
        fp.writelines("<script type=\"text/javascript\" language=\"javascript1.5\">\n")
        fp.writelines("    <!-- Hide script from old browsers\n")
        fp.writelines("\n")
        fp.writelines("        window.onload = theClock\n")
        fp.writelines("\n")
        fp.writelines("        function theClock() {\n")
        fp.writelines("            var one_minute=1000*60;\n")
        fp.writelines("            var one_hour=one_minute*60;\n")
        fp.writelines("            var one_day=one_hour*24;\n")
        # format: 20150410T080041
        date32 = mostrecenttimes.get('32Ch_Head', '19700101T000000')
        date12 = mostrecenttimes.get('HeadMatrix', '19700101T000000')
        datecp = mostrecenttimes.get('TxRx_Head', '19700101T000000')
        args32 = date32[0:4] + "," + str(int(date32[4:6]) - 1) + "," + date32[6:8] + "," + date32[9:11] + "," + date32[
                                                                                                                11:13] + "," + date32[
                                                                                                                               13:15]
        args12 = date12[0:4] + "," + str(int(date12[4:6]) - 1) + "," + date12[6:8] + "," + date12[9:11] + "," + date12[
                                                                                                                11:13] + "," + date12[
                                                                                                                               13:15]
        argscp = datecp[0:4] + "," + str(int(datecp[4:6]) - 1) + "," + datecp[6:8] + "," + datecp[9:11] + "," + datecp[
                                                                                                                11:13] + "," + datecp[
                                                                                                                               13:15]
        fp.writelines("            head32time = new Date(" + args32 + ");\n")
        fp.writelines("            head12time = new Date(" + args12 + ");\n")
        fp.writelines("            headcptime = new Date(" + argscp + ");\n")
        fp.writelines("            now = new Date();\n")
        fp.writelines("            var time32=now-head32time;\n")
        fp.writelines("            var time12=now-head12time;\n")
        fp.writelines("            var timecp=now-headcptime;\n")
        fp.writelines("            var days32=Math.floor(time32/one_day);\n")
        fp.writelines("            var days12=Math.floor(time12/one_day);\n")
        fp.writelines("            var dayscp=Math.floor(timecp/one_day);\n")
        fp.writelines("            var hours32=Math.floor((time32-days32*one_day)/one_hour);\n")
        fp.writelines("            var hours12=Math.floor((time12-days12*one_day)/one_hour);\n")
        fp.writelines("            var hourscp=Math.floor((timecp-dayscp*one_day)/one_hour);\n")
        fp.writelines(
            "            var minutes32=Math.floor((time32-days32*one_day-hours32*one_hour)/one_minute);\n")
        fp.writelines(
            "            var minutes12=Math.floor((time12-days12*one_day-hours12*one_hour)/one_minute);\n")
        fp.writelines(
            "            var minutescp=Math.floor((timecp-dayscp*one_day-hourscp*one_hour)/one_minute);\n")
        fp.writelines("            var mindays=Math.min(dayscp,days12,days32);\n")
        fp.writelines("\n")
        fp.writelines("            fmttime32=days32+\"d \"+hours32+\"h \"+minutes32+\"m\";\n")
        fp.writelines("            fmttime12=days12+\"d \"+hours12+\"h \"+minutes12+\"m\";\n")
        fp.writelines("            fmttimecp=dayscp+\"d \"+hourscp+\"h \"+minutescp+\"m\";\n")
        fp.writelines("\n")
        fp.writelines("            if(mindays>4)\n")
        fp.writelines("                {\n")
        fp.writelines("                document.body.style.backgroundColor=\'#ffaaaa\';\n")
        fp.writelines("                }\n")
        fp.writelines("            else\n")
        fp.writelines("                {\n")
        fp.writelines("                document.body.style.backgroundColor=\'#cccccc\';\n")
        fp.writelines("                }\n")
        fp.writelines("\n")
        fp.writelines("            time32Span = document.getElementById(\"formtime32\");\n")
        fp.writelines(
            "            time32Span.replaceChild(document.createTextNode(fmttime32), time32Span.firstChild);\n")
        fp.writelines("\n")
        fp.writelines("            time12Span = document.getElementById(\"formtime12\");\n")
        fp.writelines(
            "            time12Span.replaceChild(document.createTextNode(fmttime12), time12Span.firstChild);\n")
        fp.writelines("\n")
        fp.writelines("            timecpSpan = document.getElementById(\"formtimecp\");\n")
        fp.writelines(
            "            timecpSpan.replaceChild(document.createTextNode(fmttimecp), timecpSpan.firstChild);\n")
        fp.writelines("\n")
        fp.writelines("            setTimeout(\"theClock()\",1000);\n")
        fp.writelines("        }\n")
        fp.writelines("\n")
        fp.writelines("    // End hiding script from old browsers -->\n")
        fp.writelines("</script>\n")

        fp.writelines("<title>Stability summary generated on " + thisdate + "</title>\n")
        fp.writelines("<style type=\"text/css\">\n")
        fp.writelines("<meta http-equiv=\"Refresh\" content=\"600\">\n")
        fp.writelines("h1 {font-family:courier new;text-decoration:underline;}\n")
        fp.writelines("h2 {font-family:courier new;color: teal; text-decoration:underline;}\n")
        fp.writelines("h3 {font-family:courier new;color: maroon; text-decoration:none;}\n")
        fp.writelines("h4 {font-family:courier new;text-decoration:none;}\n")
        fp.writelines("p {font-family:courier new;color:black; font-size:16px;text-decoration:none;}\n")
        fp.writelines("td {font-family:courier new;color:black; font-size:12px;text-decoration:none;}\n")
        fp.writelines("</style>\n")
        fp.writelines("</head>\n\n")

        fp.writelines("<body>\n")

        # Compose the image table
        myimwidth = 400
        imagehdrstr = bigheadertag(
            "{} phantom stability tests".format('BIRN' if TargetisBIRNphantom else 'NONBIRN')) + headertag("Summary images")

        ROIstab_headerstr = headertag("ROI p-p% variation:")
        ROIdrift_headerstr = headertag("ROI lin-quad drift % amplitude:")
        snrsfnr_headerstr = headertag("SNR, SFNR:")
        ghostamp_headerstr = headertag("Ghost amplitude:")
        objradius_headerstr = headertag("Object radius:")
        objshape_headerstr = headertag("Object shape:")
        weissrdc_headerstr = headertag("Weisskoff RDC:")
        centralsignal_headerstr = headertag("Mean central signal intensity:")

        thecpinfostring = bigheadertag("CP TxRx") + headertag(
            "most recent scan: " + breaktag("<span id=\"formtimecp\">?</span>\n"))
        thecproi_imagestr = imagetag("TxRx_Head_roistab.jpg", myimwidth)
        thecproidrift_imagestr = imagetag("TxRx_Head_roidrift.jpg", myimwidth)
        thecpsnrsfnr_imagestr = imagetag("TxRx_Head_snrsfnr.jpg", myimwidth)
        thecpghost_imagestr = imagetag("TxRx_Head_ghost.jpg", myimwidth)
        thecpobjradius_imagestr = imagetag("TxRx_Head_objradius.jpg", myimwidth)
        thecpobjshape_imagestr = imagetag("TxRx_Head_objshape.jpg", myimwidth)
        thecpweissrdc_imagestr = imagetag("TxRx_Head_weissrdc.jpg", myimwidth)
        thecpcentralsignal_imagestr = imagetag("TxRx_Head_centralsignal.jpg", myimwidth)

        the12chinfostring = bigheadertag("12 channel PA") + headertag(
            "most recent scan: " + breaktag("<span id=\"formtime12\">?</span>\n"))
        the12chroi_imagestr = imagetag("HeadMatrix_roistab.jpg", myimwidth)
        the12chroidrift_imagestr = imagetag("HeadMatrix_roidrift.jpg", myimwidth)
        the12chsnrsfnr_imagestr = imagetag("HeadMatrix_snrsfnr.jpg", myimwidth)
        the12chghost_imagestr = imagetag("HeadMatrix_ghost.jpg", myimwidth)
        the12chobjradius_imagestr = imagetag("HeadMatrix_objradius.jpg", myimwidth)
        the12chobjshape_imagestr = imagetag("HeadMatrix_objshape.jpg", myimwidth)
        the12chweissrdc_imagestr = imagetag("HeadMatrix_weissrdc.jpg", myimwidth)
        the12chcentralsignal_imagestr = imagetag("HeadMatrix_centralsignal.jpg", myimwidth)

        the32chinfostring = bigheadertag("32 channel PA") + headertag(
            "most recent scan: " + breaktag("<span id=\"formtime32\">?</span>\n"))
        the32chroi_imagestr = imagetag("32Ch_Head_roistab.jpg", myimwidth)
        the32chroidrift_imagestr = imagetag("32Ch_Head_roidrift.jpg", myimwidth)
        the32chsnrsfnr_imagestr = imagetag("32Ch_Head_snrsfnr.jpg", myimwidth)
        the32chghost_imagestr = imagetag("32Ch_Head_ghost.jpg", myimwidth)
        the32chobjradius_imagestr = imagetag("32Ch_Head_objradius.jpg", myimwidth)
        the32chobjshape_imagestr = imagetag("32Ch_Head_objshape.jpg", myimwidth)
        the32chweissrdc_imagestr = imagetag("32Ch_Head_weissrdc.jpg", myimwidth)
        the32chcentralsignal_imagestr = imagetag("32Ch_Head_centralsignal.jpg", myimwidth)

        row0str = tablerowtag(
            tableentrytag("") + tableentrytag(thecpinfostring) + tableentrytag(the12chinfostring) + tableentrytag(
                the32chinfostring))
        row1str = tablerowtag(tableentrytag(ROIstab_headerstr) + tableentrytag(thecproi_imagestr) + tableentrytag(
            the12chroi_imagestr) + tableentrytag(the32chroi_imagestr))
        row2str = tablerowtag(tableentrytag(ROIdrift_headerstr) + tableentrytag(thecproidrift_imagestr) + tableentrytag(
            the12chroidrift_imagestr) + tableentrytag(the32chroidrift_imagestr))
        row3str = tablerowtag(tableentrytag(snrsfnr_headerstr) + tableentrytag(thecpsnrsfnr_imagestr) + tableentrytag(
            the12chsnrsfnr_imagestr) + tableentrytag(the32chsnrsfnr_imagestr))
        row4str = tablerowtag(tableentrytag(ghostamp_headerstr) + tableentrytag(thecpghost_imagestr) + tableentrytag(
            the12chghost_imagestr) + tableentrytag(the32chghost_imagestr))
        row5str = tablerowtag(tableentrytag(objradius_headerstr) + tableentrytag(thecpobjradius_imagestr) + tableentrytag(
            the12chobjradius_imagestr) + tableentrytag(the32chobjradius_imagestr))
        row6str = tablerowtag(tableentrytag(objshape_headerstr) + tableentrytag(thecpobjshape_imagestr) + tableentrytag(
            the12chobjshape_imagestr) + tableentrytag(the32chobjshape_imagestr))
        row7str = tablerowtag(tableentrytag(weissrdc_headerstr) + tableentrytag(thecpweissrdc_imagestr) + tableentrytag(
            the12chweissrdc_imagestr) + tableentrytag(the32chweissrdc_imagestr))
        row8str = tablerowtag(
            tableentrytag(centralsignal_headerstr) + tableentrytag(thecpcentralsignal_imagestr) + tableentrytag(
                the12chcentralsignal_imagestr) + tableentrytag(the32chcentralsignal_imagestr))
        fp.writelines(tablepropstag(
            imagehdrstr + row0str + row1str + row2str + row3str + row4str + row5str + row6str + row7str + row8str,
            int(3.5 * myimwidth), "left"))

        # put in a key to help interpret notations on the individual entries
        specs = sf.getlimits('TxRx_Head')
        row0 = headertag("Key:")
        row1 = tableentrytag(yellowtag("Warning, ") + redtag("Out of spec")) + tableentrytag("")
        keyrows = ""
        for theentry in specs:
            entrydata = specs[theentry]
            # TODO this is a hack to get it to run but is not what i want. dmd.
            try:
                if entrydata[2][1] == 1:
                    keyrows = keyrows + tablerowtag(tableentrytag(entrydata[4]) + tableentrytag(entrydata[3]))
            except KeyError:
                pass
        fp.writelines(tablepropstag(bigtag(row0 + row1 + keyrows), int(1.0 * myimwidth), "left"))

        # copy and link to individual reports
        fp.writelines(headertag("Links to individual reports"))

        fp.writelines(tablepropsstarttag(int(3.5 * myimwidth), "left"))
        fp.writelines(tablerowstarttag())
        fp.writelines(tableentryopttag('', widthopt(14.44444)))

        for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
            specs = sf.getlimits(targetcoil)
            coildirlist = bigtag(smallheadertag(targetcoil))
            for i in range(filenumber_TARGET - 1, -1, -1):
                if datadict[i]['Coil'] == targetcoil:
                    # copy the individual scan data if necessary
                    dat_procresults = pjoin(datadirectory, datadict[i]['datadir'])
                    out_procresults = pjoin(outputdirectory, whichscan, datadict[i]['datadir'])
                    if os.path.exists(out_procresults):
                        copypreamble = out_procresults + " exists..."
                        desttime = os.path.getmtime(out_procresults)
                        sourcetime = os.path.getmtime(pjoin(datadirectory, datadict[i]['datadir']))
                        if sourcetime >= desttime:
                            logging.debug(copypreamble + "and is modified - copying " + dat_procresults)
                            logging.debug('time difference={}'.format(desttime - sourcetime))
                            shutil.rmtree(out_procresults)
                            shutil.copytree(dat_procresults, out_procresults)
                        else:
                            logging.debug(copypreamble + "and is current - not copying")
                    else:
                        logging.debug(out_procresults + " does not already exist... copying")
                        shutil.copytree(dat_procresults, out_procresults)

                    # check the data quality
                    thedataquality = stability_eval(specs, datadict[i])
                    themarker = ""
                    if thedataquality != {}:
                        flag = 0
                        for theentry in thedataquality:
                            if thedataquality[theentry]['critical']:
                                if thedataquality[theentry]['quality'] > flag:
                                    flag = thedataquality[theentry]['quality']
                                if thedataquality[theentry]['quality'] > 0:
                                    themarker = themarker + sf.qualitytag(thedataquality[theentry]['symbol'], flag)
                                else:
                                    themarker += " "

                    # generate the output line
                    coildirlist = coildirlist + "<p><a href=" + datadict[i]['datadir'] + "/output.html>" + datadict[i][
                        'Date'] + " " + datadict[i]['Time'] + "</a> " + themarker + "</p>\n"
            fp.writelines(tableentryopttag(coildirlist, widthopt(28.88888) + valignopt("baseline")))
        fp.writelines(tablerowendtag())
        fp.writelines(tablepropsendtag())

        fp.writelines("</body>\n")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create the stability summary.')
    parser.add_argument('--nonbirn', action='store_true', help='process as non-BIRN data')
    parser.add_argument('datadirectory')
    parser.add_argument('outputdirectory')
    parser.add_argument('whichscan')
    args = parser.parse_args()

    stabilitysummary(args.datadirectory, args.outputdirectory, args.whichscan, not args.nonbirn)
