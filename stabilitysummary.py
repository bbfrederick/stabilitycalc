#!/usr/bin/env python

import os
from os.path import join as pjoin
import shutil
import time
import logging
import subprocess
from collections import OrderedDict

from mako.lookup import TemplateLookup

makolookup = TemplateLookup(directories=['./tpl'])
from htmltagutils import *
import stabilityfuncs as sf


# noinspection PyPep8Naming
def stabilitysummary(datadirectory, outputdirectory, whichscan, TargetisBIRNphantom):
    logging.debug('stabilitysummary: running as stabilitysummary {} {} {} {}'.format(datadirectory,
                                                                                     outputdirectory,
                                                                                     whichscan,
                                                                                     '' if TargetisBIRNphantom else '--nonbirn'))

    # initialize the outut directory if need be
    if not os.path.exists(pjoin(outputdirectory, whichscan)):
        os.makedirs(pjoin(outputdirectory, whichscan))

    # scan the data directory for stability scans

    stabilitydirs = os.listdir(datadirectory)
    stabilitydirs = sorted([filename for filename in stabilitydirs if filename.startswith("stability_")])

    # pull all the data files into a dictionary array

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
                datadict[filenumber_TARGET].update(sf.dict_from_tsvfile(
                    pjoin(datadirectory, datadict[filenumber_TARGET]['datadir'], 'analysissummary.txt')))
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
    phantomtype = 'BIRN' if TargetisBIRNphantom else 'NONBIRN'
    logging.debug("{} CP coil runs ({} phantom)".format(num_cp_TARGET, phantomtype))
    logging.debug("{} 12 channel coil runs ({} phantom)".format(num_12_TARGET, phantomtype))
    logging.debug("{} 32 channel coil runs ({} phantom)".format(num_32_TARGET, phantomtype))

    # sort the data up by coil and write to files

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

    # generate plot control files to graph all interesting stability parameters
    # central and peripheral SNR and SFNR

    outplotnames = ('plotcmds_centralsignal',
                    'plotcmds_ghost',
                    'plotcmds_objradius',
                    'plotcmds_objshape',
                    'plotcmds_roidrift',
                    'plotcmds_roistab',
                    'plotcmds_snrsfnr',
                    'plotcmds_weissrdc')
    for outplotname in outplotnames:
        tpl = makolookup.get_template(outplotname)
        outplotfile = pjoin(outputdirectory, whichscan, outplotname)
        with open(outplotfile, 'w') as fp:
            fp.write(tpl.render(outputdirectory=outputdirectory,
                                whichscan=whichscan,
                                wlp='points' if TargetisBIRNphantom else 'linespoints'))

    for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
        outscandir = pjoin(outputdirectory, whichscan)
        shutil.copyfile(pjoin(outscandir, targetcoil + '_vals.txt'), pjoin(outscandir, 'graphtemp'))
        for plottype in ['snrsfnr', 'roistab', 'roidrift', 'ghost', 'objradius', 'objshape', 'weissrdc',
                         'centralsignal']:
            subprocess.Popen(['gnuplot', pjoin(outscandir, 'plotcmds_' + plottype)],
                             stdout=open(pjoin(outscandir, '{}_{}.jpg'.format(targetcoil, plottype)), 'wb'))

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

    # generate a report file

    thisdate = time.strftime("%m/%d/%Y %H:%M:%S", time.localtime())
    date32 = mostrecenttimes.get('32Ch_Head', '19700101T000000')
    date12 = mostrecenttimes.get('HeadMatrix', '19700101T000000')
    datecp = mostrecenttimes.get('TxRx_Head', '19700101T000000')
    args32 = ','.join((date32[0:4], str(int(date32[4:6]) - 1), date32[6:8], date32[9:11], date32[11:13], date32[13:15]))
    args12 = ','.join((date12[0:4], str(int(date12[4:6]) - 1), date12[6:8], date12[9:11], date12[11:13], date12[13:15]))
    argscp = ','.join((datecp[0:4], str(int(datecp[4:6]) - 1), datecp[6:8], datecp[9:11], datecp[11:13], datecp[13:15]))

    tpl = makolookup.get_template('stabilityreport.html')
    tplindividual = makolookup.get_template('stabilityreport_individual.html')
    with open(pjoin(outputdirectory, whichscan, 'stabilityreport.html'), 'w') as fp:
        coiltemplatedata = OrderedDict()
        for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
            coiltemplatedata[targetcoil] = []
            specs = sf.getlimits(targetcoil)
            for i in range(filenumber_TARGET - 1, -1, -1):
                if datadict[i]['Coil'] == targetcoil:
                    themarker = ""
                    flag = 0
                    for specid, spec in specs.items():
                        if spec['critical']:
                            if sf.limitcheck(datadict[i][specid], spec) > flag:
                                flag = sf.limitcheck(datadict[i][specid], spec)
                            if sf.limitcheck(datadict[i][specid], spec) > 0:
                                themarker = themarker + qualitytag(spec['flag'], flag)
                            else:
                                themarker += " "

                    # populate the output db
                    coiltemplatedata[targetcoil].append({
                        'path': datadict[i]['datadir'],
                        'datetime': datadict[i]['Date'] + ' ' + datadict[i]['Time'],
                        'marker': themarker})

        fp.write(tpl.render(**locals()))

if __name__ == '__main__':
    import argparse

parser = argparse.ArgumentParser(description='Create the stability summary.')
parser.add_argument('--nonbirn', action='store_true', help='process as non-BIRN data')
parser.add_argument('datadirectory')
parser.add_argument('outputdirectory')
parser.add_argument('whichscan')
args = parser.parse_args()

stabilitysummary(args.datadirectory, args.outputdirectory, args.whichscan, not args.nonbirn)
