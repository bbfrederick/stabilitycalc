#!/usr/bin/env python

import os
from os.path import join as pjoin
import shutil
import time
import logging
import matplotlib
matplotlib.use('Agg', warn=False)
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style("dark")
import pandas as pd
from collections import OrderedDict
from datetime import datetime
import csv
import itertools
import stabilityfuncs as sf
from mako.lookup import TemplateLookup
makolookup = TemplateLookup(directories=['./tpl'])


# noinspection PyPep8Naming
def stabilitysummary(datadirectory, outputdirectory, whichscan, TargetisBIRNphantom):
    logging.debug('stabilitysummary: running as stabilitysummary {} {} {} {}'.format(datadirectory,
                                                                                     outputdirectory,
                                                                                     whichscan,
                                                                                     '' if TargetisBIRNphantom else '--nonbirn'))

    # initialize the outut directory if need be
    if not os.path.exists(pjoin(outputdirectory, whichscan)):
        os.makedirs(pjoin(outputdirectory, whichscan))

    phantomtype = 'BIRN' if TargetisBIRNphantom else 'NONBIRN'

    # scan the data directory for stability scans and populate the dict
    stabilitydirs = os.listdir(datadirectory)
    stabilitydirs = sorted([filename for filename in stabilitydirs if filename.startswith("stability_")])

    datadict = {}
    filenumber_TARGET = 0
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
                    filenumber_TARGET += 1
            except IOError:
                pass
        except KeyError:
            pass

    plot_epoch = pd.Timestamp(sf.stabilityparms('epoch', 'plots'))
    plot_timelims = (plot_epoch, pd.Timestamp(datetime.now()))
    df = pd.DataFrame.from_dict(datadict, orient='index', dtype=float)
    df.DateTime = pd.to_datetime(df.DateTime)

    # normalize the min and max error values by the means
    df.odd_ghost_min -= df.odd_ghost_mean
    df.odd_ghost_max -= df.odd_ghost_mean
    df.even_ghost_min -= df.even_ghost_mean
    df.even_ghost_max -= df.even_ghost_mean

    # read the plot config file
    plotconfig = {}
    for row in csv.DictReader(open('config/plots.csv')):
        key = row.pop('plot')
        row['variables'] = row['variables'].split(';')
        row['legends'] = row['legends'].split(';')
        plotconfig[key] = row

    mostrecenttimes = {}

    for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
        mostrecenttimes[targetcoil] = df[df.Coil == targetcoil].DateTime.max()

        dftc = df[df.Coil == targetcoil].groupby('DateTime')

        # plot ghost. this one is a little special because of the error bars, so we do it manually
        even_errors = np.array(df.ix[df.Coil == targetcoil, ['even_ghost_min', 'even_ghost_max']]).T
        odd_errors = np.array(df.ix[df.Coil == targetcoil, ['odd_ghost_min', 'odd_ghost_max']]).T
        plt.hold(True)
        dftc['even_ghost_mean'].mean().plot(yerr=even_errors, marker='*', linestyle='None', label='Evens')
        dftc['odd_ghost_mean'].mean().plot(yerr=odd_errors, marker='.', linestyle='None', label='Odds')
        plt.xlim(plot_timelims)
        plt.ylim(0, 15)
        plt.title('Ghost percentage')
        plt.xlabel('Date')
        plt.ylabel('Ghost amplitude (%)')
        plt.legend()
        plt.savefig(pjoin(outputdirectory, whichscan, '{}_ghost.png'.format(targetcoil)), format='png')
        plt.hold(False)
        plt.close()

        # now let's do the rest from the config file.
        for plotfilename, config in plotconfig.items():
            plt.hold(True)
            marker = itertools.cycle(('*', '.', 'o', '+', 'h'))
            for i, variable in enumerate(config['variables']):
                dftc[variable].mean().plot(marker=marker.next(), linestyle='None', label=config['legends'][i])
            plt.xlim(plot_timelims)
            plt.ylim(float(config['ymin']), float(config['ymax']))
            plt.xlabel('Date')
            plt.ylabel(config['ylabel'])
            plt.title(config['title'])
            plt.legend()
            plt.hold(False)
            plt.savefig(pjoin(outputdirectory, whichscan, '{}_{}.png'.format(targetcoil, plotfilename)), format='png')
            plt.close()

    for targetcoil in ['TxRx_Head', 'HeadMatrix', '32Ch_Head']:
        outscandir = pjoin(outputdirectory, whichscan)

        # TODO there's a simpler way to do this.
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
    args32 = str(mostrecenttimes.get('32Ch_Head', '1970')).replace(' ', 'T')
    args12 = str(mostrecenttimes.get('HeadMatrix', '1970')).replace(' ', 'T')
    argscp = str(mostrecenttimes.get('TxRx_Head', '1970')).replace(' ', 'T')

    tpl = makolookup.get_template('stabilityreport.html')
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
                                themarker = themarker + sf.qualitytag(spec['flag'], flag)
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
