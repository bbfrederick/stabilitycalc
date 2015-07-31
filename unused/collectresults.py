#!/usr/bin/env python

from stabilityfuncs import stabilityparms, dict_from_tsvfile
import os
from os.path import join as pjoin
from glob import glob
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def collectresults(files, output):
    keys = dict_from_tsvfile(files[0]).keys()

    tsvfile = open(output, 'w')
    # write the header line
    tsvfile.write('\t'.join(keys))
    tsvfile.write('\n')
    for fn in files:
        d = dict_from_tsvfile(fn)
        if d.keys() != keys:
            logging.critical('{} contained different keys from the first file. Aborting write, removing.'.format(file))
            tsvfile.close()
            os.remove(output)
            exit(1)
        for k in keys:
            tsvfile.write(d[k])
            tsvfile.write('\t')
        tsvfile.write('\n')

if __name__ == '__main__':
    processedscandir = stabilityparms('processedscandir')
    analysissums = glob(pjoin(processedscandir, 'stability*', 'epi_pace', 'procresults', 'analysissummary.txt'))

    collectresults(analysissums, pjoin(processedscandir, 'stabilitytable.txt'))
