#!/usr/bin/env python

"""processstability.py: convert dicoms, run stabilitycalc and stabilitysummary"""

from os.path import join as pjoin
import os
from glob import glob
import re
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

from stabilityfuncs import stabilityparms, dict_from_tsvfile
from dicom2nifti import dicom2nifti
from stabilitycalc import stabilitycalc
from stabilitysummary import stabilitysummary

rawdicomdir = stabilityparms('rawdicomdir')
sorteddicomdir = stabilityparms('sorteddicomdir')
processedscandir = stabilityparms('processedscandir')
outputdir = stabilityparms('outputdir')

scantypes = {'epi_bold': 'ep2d_bold_normalize_*',
             'epi_multiflip10': 'ep2d_multiflip_10_*',
             'epi_multiflip77': 'ep2d_multiflip_77_*',
             'epi_pace': 'ep2d_pace_normalize_*',
             'epi_singleslice': 'epi_stability_single*sPlot*'}


def processstability(session, dicom=True, summ=True, recalc=True):
    """main"""

    def name(pattern):
        n = [x for x in glob(pjoin(sorteddicomdir, session, pattern)) if 'adaptive' not in x]
        if not n:
            return ''
        else:
            return n[0]

    for scantype, pattern in scantypes.items():
        scan = name(pattern)
        if scan:
            logging.debug('processstability: {}, {}'.format(scantype, scan))
            if recalc:
                _dicom_to_stabilitycalc(dicomseries=scan, niftiname=scantype, dicom=dicom)

            if scantype == 'epi_pace':
                analsum = dict_from_tsvfile(pjoin(processedscandir, session, scantype, 'procresults', 'analysissummary.txt'))

                for unc in glob(pjoin(sorteddicomdir, session, 'uncombined*')):
                    bunc = os.path.basename(unc)
                    m = re.search(r'uncombined(H\d+)_', bunc)
                    try:
                        ele = m.group(1)
                    except AttributeError:
                        ele = bunc

                    ele = 'element_' + ele

                    logging.debug('processstability: PACE: doing scan={}, unc={}, ele={}'.format(scan, bunc, ele))
                    if recalc:
                        _dicom_to_stabilitycalc(dicomseries=bunc, niftiname=ele, dicom=dicom, starttime=0,
                                               initxcenter=analsum['center_of_mass_x'],
                                               initycenter=analsum['center_of_mass_y'],
                                               initzcenter=analsum['center_of_mass_z'])

            if summ:
                stabilitysummary(processedscandir, pjoin(outputdir, '3T', 'birn'), scantype, TargetisBIRNphantom=True)

                if scantype != 'epi_pace':
                    stabilitysummary(processedscandir, pjoin(outputdir, '3T', 'nonbirn'), scantype, TargetisBIRNphantom=False)

    logging.info('processstability done.')

def _dicom_to_stabilitycalc(dicomseries, niftiname, dicom=True, starttime=10, initxcenter=None, initycenter=None, initzcenter=None):

    if dicom:
        dicom2nifti(pjoin(sorteddicomdir, session, dicomseries),
                    pjoin(processedscandir, session), niftiname)

    stabilitycalc(pjoin(processedscandir, session, niftiname), niftiname + '.nii.gz', starttime, initxcenter, initycenter, initzcenter)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Convert DICOMs, run stability{calc,summary}.')
    parser.add_argument('session', help='The session to work on, relative to SORTEDDICOMDIR or PROCESSEDSCANDIR.')
    parser.add_argument('--nodicom', help='Do not reprocess DICOM files; assume they are already done.', action='store_true')
    parser.add_argument('--norecalc', help='Do not reprocess stability; assume it is already done. (Just do summary.)', action='store_true')
    parser.add_argument('--nosumm', help='Do not reprocess summaries.', action='store_true')
    args = parser.parse_args()

    session = args.session

    processstability(session, not args.nodicom, not args.nosumm, not args.norecalc)
