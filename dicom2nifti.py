#!/usr/bin/env python

import sys
import os
from os.path import join as pjoin
import shutil
from glob import glob
import logging
from nipype.interfaces.dcm2nii import Dcm2nii
import studyinfo as si

def dicom2nifti(dicomdir, niftidir, niftiname, studyinfo=True):
    # clean up anything we did earlier
    shutil.rmtree(pjoin(niftidir, niftiname), ignore_errors=True)
    os.makedirs(pjoin(niftidir, niftiname))

    # do the conversion
    c = Dcm2nii()
    c.inputs.source_dir = dicomdir
    c.inputs.output_dir = pjoin(niftidir, niftiname)
    c.inputs.gzip_output = True
    c.inputs.date_in_filename = False
    c.run()

    logging.info('dicom2nifti: wrote {}/{}'.format(niftidir, niftiname))

    # but, we're weird and we want a standardized name
    dest = pjoin(niftidir, niftiname, niftiname + '.nii.gz')
    os.rename(glob(pjoin(niftidir, niftiname, '*nii.gz'))[0], dest)

    if studyinfo:
        si.embed_studyinfo(dest, si.studyinfo_from_dicom(glob(pjoin(dicomdir, 'IM*000[1-9].dcm'))[0]))

    return dest

if __name__ == '__main__':
    import argparse
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    parser = argparse.ArgumentParser(description='Convert dicom to nifti.')
    parser.add_argument('dicomdir', help='The dicomdir to convert.')
    parser.add_argument('niftidir', help='Destination NIFTI directory.')
    parser.add_argument('niftiname', help='Name of NIFTI file to write.')

    args = parser.parse_args()

    dicom2nifti(args.dicomdir, args.niftidir, args.niftiname)
    
    try:
        info = si.studyinfo_from_dicom(glob(pjoin(args.dicomdir, 'IM*000[1-9].dcm'))[0])
        si.studyinfo_write(pjoin(args.niftidir, args.niftiname, 'studyinfo'), info)
    except IndexError:
        logging.warning('dicom2nifti: no IM*0001.dcm found, writing empty studyinfo file')
        open(pjoin(args.niftidir, args.niftiname, 'studyinfo'), 'a')


