#!/usr/bin/env python

import os
import dicom
import logging


def studyinfo_from_dicom(dicomfilename):
    """generate a studyinfo dict from a dicom file"""

    def getsiemensmrheader(theplan):
        siemens_csa_header2 = theplan[0x0029, 0x1020].value
        startposition = siemens_csa_header2.find('### ASCCONV BEGIN ###') + len('### ASCCONV BEGIN ###')
        partialheader = siemens_csa_header2[startposition:]
        endposition = partialheader.find('### ASCCONV END ###')
        splitheader = partialheader[:endposition].splitlines()
        d = {}
        for line in splitheader[1:]:
            pair = line.split()
            d[pair[0]] = pair[2]
        return d

    plan = dicom.read_file(dicomfilename)
    siemensheader = getsiemensmrheader(plan)

    info = {'Coil': siemensheader['asCoilSelectMeas[0].asList[0].sCoilElementID.tCoilID'].strip('"'),
            'StudyDate': plan.StudyDate,
            'StudyTime': plan.StudyTime}

    try:
        info['ElementName'] = plan[0x0051, 0x100f].value
    except KeyError:
        info['ElementName'] = 'UNKNOWN'

    return info


def studyinfo(studyinfo_file):
    """read a studyinfo file"""

    if not os.path.exists(studyinfo_file):
        logging.info('studyinfo: studyinfo file {} not found, using empty one'.format(studyinfo_file))
        return {'Coil': '', 'StudyDate': '', 'StudyTime': '', 'ElementName': ''}

    info = {}
    for line in open(studyinfo_file):
        k, v = line.split(': ', 1)
        info[k] = v.strip()

    return info


def studyinfo_write(studyinfo_file, info):
    """write a studyinfo file"""

    with open(studyinfo_file, 'w') as fp:
        for pair in info.items():
            fp.write(': '.join(pair) + '\n')


if __name__ == '__main__':
    import argparse

    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    parser = argparse.ArgumentParser(description='Pull header data from a DICOM file.')
    parser.add_argument('dicomfile', help='The dicom to examine.')
    args = parser.parse_args()

    if not os.path.exists(args.dicomfile):
        logging.critical('studyinfo: dicomfile {} not found'.format(args.dicomfile))
        exit(1)

    for k, v in studyinfo_from_dicom(args.dicomfile).items():
        print('{}: {}'.format(k, v))
