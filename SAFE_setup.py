#!/usr/bin/env python

from distutils.core import setup

setup(name='stabilitycalc',
      description='evaluate fMRI scanner stability',
      long_description="""
Command-line tools to calculate numerous fMRI scanner stability metrics, based
on the FBIRN quality assurance test protocal. Any 4D volumetric timeseries
image in NIfTI format is support input. Output is a rich HTML report.
""",
      url='https://github.com/bbfrederick/stabilitycalc',
      license='BSD',
      author='Blaise Frederick',
      author_email='blaise.frederick@gmail.com',
      version='0.1',
      py_modules=['stabilityfuncs',
                   'htmltagutils', ],
      scripts=['getdicominfo',
               'stabilitycalc'],
      )
