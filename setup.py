#!/usr/bin/env python

from distutils.core import setup

setup(name='stabilitycalc',
      description='DESCRIPTION',
      long_description='LONG_DESCRIPTION',
      url='URL',
      download_url='DOWNLOAD_URL',
      license='BSD',
      author='AUTHOR',
      author_email='AUTHOR_EMAIL',
      version='0.1',
      py_modules   = ['stabilityfuncs',
                      'htmltagutils',],
      scripts      = ['getdicominfo',
                      'stabilitycalc'],
     )
