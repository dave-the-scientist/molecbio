#!/usr/bin/python
from distutils.core import setup, Extension
import os.path
from __init__ import __version__

setup(
    name='molecbio', version=__version__,
    author='Dave Curran', author_email='curran.dave.m@gmail.com',
    url='https://github.com/dave-the-scientist/molecbio',
    description='',
    license='LICENSE',
    package_dir={'molecbio':''},
    packages=['molecbio', 'molecbio.blosum', 'molecbio.rosetta', 'molecbio.aligners'],
    ext_modules = [Extension('molecbio.aligners.nwmodule', [os.path.join('aligners', 'nwmodule.c')],
                             extra_compile_args=['-std=c99'])]
    )

