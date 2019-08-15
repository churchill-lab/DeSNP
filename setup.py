#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from glob import glob
with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

import os
on_rtd = os.environ.get('READTHEDOCS', None)

requirements = []
test_requirements = []

if not on_rtd:
    requirements.append('numpy>=1.17')
    requirements.append('pysam>=0.15')

setup(
    name='desnp',
    version='1.0.0',
    description="DeSNP",
    long_description=readme + '\n\n' + history,
    author='Matthew J. Vincent and David Walton, The Jackson Laboratory',
    author_email='matt.vincent@jax.org',
    url='http://churchill-lab.github.io/DeSNP/',
    packages=[
        'desnp',
    ],
    package_dir={'desnp':
                 'desnp'},
    scripts=glob("bin/*"),
    include_package_data=True,
    install_requires=requirements,
    license="MIT",
    zip_safe=False,
    keywords=['desnp', 'snp'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
    ]
)
