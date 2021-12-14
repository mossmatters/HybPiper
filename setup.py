#!/usr/bin/env python

import setuptools

hybpiper_scripts = ['hybpiper']
hybpiper_description = 'Recovery of target gene sequences from bait-capture data'
hybpiper_url = 'https://github.com/chrisjackson-pellicle/HybPiper.git'
hybpiper_entry_points = {'console_scripts': ['reads_first = hybpiper.reads_first:main']}

setuptools.setup(name='HybPiper',
                 version='1.4',
                 packages=setuptools.find_packages(),
                 author='Chris Jackson, Matt Johnson',
                 author_email='chris.jackson@rbg.vic.gov.au',
                 description=hybpiper_description,
                 keywords='target-capture phylogeny',
                 url=hybpiper_url,
                 entry_points=hybpiper_entry_points)
