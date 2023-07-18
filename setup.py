#!/usr/bin/env python

import setuptools

# Read version number into a dictionary
version = {}
with open('hybpiper/version.py') as fp:
    exec(fp.read(), version)


hybpiper_scripts = ['hybpiper']
hybpiper_description = 'Recovery of target gene sequences from bait-capture data'
hybpiper_url = 'https://github.com/mossmatters/HybPiper'
hybpiper_entry_points = {'console_scripts': ['hybpiper = hybpiper.assemble:main']}

setuptools.setup(name='hybpiper',
                 version=version['__version__'],
                 packages=setuptools.find_packages(),
                 author='Chris Jackson, Matt Johnson',
                 author_email='chris.jackson@rbg.vic.gov.au',
                 description=hybpiper_description,
                 keywords='target-capture phylogeny',
                 url=hybpiper_url,
                 entry_points=hybpiper_entry_points)
