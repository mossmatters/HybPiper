#!/usr/bin/env python

import setuptools

hybpiper_scripts = ['hybpiper']
hybpiper_description = 'Recovery of target gene sequences from bait-capture data'
hybpiper_url = 'https://github.com/chrisjackson-pellicle/HybPiper.git'
hybpiper_entry_points = {'console_scripts': ['reads_first = hybpiper.reads_first:main',
                                             'retrieve_sequences = hybpiper.retrieve_sequences:main',
                                             'get_seq_lengths = hybpiper.retrieve_sequences:main',
                                             'paralog_retriever = hybpiper.paralog_retriever:main']}

setuptools.setup(name='hybpiper',
                 version='1.4',
                 packages=setuptools.find_packages(),
                 scripts=['hybpiper/gene_recovery_heatmap_ggplot.R',
                          'hybpiper/reads_first.py',
                          'hybpiper/distribute_reads_to_targets_bwa.py',
                          'hybpiper/distribute_reads_to_targets.py',
                          'hybpiper/distribute_targets.py',
                          'hybpiper/spades_runner.py',
                          'hybpiper/exonerate_hits.py'],
                 author='Chris Jackson, Matt Johnson',
                 author_email='chris.jackson@rbg.vic.gov.au',
                 description=hybpiper_description,
                 keywords='target-capture phylogeny',
                 url=hybpiper_url,
                 entry_points=hybpiper_entry_points)
