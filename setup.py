#!/usr/bin/env python

import setuptools

hybpiper_scripts = ['hybpiper']
hybpiper_description = 'Recovery of target gene sequences from bait-capture data'
hybpiper_url = 'https://github.com/chrisjackson-pellicle/HybPiper.git'
hybpiper_entry_points = {'console_scripts': ['hybpiper = hybpiper.reads_first:main']}

setuptools.setup(name='hybpiper',
                 version='1.4',
                 packages=setuptools.find_packages(),
                 scripts=['hybpiper/gene_recovery_heatmap.py',
                          'hybpiper/hybpiper.py',
                          'hybpiper/distribute_reads_to_targets_bwa.py',
                          'hybpiper/distribute_reads_to_targets.py',
                          'hybpiper/distribute_targets.py',
                          'hybpiper/spades_runner.py',
                          'hybpiper/exonerate_hits.py',
                          'hybpiper/get_seq_lengths.py',
                          'hybpiper/retrieve_sequences.py',
                          'hybpiper/paralog_retriever.py',
                          'hybpiper/hybpiper_stats.py'],
                 author='Chris Jackson, Matt Johnson',
                 author_email='chris.jackson@rbg.vic.gov.au',
                 description=hybpiper_description,
                 keywords='target-capture phylogeny',
                 url=hybpiper_url,
                 entry_points=hybpiper_entry_points)
