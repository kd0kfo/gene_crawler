#!/usr/bin/env python

from distutils.core import setup, Extension

the_scripts = ['scripts/gene_extractor','scripts/regex_search','scripts/gene_indexer']

setup (name ='gene_crawler',
       version = '0.3',
       url = 'http://code.davecoss.com',
       license = 'GPL v3',
       description = 'Utility functions for searching and manipulating gene data',
       author_email='David.Coss@stjude.org',
       packages = ['gene_crawler'],
       scripts = the_scripts,
       )

