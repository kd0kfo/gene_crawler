#!/usr/bin/env python

from distutils.core import setup, Command


class Tester(Command):
    user_options = []

    def initialize_options(self):
        import os
        self._dir = os.getcwd()

    def finalize_options(self):
        pass

    def run(self):
        from test import unit_test
        unit_test.run()

the_scripts = ['scripts/gene_extractor', 'scripts/regex_search',
               'scripts/gene_indexer']

setup(name='gene_crawler',
       version='1.0',
       url='http://code.davecoss.com',
       license='GPL v3',
       description=('Utility functions for searching and manipulating'
                    ' gene data'),
       author_email='David.Coss@stjude.org',
       packages=['gene_crawler'],
       scripts=the_scripts,
       cmdclass={'test': Tester}
       )
