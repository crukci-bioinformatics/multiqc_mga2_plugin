#!/usr/bin/env python
"""
Multi Genome Alignment 2 plugin for MultiQC.

For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '2.0'

setup(
    name = 'multiqc_mga2_plugin',
    version = version,
    author = 'Richard Bowers',
    author_email = 'richard.bowers@cruk.cam.ac.uk',
    description = "Multi Genome Alignment 2 MultiQC plugin",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/crukci-bioinformatics/multiqc_mga2_plugin',
    download_url = 'https://github.com/crukci-bioinformatics/multiqc_mga2_plugin/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    python_requires = '>=3.6',
    install_requires = [
        'multiqc>=1.9', 'natsort'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'mga2 = multiqc_mga2.modules.mga2:MultiqcModule',
        ],
        'multiqc.cli_options.v1': [
            'mga2_sequencing_run = multiqc_mga2.cli:mga2_sequencing_run',
            'mga2_title = multiqc_mga2.cli:mga2_title'
        ],
        'multiqc.hooks.v1': [
            'execution_start = multiqc_mga2.custom_code:mga2_plugin_execution_start'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
