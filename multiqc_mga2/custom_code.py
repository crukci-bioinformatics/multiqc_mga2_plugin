#!/usr/bin/env python
""" MultiQC MGA2 plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.mga2_plugin_version = get_distribution("multiqc_mga2_plugin").version


# Add default config options for the things that are used in MultiQC_NGI
def mga2_plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    log.debug("Running Multi Genome Alignment 2 MultiQC Plugin v{}".format(config.mga2_plugin_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    # Add to the search patterns used by modules
    if 'mga2' not in config.sp:
        config.update_dict( config.sp, { 'mga2': { 'fn': 'genome_alignment_summary.csv' } } )
