#!/usr/bin/env python
"""
MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters.

See the Click documentation for more command line flag types:
https://click.palletsprojects.com/en/7.x/
"""

import click

# Sets config.kwargs['mga2_sequencing_run'] to True if specified (will be False otherwise)
mga2_sequencing_run = click.option('--mga-sequencing-run', 'mga2_sequencing_run', is_flag = True,
    help = "Indicates to the MGA plugin that this run of MultiQC is for a sequencing run folder."
)

# Sets the title for MGA plots.
# Sets config.kwargs['mga2_title']
mga2_title = click.option("--mga-title", "mga2_title", type = str, show_default = False,
    help = "The title to apply to the MGA plots and tables."
)

# Sets which datasets to report on.
# Sets config.kwargs['mga2_datasets']
mga2_datasets = click.option("--mga-datasets", "mga2_datasets", type = str, show_default = False,
    help = "The datasets to limit the MGA2 report to."
)

