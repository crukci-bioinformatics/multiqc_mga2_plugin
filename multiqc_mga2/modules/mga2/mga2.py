#!/usr/bin/env python

""" MultiQC Multi Genome Alignment 2 module """

import csv
import locale
import logging
import operator
import os
import re

from collections import OrderedDict
from functools import cmp_to_key
from natsort import natsorted
from types import SimpleNamespace

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

from .colour import Colour
from .mga_structures import MGAData, MGADataset, MGAAssignment, MGADatasetSummary

log = config.logger

# See https://stackoverflow.com/a/16279578
bar_colours = SimpleNamespace(**{
    'reference': Colour.fromBytes(0, 255, 0),
    'control': Colour.fromBytes(255, 200, 0),
    'contaminant': Colour.fromBytes(255, 0, 0),
    'unmapped': Colour.fromBytes(224, 224, 224),
    'adapter': Colour.fromBytes(255, 102, 255),
    'noadapter': Colour.fromBytes(244, 244, 244)
})

max_alpha = 1.0
min_alpha = 0.2 # Is 0.1 in the Java version.
max_error = 0.015
min_error = 0.0025

assigned_fraction_threshold = 0.005
aligned_fraction_threshold = 0.01
error_rate_threshold = 0.0125
adapter_threshold_multiplier = 0.005


# Based on https://github.com/MultiQC/example-plugin

class MultiqcModule(BaseMultiqcModule):
    '''
    A MultiQC module for MGA2.
    '''

    def __init__(self):
        '''
        Constructor. Causes the processing to run too.
        '''

        super(MultiqcModule, self).__init__(
            name = 'Multi Genome Alignment',
            target = "MGA",
            anchor = 'mga',
            href = "https://github.com/crukci-bioinformatics/MGA2",
            info = """(multi-genome alignment) is a quality control tool for high-throughput sequence data
                   developed by the Bioinformatics Core at the Cancer Research UK Cambridge Institute."""
        )

        # Add to self.css to be included in template
        self.css = { 'assets/css/multiqc_mga2.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_mga2.css') }

        mga_data = self.load_mga_files()

        title = config.kwargs.get("mga2_title")
        if title is not None: mga_data.title = title

        mga_data.from_sequencing = config.kwargs.get('mga2_sequencing_run', False)

        datasets_arg = config.kwargs.get("mga2_datasets")
        if datasets_arg is not None:
            datasets_of_interest = map(str.strip, re.split(r"[ ,;:|/]+", datasets_arg))
            datasets_of_interest = [ds for ds in datasets_of_interest if len(ds) > 0]
            mga_data.filter_datasets(datasets_of_interest)

        self.create_mga_reports(mga_data)


    def load_mga_files(self) -> MGAData:
        """
        Find and load the content of MGA2 mga_alignment_summary.csv files.

        :return A MGAData structure populated with the information from the alignment
        summary file and the MGA summary file..
        :rtype MGAData
        """

        mga_data = MGAData()

        for mgafile in self.find_log_files('mga2', filecontents=False, filehandles=False):
            self._read_mga2_csv_file(mga_data, mgafile)

        if len(mga_data.datasets) == 0:
            raise UserWarning

        log.info("Found {} MGA2 datasets".format(len(mga_data.datasets)))

        mga_data.calculate_max_sequences()

        return mga_data


    def _read_mga2_csv_file(self, mga_data: MGAData, mgafile: dict):
        """
        Read a MGA2 alignment summary file.

        The file given as "mgafile" is the alignment summary file. We expect to find a file
        for the overall MGA summary next to it, such if the original is prefix.mga_alignment_summary.csv
        the summary file is prefix.mga_summary.csv

        :param MGAData mga_data: The MGAData structure to read into.
        :param dict mgafile: The CSV file to load.
        """

        mgafile_path = os.path.join(mgafile["root"], mgafile["fn"])
        summaryfile_name = mgafile['fn'].replace("mga_alignment_summary", "mga_summary")
        summaryfile_path = os.path.join(mgafile["root"], summaryfile_name)

        if not os.path.exists(summaryfile_path):
            log.warn("Cannot find {} corresponding to {} in {}".format(summaryfile_name, mgafile['fn'], mgafile['root']))
            return

        # See https://docs.python.org/3/library/csv.html#csv.DictReader

        with open(mgafile_path, newline='') as fh:
            log.debug(f"Found file {mgafile_path}")
            reader = csv.DictReader(fh, dialect = 'unix')
            for assignment in reader:
                mga_data.add_assignment_from_csv(assignment)

        with open(summaryfile_path, newline='') as fh:
            log.debug(f"Found file {summaryfile_path}")
            reader = csv.DictReader(fh, dialect = 'unix')
            for summary in reader:
                mga_data.set_summary_from_csv(summary)


    def create_mga_reports(self, mga_data: MGAData):
        '''
        Build the sections for the overall report

        :param MGAData mga_data: The data loaded from the MGA files.
        '''

        # The plot.

        plot_data, plot_categories = self._plot_data(mga_data)

        self.add_section(
            name = 'Summary',
            anchor = 'mga_summary',
            # description = self._plot_description(mga_data),
            content = self._plot_key(),
            helptext = self._plot_help(mga_data),
            plot = bargraph.plot(plot_data, plot_categories, self._plot_config(mga_data))
        )

        # The summary table (one per dataset).

        for dataset_id, mga_dataset in mga_data.datasets.items():

            self.add_section(
                name = f'Lane {dataset_id} Statistics' if mga_data.from_sequencing else f'Dataset "{dataset_id}" Statistics',
                anchor = f"mga_stats_{dataset_id}",
                plot = table.plot(self._main_table_data(mga_dataset), self._main_table_headers(), self._main_table_config(dataset_id)),
                description = f"""
                    <table>
                        <tr><td>Sequences:</td><td>{mga_dataset.summary.sequences:,}</td></tr>
                        <tr><td>Sampled:</td><td>{mga_dataset.summary.sampled:,}</td></tr>
                    </table>
                    <br/>
                """
            )


    def _plot_data(self, mga_data: MGAData) -> tuple:
        '''
        Assemble the plot information.

        :param MGAData mga_data: The data loaded from the MGA files.

        :return A tuple of bargraph plot data and categories. Both values are dictionaries keyed
        by dataset id.
        :rtype tuple
        '''

        bar_data = OrderedDict()
        categories = OrderedDict()

        for dataset_id, mga_dataset in mga_data.datasets.items():
            dataset_summary = mga_dataset.summary

            sequence_count = dataset_summary.sequences
            sampled_count = dataset_summary.sampled

            sampled_to_sequenced = 0 if sequence_count == 0 else float(sequence_count) / float(sampled_count)

            # Sort assignments first.
            sorted_assignments = sorted(mga_dataset.assignments.values())

            dataset_bar_data = OrderedDict()
            dataset_categories = dict()

            for assignment in sorted_assignments:
                if self._accept_genome(assignment):
                    category_id = f"{dataset_id}.{assignment.genome}"

                    colour = bar_colours.contaminant
                    if assignment.control:
                        colour = bar_colours.control
                    elif assignment.expected:
                        colour = bar_colours.reference

                    # log.debug("Dataset {} - {}".format(dataset_id, reference_genome_id))
                    # log.debug("max_alpha - (max_alpha - min_alpha) * (assigned_error - min_error) / (max_error - min_error)")
                    # log.debug("{} - ({} - {}) * ({} - {}) / ({} - {})".format(max_alpha, max_alpha, min_alpha, assigned_error, min_error, max_error, min_error))
                    # log.debug("{} - {} * {} / {}".format(max_alpha, max_alpha - min_alpha, assigned_error - min_error, max_error - min_error))

                    error_rate = assignment.assigned_error_rate if assignment.assigned > 0 else assignment.error_rate
                    alpha = max_alpha - (max_alpha - min_alpha) * (error_rate - min_error) / (max_error - min_error)

                    # log.debug("alpha = {}".format(alpha))

                    alpha = max(min_alpha, min(max_alpha, alpha))

                    # log.debug("capped alpha = {}".format(alpha))

                    if assignment.assigned >= 100:
                        log.debug("{}\t{}\t{}\t{}".format(assignment.genome, assignment.assigned, assignment.error_rate * 100.0, alpha))

                    dataset_bar_data[category_id] = int(assignment.assigned * sampled_to_sequenced)

                    dataset_categories[category_id] = {
                        'name': assignment.species,
                        'color': colour.applyAlpha(alpha).toHtml()
                    }

            # See hpw many genomes were not accepted, and if there are any add a section for "other"
            number_of_others, other_aligned_count, other_assigned_count = self._count_others(mga_dataset)

            if number_of_others > 0:
                log.debug("{} others with {} assigned reads".format(number_of_others, other_assigned_count))

                category_id = f"{dataset_id}.other"

                dataset_bar_data[category_id] = int(other_assigned_count * sampled_to_sequenced)

                dataset_categories[category_id] = {
                    'name': f"{number_of_others} others",
                    'color': bar_colours.contaminant.applyAlpha(min_alpha).toHtml()
                }

            bar_data[dataset_id] = dataset_bar_data

            log.debug(f"Unmapped count: {dataset_summary.unmapped} / {sampled_count}")

            # Add the unmapped count.
            category_id = f"{dataset_id}.unmapped"
            dataset_bar_data[category_id] = int(dataset_summary.unmapped * sampled_to_sequenced)
            dataset_categories[category_id] = {
                'name': 'Unmapped',
                'color': bar_colours.unmapped.toHtml()
            }

            # The order of categories matters for the order in the plot, so add them to the
            # overall categories dictionary in the order they are in for the bar data.
            for category_id in dataset_bar_data.keys():
                categories[category_id] = dataset_categories[category_id]

            log.debug(f"Adapter count: {dataset_summary.adapter} / {sampled_count}")

            if dataset_summary.adapter > 0 and dataset_summary.adapter >= sampled_count * adapter_threshold_multiplier:
                dataset_adapter_id = f"{dataset_id}A" if mga_data.from_sequencing else f"{dataset_id} adapter"
                dataset_bar_data = OrderedDict()
                category_id = f"{dataset_id}.adapter"
                dataset_bar_data[category_id] = int(dataset_summary.adapter * sampled_to_sequenced)
                categories[category_id] = {
                    'name': 'Adapter',
                    'color': bar_colours.adapter.toHtml()
                }
                category_id = f"{dataset_id}.noadapter"
                dataset_bar_data[category_id] = int(sequence_count - dataset_summary.adapter * sampled_to_sequenced)
                categories[category_id] = {
                    'name': 'No adapter match',
                    'color': bar_colours.noadapter.toHtml()
                }
                bar_data[dataset_adapter_id] = dataset_bar_data


        log.debug(f"Bar data = {bar_data}")
        log.debug(f"Categories = {categories}")

        return bar_data, categories


    def _plot_config(self, mga_data: MGAData) -> dict:
        '''
        Create the plot configuration for the main MGA plot.

        :param MGAData mga_data: The data loaded from the MGA files.

        :return A dictionary of plot configuration parameters.
        :rtype dict
        '''
        return {
            'id': "mga_plot_{}".format(mga_data.title.replace(' ', '_')),
            'title': f"Multi Genome Alignment: {mga_data.title}",
            'cpswitch_counts_label': 'Number of reads',
            'xlab': "Lane" if mga_data.from_sequencing else "Data set",
            'ylab': "Number of reads",
            'ymin': 0,
            'ymax': mga_data.max_sequence_count,
            'use_legend': False,
            'tt_percentages': True,
            'hide_zero_cats': False
        }


    def _plot_key(self) -> str:
        '''
        Create a key for the MGA plot. This cannot easily be done through the standard
        mechanisms, so an HTML snippet is created here with little coloured spans that
        form a key to the plot.

        :return An HTML snippet to insert as custom content to the plot.
        :rtype str
        '''
        key_elem_span = '<span class="multiqc_mga_key_element" style="background-color:{}">&nbsp;&nbsp;&nbsp;&nbsp;</span>&nbsp;{}\n'

        key_desc = '<div id="mga_plot_key">\n'
        key_desc += key_elem_span.format(bar_colours.reference.toHtml(), 'Sequenced&nbsp;species/genome')
        key_desc += key_elem_span.format(bar_colours.control.toHtml(), 'Control')
        key_desc += key_elem_span.format(bar_colours.contaminant.toHtml(), 'Contaminant')
        key_desc += key_elem_span.format(bar_colours.unmapped.toHtml(), 'Unmapped')
        key_desc += key_elem_span.format(bar_colours.adapter.toHtml(), 'Adapter')
        key_desc += '</div>\n'

        return key_desc


    def _plot_help(self, mga_data: MGAData) -> str:
        '''
        Build the help text for the plot.

        :param MGAData mga_data: The data loaded from the MGA files.

        :return Markdown help text.
        :rtype str
        '''

        species = set()
        for dataset in mga_data.datasets.values():
            for assignment in dataset.assignments.values():
                if assignment.species is None:
                    species.add(f"{assignment.genome} (unlisted species)")
                else:
                    species.add(assignment.species)

        # See https://stackoverflow.com/questions/36139/how-to-sort-a-list-of-strings
        genomes = natsorted(species, key = cmp_to_key(locale.strcoll))
        number_of_genomes = len(genomes)

        help = f"""
            MGA is a quality control tool for high-throughput sequencing data. It screens
            for contaminants by aligning sequence reads in FASTQ format to a series of
            reference genomes with Bowtie and to a set of adapter sequences using Exonerate.

            MGA samples a subset of the reads, by default 100000, prior to alignment against
            reference genome sequences and adapters, providing a fast screen for unexpected
            sequence content by limiting the computational effort. Sequence reads are trimmed
            to 36 bases by default, prior to alignment to the reference genomes, with the
            aim of minimizing the run time and ensuring consistency of the resulting mapping
            and error rates across sequencing runs with differing read lengths. Full length
            reads are used for matching adapter sequences.

            #### Summary plot

            The summary bar chart displays separate bars for each sample or dataset, one
            representing the genome alignments and the other adapter matches. The genome bar
            contains separate segments for each genome to which sampled reads have been
            assigned and a segment for unmapped reads. The sizes of the segments are based
            on the sampled reads and have been scaled for the total number of sequences. The
            transparency of each segment represents the error or mismatch rate for alignments
            to that genome with segments appearing in bolder colours if the error rate is low
            indicating that these are more likely to be real contaminants; very light or
            transparent segments are of less concern as these largely consist of reads that
            don't align well to any of the available genomes.

            Results are only shown for genomes for which at least
            {assigned_fraction_threshold * 100:.4g}% of reads have assigned or to which at
            least {aligned_fraction_threshold * 100:.4g}% reads align with an average mismatch
            rate of below {error_rate_threshold * 100:.4g}%. An additional bar for adapter
            sequences is only displayed if at least {adapter_threshold_multiplier * 100:.4g}%
            of reads contain adapter sequence matches.

            #### Assigning reads to genomes

            Sequence reads will often align to more than one reference genome. MGA assigns
            reads to the genome with the fewest mismatches. If a read aligns with the same
            smallest number of mismatches to multiple genomes, an expected genome will be
            selected preferentially. Tie-breaks are decided by ranking genomes based on how
            many reads align with fewest mismatches to each.

            #### Reference Genomes

            Sequences were aligned to the following reference genomes ({number_of_genomes} in total).

        """

        # See https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros

        for genome in genomes:
            help += f"* {genome}\n"

        return self._strip_from_lines(help)


    def _main_table_data(self, mga_dataset: MGADataset) -> OrderedDict:
        '''
        Assemble the data for a summary table for a given dataset.

        :param MGADataset mga_dataset: The MGA dataset information.

        :return The table data as required for table.plot.
        :rtype dict
        '''
        dataset_id = mga_dataset.id
        summary = mga_dataset.summary

        number_of_others, other_aligned_count, other_assigned_count = self._count_others(mga_dataset)

        table_data = dict()

        for reference_genome_id, assignment in mga_dataset.assignments.items():
            if number_of_others < 2 or self._accept_genome(assignment):
                table_data[reference_genome_id] = {
                    'species': assignment.species,
                    'aligned_count': assignment.aligned,
                    'aligned_perc': assignment.aligned_frac,
                    'aligned_error': assignment.error_rate,
                    'assigned_count': assignment.assigned,
                    'assigned_perc': assignment.assigned_frac,
                    'assigned_error': assignment.assigned_error_rate
                }

        # Sort into decreasing order of assigned count.
        table_data = OrderedDict(sorted(table_data.items(), key = lambda x: -x[1]['assigned_count']))

        if number_of_others >= 2:
            aligned_frac = 0 if other_assigned_count == 0 else float(other_assigned_count) / float(summary.sampled)
            table_data['Other'] = {
                'species': f"{number_of_others} others",
                'aligned_count': other_assigned_count,
                'aligned_perc': aligned_frac,
                'assigned_count': other_assigned_count,
                'assigned_perc': aligned_frac
            }

        table_data['Unmapped'] = {
            'species': '',
            'aligned_count': summary.unmapped,
            'aligned_perc': summary.unmapped_frac
        }

        table_data['Adapter'] = {
            'species': '',
            'aligned_count': summary.adapter,
            'aligned_perc': summary.adapter_frac
        }

        return table_data


    def _main_table_headers(self) -> OrderedDict:
        '''
        Create the headers for an MGA summary table.

        :return A dictionary of table header information.
        :rtype dict
        '''
        headers = OrderedDict()
        headers['species'] = {
            'title': 'Species/Reference Genome',
            'description': 'Reference genome species',
            'scale': False
        }
        headers['aligned_count'] = {
            'title': 'Aligned',
            'description': 'Number of reads aligned',
            'min': 0,
            'format': '{:d}',
            'scale': False,
            #'shared_key': 'read_count'
        }
        headers['aligned_perc'] = {
            'title': 'Aligned %',
            'description': 'Percentage of reads aligned',
            'min': 0,
            'max': 100,
            'format': '{:,.1%}',
            'scale': False,
            #'shared_key': 'percent_aligned'
        }
        headers['aligned_error'] = {
            'title': 'Error rate',
            'description': 'Aligned error rate',
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
            'scale': False
        }
        headers['assigned_count'] = {
            'title': 'Assigned',
            'description': 'Number of reads assigned',
            'min': 0,
            'format': '{:d}',
            'scale': False,
            #'shared_key': 'read_count'
        }
        headers['assigned_perc'] = {
            'title': 'Assigned %',
            'description': 'Percentage of reads assigned',
            'min': 0,
            'max': 100,
            'format': '{:,.1%}',
            'scale': False,
            #'shared_key': 'percent_aligned'
        }
        headers['assigned_error'] = {
            'title': 'Error rate',
            'description': 'Assigned error rate',
            'min': 0,
            'max': 100,
            'format': '{:,.2%}',
            'scale': False
        }
        return headers


    def _main_table_config(self, dataset_id: str) -> dict:
        '''
        Create the table configuration for a summary table.

        :param str dataset_id: The id of the dataset the table contains.

        :return A dictionary of table configuration parameters.
        :rtype dict
        '''
        return {
            'namespace': 'mga',
            'id': f'mga_stats_table_{dataset_id}',
            'table_title': dataset_id,
            'col1_header': 'Reference ID',
            'no_beeswarm': True,
            'sortRows': False
        }


    def _count_others(self, mga_dataset: MGADataset) -> tuple:
        '''
        Count how many genomes are to go into an "other" category, and count
        their alignments and assignments.

        The "other" classification is those genomes that don't fulfil the acceptance
        criteria as defined by "_accept_genome".

        :param MGADataset mga_dataset: The MGA dataset being considered.

        :return A tuple of the number of "other" genomes, and the total aligned
        and assigned counts for those genomes.
        :rtype tuple
        '''
        number_of_others = 0
        other_aligned_count = 0
        other_assigned_count = 0

        for assignment in mga_dataset.assignments.values():
            if not self._accept_genome(assignment):
                number_of_others = number_of_others + 1
                other_aligned_count = other_aligned_count + assignment.aligned
                other_assigned_count = other_assigned_count + assignment.assigned

        return number_of_others, other_aligned_count, other_assigned_count


    def _accept_genome(self, assignment: MGAAssignment) -> bool:
        '''
        Determine whether an alignment summary qualifies as a genome that is acceptable for
        explicit listing in the plot or table. This is essentially either that the genome is one
        of the expected genomes, or one where the counts are high enough to be significant.

        :param MGAAssignment assignment: The structure containing the genome's assignment information.

        :return if the alignment summary needs to appear in the plot or table, false if it is
        unnecessary.
        :rtype bool
        '''
        if assignment.expected: return True

        return assignment.assigned_frac >= assigned_fraction_threshold or \
               assignment.aligned_frac >= aligned_fraction_threshold and \
               assignment.error_rate < error_rate_threshold


    def _strip_from_lines(self, str: str) -> str:
        '''
        Remove leading and trailing white space from the given multi-line string,
        preserving line breaks.

        :param str: The string to strip.

        :return The string with each line's leading and trailing white space removed.
        :rtype str
        '''
        lines = [ l.strip() for l in str.splitlines() ]
        return "\n".join(lines)

