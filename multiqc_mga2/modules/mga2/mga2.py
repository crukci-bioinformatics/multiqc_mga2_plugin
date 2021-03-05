#!/usr/bin/env python

""" MultiQC Multi Genome Alignment 2 module """

import csv
import locale
import logging
import operator
import os

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
    'adapter': Colour.fromBytes(255, 102, 255)
})

max_alpha = 1.0
min_alpha = 0.4 # Is 0.1 in the Java version.
max_error = 0.01
min_error = 0.0025

assigned_fraction_threshold = 0.01
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
                   written by Matthew Eldridge at the Cancer Research UK Cambridge Institute."""
        )

        # TODO - work out if we're part of a sequencing run.

        # Add to self.css to be included in template
        self.css = { 'assets/css/multiqc_mga2.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_mga2.css') }

        mga_data = self.load_mga_files()

        self.create_mga_reports(mga_data)


    def load_mga_files(self):
        """
        Find and load the content of MGA2 mga_alignment_summary.csv files.

        :return A dictionary of MGA dataset ids to further dictionaries of summary
        dictionaries (the rows from the files) keyed by genome ids. Or
        dataset id -> genome id -> stats for that genome in the dataset.
        :rtype dict
        """

        mga_data = MGAData()

        for mgafile in self.find_log_files('mga2', filecontents=False, filehandles=False):
            self._read_mga2_csv_file(mga_data, mgafile)

        if len(mga_data.datasets) == 0:
            raise UserWarning

        log.info("Found {} MGA2 datasets".format(len(mga_data.datasets)))

        mga_data.calculate_max_sequences()

        return mga_data


    def _read_mga2_csv_file(self, mga_data, mgafile):
        """
        Read a MGA2 alignment summary file.

        :param dict mga_data: A dictionary into which found content can be added, keyed by dataset id.
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
                dataset_id = assignment['id']
                dataset = mga_data.datasets.get(dataset_id)
                if dataset is None:
                    dataset = MGADataset(dataset_id)
                    mga_data.datasets[dataset_id] = dataset

                genome = assignment['genome']
                if genome == 'unmapped':
                    dataset.unmapped = MGAAssignment(assignment)
                else:
                    dataset.assignments[genome] = MGAAssignment(assignment)

        with open(summaryfile_path, newline='') as fh:
            log.debug(f"Found file {summaryfile_path}")
            reader = csv.DictReader(fh, dialect = 'unix')
            for summary in reader:
                dataset_id = summary['id']
                mga_data.datasets[dataset_id].summary = MGADatasetSummary(summary)


    def create_mga_reports(self, mga_data):
        '''
        Build the sections for the overall report

        :param mga_data: The data loaded from the MGA files. It is a dictionary with the
        dataset id as the key and a further dictionary keyed by genome id to that pairing's
        full data.

        :param mga_summaries: A dictionary of the summary information, keyed by dataset id.
        '''

        # The plot.

        plot_data, plot_categories = self._plot_data(mga_data)

        self.add_section(
            name = 'Sequencing Results',
            anchor = 'mga_sequencing_results',
            # description = self._plot_description(mga_data),
            content = self._plot_key(),
            helptext = self._plot_help(mga_data),
            plot = bargraph.plot(plot_data, plot_categories, self._plot_config(mga_data))
        )

        # The summary table (one per dataset).
        
        for dataset_id, mga_dataset in mga_data.datasets.items():

            self.add_section(
                name = f'Lane {dataset_id} Statistics' if mga_data.from_sequencing else f'Dataset "{mga_summary}" Statistics',
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


    def _plot_data(self, mga_data):
        '''
        Assemble the plot information.

        :param mga_summaries: A dictionary of the summary information, keyed by dataset id.

        :return A tuple of bargraph plot data and categories. Both values are dictionaries keyed
        by dataset id.
        '''

        bar_data = OrderedDict()
        categories = OrderedDict()

        for dataset_id, mga_dataset in mga_data.datasets.items():
            dataset_summary = mga_dataset.summary

            sequence_count = dataset_summary.sequences
            sampled_count = dataset_summary.sampled
            adapter_count = dataset_summary.adapter
            unmapped_count = dataset_summary.unmapped

            sampled_to_sequenced = float(sequence_count) / float(sampled_count)

            dataset_bar_data = dict()
            dataset_categories = dict()

            for reference_genome_id, assignment in mga_dataset.assignments.items():
                if self._accept_genome(assignment):
                    aligned_count = assignment.aligned
                    aligned_error = assignment.error_rate
                    assigned_count = assignment.assigned
                    assigned_error = assignment.assigned_error_rate

                    reference_genome_name = assignment.species

                    category_id = f"{dataset_id}.{reference_genome_id}"

                    colour = bar_colours.contaminant
                    if assignment.control:
                        colour = bar_colours.control
                    elif assignment.expected:
                        colour = bar_colours.reference

                    # log.debug("Dataset {} - {}".format(dataset_id, reference_genome_id))
                    # log.debug("max_alpha - (max_alpha - min_alpha) * (assigned_error - min_error) / (max_error - min_error)")
                    # log.debug("{} - ({} - {}) * ({} - {}) / ({} - {})".format(max_alpha, max_alpha, min_alpha, assigned_error, min_error, max_error, min_error))
                    # log.debug("{} - {} * {} / {}".format(max_alpha, max_alpha - min_alpha, assigned_error - min_error, max_error - min_error))

                    alpha = max_alpha - (max_alpha - min_alpha) * (assigned_error - min_error) / (max_error - min_error)

                    # log.debug("alpha = {}".format(alpha))

                    alpha = max(min_alpha, min(max_alpha, alpha))

                    # log.debug("capped alpha = {}".format(alpha))

                    if assigned_count >= 100:
                        log.debug("{}\t{}\t{}\t{}".format(reference_genome_id, assigned_count, aligned_error * 100.0, alpha))

                    dataset_bar_data[category_id] = int(assigned_count * sampled_to_sequenced)

                    dataset_categories[category_id] = {
                        'name': reference_genome_name,
                        'color': colour.applyAlpha(alpha).toHtml()
                    }

            # Sort into decreasing order of assigned count.
            dataset_bar_data = OrderedDict(sorted(dataset_bar_data.items(), key = lambda x: -x[1]))

            bar_data[dataset_id] = dataset_bar_data

            log.debug(f"Unmapped count: {unmapped_count} / {sampled_count}")

            # Add the unmapped count.
            category_id = f"{dataset_id}.unmapped"
            dataset_bar_data[category_id] = int(unmapped_count * sampled_to_sequenced)
            dataset_categories[category_id] = {
                'name': 'Unmapped',
                'color': bar_colours.unmapped.toHtml()
            }

            # The order of categories matters for the order in the plot, so add them to the
            # overall categories dictionary in the order they are in for the bar data.
            for category_id in dataset_bar_data.keys():
                categories[category_id] = dataset_categories[category_id]

            log.debug(f"Adapter count: {adapter_count} / {sampled_count}")

            if adapter_count >= sampled_count * adapter_threshold_multiplier:
                dataset_adapter_id = f"{dataset_id}A" if mga_data.from_sequencing else f"{dataset_id} adapter"
                category_id = f"{dataset_id}.adapter"
                dataset_bar_data = { category_id: int(adapter_count * sampled_to_sequenced) }
                bar_data[dataset_adapter_id] = dataset_bar_data
                categories[category_id] = {
                    'name': 'Adapter',
                    'color': bar_colours.adapter.toHtml()
                }


        log.debug(f"Bar data = {bar_data}")
        log.debug(f"Categories = {categories}")

        return bar_data, categories


    def _plot_config(self, mga_data):
        '''
        Create the plot configuration for the main MGA plot.

        :param run_info: The information for MGA data with one run id, as created by
        _merge_run_summaries

        :return A dictionary of plot configuration parameters.
        '''
        return {
            'id': "mga_plot_{}".format(mga_data.run_id.replace(' ', '_')),
            'title': f"Multi Genome Alignment: {mga_data.run_id}",
            'cpswitch_counts_label': 'Number of reads',
            'xlab': "Lane" if mga_data.from_sequencing else "Data set",
            'ylab': "Number of reads",
            'ymin': 0,
            'ymax': mga_data.max_sequence_count,
            'use_legend': False,
            'tt_percentages': True,
            'hide_zero_cats': False
        }


    def _plot_key(self):
        '''
        Create a key for the MGA plot. This cannot easily be done through the standard
        mechanisms, so an HTML snippet is created here with little coloured spans that
        form a key to the plot.

        :return An HTML snippet to insert as custom content to the plot.
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


    def _plot_help(self, mga_data):
        '''
        Build the help text for the plot.

        :param run_info: The information for MGA data with one run id, as created by
        _merge_run_summaries

        :return Markdown help text.
        '''

        species = set()
        for dataset in mga_data.datasets.values():
            for assignment in dataset.assignments.values():
                species.add(assignment.species)

        # See https://stackoverflow.com/questions/36139/how-to-sort-a-list-of-strings
        genomes = natsorted(species, key = cmp_to_key(locale.strcoll))
        number_of_genomes = len(genomes)

        help = f"""
            Sequences were sampled, trimmed to ##TRIM LENGTH## bases starting from position ##TRIM START##,
            and mapped to {number_of_genomes} reference genomes (see list below) using Bowtie. Sequences containing
            adapters were found by ungapped alignment of the full length sequence to a set of known adapter and
            primer sequences using Exonerate.

        """

        if mga_data.from_sequencing:
            help += f"""
                Some lanes may have an addition bar with the suffix "A" (e.g. "1A" for lane 1).
                This is the adapter when the number of adapter reads is significant, equal to or
                above {adapter_threshold_multiplier * 100:.4g}% of the number of sampled reads.
            """
        else:
            help += f"""
                Some data sets may have an addition bar with the suffix "adapter".
                This is the adapter when the number of adapter reads is significant for the
                data set, equal to or above {adapter_threshold_multiplier * 100:.4g}% of the number of sampled reads.
            """

        # See https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros

        help += f"""

            #### Alignment Details

            Reference genomes are sorted according to how many sequence reads have been assigned to each.
            Separate entries are given for reference genomes for which at least {assigned_fraction_threshold * 100:.4g}%
            of reads have been assigned or for which at least {aligned_fraction_threshold * 100:.4g}% of reads align
            with an average mismatch or error rate of below {error_rate_threshold * 100:.4g}%.

            In addition to the total number of reads aligning to each reference genome and the average error
            rate for those alignments, details are also provided for the the number of reads aligning uniquely
            to the reference genome and and the associated error rate for those unique reads.

            The 'Best' column and accompanying error rate refer to those reads that align
            to a given reference genome with fewer mismatches than to other genomes.
            Included are reads that align uniquely to the genome and those for which there
            is a tie-break, aligning equally well to more than one genome. In the latter
            case, a read will contribute to the 'Best' column for each of the genomes to
            which it aligns the best. The sum of read counts in the 'Best' column can
            therefore exceed the total number of sampled reads.

            Reads that align uniquely to a genome are assigned unambiguously to that genome.
            Tie-breaks in which reads align equally well to multiple genomes are assigned
            preferentially to one of the expected species for this sequence dataset and/or
            to the genome with the highest read count in the 'Best' column.

            Note that because reads are trimmed prior to alignment with Bowtie, it is possible for a read to be
            counted both as aligned to one or more of the reference genomes and among the reads with adapter
            content. The adapter will most likely be present in the portion of the read that has been trimmed.

            #### Reference Genomes

            Sequences were aligned to the following reference genomes ({number_of_genomes} in total).

        """

        for genome in genomes:
            help += f"* {genome}\n"

        return self._strip_from_lines(help)


    def _main_table_data(self, mga_dataset):
        '''
        Assemble the data for a summary table for a given dataset.

        :param mga_summary: An XML tree for a dataset.

        :return The table data as required for table.plot.
        '''
        dataset_id = mga_dataset.id
        summary = mga_dataset.summary

        number_of_others = 0
        other_aligned_count = 0
        other_assigned_count = 0
        table_data = dict()

        for assignment in mga_dataset.assignments.values():
            if not self._accept_genome(assignment):
                number_of_others = number_of_others + 1
                other_aligned_count = other_aligned_count + assignment.aligned
                other_assigned_count = other_assigned_count + assignment.assigned

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
            table_data['Other'] = {
                'species': f"{number_of_others} others",
                'aligned_count': other_assigned_count,
                'aligned_perc': float(other_assigned_count) / float(summary.sampled),
                'assigned_count': other_assigned_count,
                'assigned_perc': float(other_assigned_count) / float(summary.sampled)
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


    def _main_table_headers(self):
        '''
        Create the headers for an MGA summary table.

        :return A dictionary of table header information.
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


    def _main_table_config(self, dataset_id):
        '''
        Create the table configuration for a summary table.

        :param dataset_id: The id of the dataset the table contains.

        :return A dictionary of table configuration parameters.
        '''
        return {
            'namespace': 'mga',
            'id': f'mga_stats_table_{dataset_id}',
            'table_title': dataset_id,
            'col1_header': 'Reference ID',
            'no_beeswarm': True,
            'sortRows': False
        }


    def _accept_genome(self, assignment):
        '''
        Determine whether an alignment summary qualifies as a genome that is acceptable for
        a sample. This is essentially either that the genome is one of the genomes that apply
        to any sample in the data set, or one where the counts are high enough to be significant.

        :param assignment: The dictionary containing the information for a dataset and genome.

        :return if the alignment summary needs to appear in the plot or table, false if it is
        unnecessary.
        '''
        if assignment.expected or assignment.control:
            return True

        return assignment.assigned_frac >= assigned_fraction_threshold or \
               assignment.aligned_frac >= aligned_fraction_threshold and \
               assignment.error_rate < error_rate_threshold


    def _strip_from_lines(self, str):
        '''
        Remove leading and trailing white space from the given multi-line string,
        preserving line breaks.

        :param str: The string to strip.

        :return The string with each line's leading and trailing white space removed.
        '''
        lines = [ l.strip() for l in str.splitlines() ]
        return "\n".join(lines)

