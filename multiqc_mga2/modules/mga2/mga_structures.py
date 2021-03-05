#!/usr/bin/env python

""" MultiQC Multi Genome Alignment 2 module - supporting data structures """

from distutils.util import strtobool

def _trim_to_none(str: str) -> str:
    '''
    Trim a string to None.

    :param str str: The string to trim.
    :return The trimmed string, or None if the string contained nothing or only whitespace.
    :rtype str
    '''
    if str is None: return None
    str = str.strip()
    return None if len(str) == 0 else str

def _to_bool(str: str) -> str:
    '''
    Convert a string to a bool.

    :param str str: The string convert.
    :return True if the string is a value indicating this, as defined by "strtobool". False otherwise.
    :rtype bool
    '''
    str = _trim_to_none(str)
    return False if str is None else bool(strtobool(str))

def _to_int(str):
    '''
    Convert a string to an integer.

    :param str str: The string convert.
    :return The string converted to a number. If the string is empty or None, return 0.
    :rtype int
    '''
    str = _trim_to_none(str)
    return 0 if str is None else int(str)

def _to_float(str):
    '''
    Convert a string to a float.

    :param str str: The string convert.
    :return The string converted to a number. If the string is empty or None, return 0.0.
    :rtype float
    '''
    str = _trim_to_none(str)
    return 0.0 if str is None else float(str)

def _to_frac(str):
    '''
    Convert a string percentage to a float decimal fraction.

    :param str str: The string convert.
    :return The string converted to a number and divided by 100. If the string is empty or None, return 0.0.
    :rtype float
    '''
    return _to_float(str) / 100.0


class MGAData(object):
    '''
    The top level object for MGA data.

    Contains a dictionary of MGADataset objects
    (keyed by dataset id), a maximum sequence count (the greatest number of reads
    in any data set), a flag indicating whether this is from a sequencing process
    or arbitrary set of files, and a name for the analysis to put in plot titles
    at the like.
    '''

    def __init__(self):
        self.datasets = dict()
        self.max_sequence_count = 0
        self.from_sequencing = False
        self.title = "Multi Genome Alignment Analysis"

    def add_assignment_from_csv(self, assignment):
        '''
        Turn a dictionary of row values from the CSV file into a MGADataset
        object and add it to our "datasets" dictionary.

        :param dict assignment: An assignment information row from the MGA alignment summary file.
        '''
        dataset_id = assignment['id']
        dataset = self.datasets.get(dataset_id)

        if dataset is None:
            dataset = MGADataset(dataset_id)
            self.datasets[dataset_id] = dataset

        dataset.add_assignment_from_csv(assignment)

    def set_summary_from_csv(self, summary):
        '''
        Turn a dictionary of row values from the CSV file into an MGADatasetSummary
        object and set it as the dataset's summary.

        :param dict summary: A summary information row from the MGA summary file.
        '''
        dataset_id = summary['id']
        dataset = self.datasets[dataset_id]
        dataset.set_summary_from_csv(summary)

    def calculate_max_sequences(self):
        '''
        Calculate the maximum number of reads from each of the data sets in this
        object.

        :return The maximum number of reads in any of the data sets, subsequently
        available through the "max_sequence_count" field.
        :rtype Integer
        '''
        self.max_sequence_count = 0
        for dataset in self.datasets.values():
            self.max_sequence_count = max(self.max_sequence_count, dataset.summary.sequences)
        return self.max_sequence_count


class MGADataset(object):
    '''
    The dataset level object for MGA.

    This is the information for a specific data set. It contains a dictionary of
    MGAAssignment objects (keyed by genome), an MGAAssignment object for the unmapped
    counts, and a MGADatasetSummary object for the summary information for this dataset.
    '''
    def __init__(self, id):
        self.id = id
        self.assignments = dict()
        self.unmapped = None
        self.summary = None

    def add_assignment_from_csv(self, assignment):
        '''
        Turn a dictionary of row values from the CSV file into a MGAAssignment
        object and add it to our "assignments" dictionary.

        :param dict assignment: An assignment information row from the MGA alignment summary file.
        '''
        genome = assignment['genome']
        if genome == 'unmapped':
            self.unmapped = MGAAssignment(assignment)
        else:
            self.assignments[genome] = MGAAssignment(assignment)

    def set_summary_from_csv(self, summary):
        '''
        Turn a dictionary of row values from the CSV file into an MGADatasetSummary
        object and set it as our summary.

        :param dict summary: A summary information row from the MGA summary file.
        '''
        self.summary = MGADatasetSummary(summary)


class MGAAssignment(object):
    '''
    Structure for the information from an alignment summary.

    Percentages are turned into decimals.
    '''
    def __init__(self, assignment):
        self.id = _trim_to_none(assignment['id'])
        self.genome = _trim_to_none(assignment['genome'])
        self.species = _trim_to_none(assignment['species'])
        self.expected = _to_bool(assignment['expected'])
        self.control = _to_bool(assignment['control'])
        self.aligned = _to_int(assignment['aligned'])
        self.aligned_frac = _to_frac(assignment['aligned %'])
        self.error_rate = _to_frac(assignment['error rate'])
        self.assigned = _to_int(assignment['assigned'])
        self.assigned_frac = _to_frac(assignment['assigned %'])
        self.assigned_error_rate = _to_frac(assignment['assigned error rate'])

    def _sort_order(self, other):
        '''
        Helper method for sorting MGAAssignment objects for the MGA plot.

        The order is expected before unexpected genomes, then real genomes before controls,
        then sort by descending number of assignment counts.

        :return A number < 0 if this assignment should be before other,
        a number > 0 if other should be before this assignment, and zero
        if no difference can be made.
        '''
        if other is None: return 0

        result = other.expected - self.expected
        if result == 0:
            result = self.control - other.control
        if result == 0:
            result = other.assigned - self.assigned
        return result

    def __eq__(self, other):
        return self._sort_order(other) == 0
    def __ne__(self, other):
        return self._sort_order(other) != 0
    def __lt__(self, other):
        return self._sort_order(other) < 0
    def __gt__(self, other):
        return self._sort_order(other) > 0
    def __le__(self, other):
        return self._sort_order(other) <= 0
    def __ge__(self, other):
        return self._sort_order(other) >= 0


class MGADatasetSummary(object):
    '''
    Structure for the information from a dataset summary.

    Percentages are turned into decimals.
    '''
    def __init__(self, summary):
        self.id = _trim_to_none(summary['id'])
        self.species = _trim_to_none(summary['species'])
        self.controls = _trim_to_none(summary['controls'])
        self.genomes = _trim_to_none(summary['genomes'])
        self.control_genomes = _trim_to_none(summary['control genomes'])
        self.sequences = _to_int(summary['sequences'])
        self.sampled = _to_int(summary['sampled'])
        self.genome_aligned = _to_int(summary['genome aligned'])
        self.genome_aligned_frac = _to_frac(summary['genome aligned %'])
        self.genome_error_rate = _to_frac(summary['genome error rate'])
        self.control_aligned = _to_int(summary['control aligned'])
        self.control_aligned_frac = _to_frac(summary['control aligned %'])
        self.control_error_rate = _to_frac(summary['control error rate'])
        self.unmapped = _to_int(summary['unmapped'])
        self.unmapped_frac = _to_frac(summary['unmapped %'])
        self.adapter = _to_int(summary['adapter'])
        self.adapter_frac = _to_frac(summary['adapter %'])
