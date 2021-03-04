#!/usr/bin/env python

""" MultiQC Multi Genome Alignment 2 module - supporting data structures """

from distutils.util import strtobool

def _to_bool(str):
    str = _trim_to_none(str)
    return False if str is None else bool(strtobool(str))

def _to_int(str):
    str = _trim_to_none(str)
    return 0 if str is None else int(str)

def _to_float(str):
    str = _trim_to_none(str)
    return 0.0 if str is None else float(str)

def _trim_to_none(str):
    if str is None:
        return None
    str = str.strip()
    return None if len(str) == 0 else str


class MGAData(object):
    def __init__(self):
        self.datasets = dict()
        self.max_sequence_count = 0
        self.from_sequencing = True
        self.run_id = "210121_A00489_0753_BHYTFKDRXX"

    def calculate_max_sequences(self):
        self.max_sequence_count = 0
        for dataset in self.datasets.values():
            self.max_sequence_count = max(self.max_sequence_count, dataset.summary.sequences)
        return self.max_sequence_count


class MGADataset(object):
    def __init__(self):
        self.assignments = dict()
        self.summary = None


class MGAAssignment(object):
    def __init__(self, assignment):
        self.id = _trim_to_none(assignment['id'])
        self.genome = _trim_to_none(assignment['genome'])
        self.species = _trim_to_none(assignment['species'])
        self.expected = _to_bool(assignment['expected'])
        self.control = _to_bool(assignment['control'])
        self.aligned = _to_int(assignment['aligned'])
        self.aligned_pc = _to_float(assignment['aligned %'])
        self.error_rate = _to_float(assignment['error rate'])
        self.assigned = _to_int(assignment['assigned'])
        self.assigned_pc = _to_float(assignment['assigned %'])
        self.assigned_error_rate = _to_float(assignment['assigned error rate'])


class MGADatasetSummary(object):
    def __init__(self, summary):
        self.id = _trim_to_none(summary['id'])
        self.species = _trim_to_none(summary['species'])
        self.controls = _trim_to_none(summary['controls'])
        self.genomes = _trim_to_none(summary['genomes'])
        self.control_genomes = _trim_to_none(summary['control genomes'])
        self.sequences = _to_int(summary['sequences'])
        self.sampled = _to_int(summary['sampled'])
        self.genome_aligned = _to_int(summary['genome aligned'])
        self.genome_aligned_pc = _to_float(summary['genome aligned %'])
        self.genome_error_rate = _to_float(summary['genome error rate'])
        self.control_aligned = _to_int(summary['control aligned'])
        self.control_aligned_pc = _to_float(summary['control aligned %'])
        self.control_error_rate = _to_float(summary['control error rate'])
        self.unmapped = _to_int(summary['unmapped'])
        self.unmapped_pc = _to_float(summary['unmapped %'])
        self.adapter = _to_int(summary['adapter'])
        self.adapter_pc = _to_float(summary['adapter %'])
