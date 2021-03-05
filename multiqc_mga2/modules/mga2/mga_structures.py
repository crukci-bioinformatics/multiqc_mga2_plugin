#!/usr/bin/env python

""" MultiQC Multi Genome Alignment 2 module - supporting data structures """

from distutils.util import strtobool

def _trim_to_none(str):
    if str is None:
        return None
    str = str.strip()
    return None if len(str) == 0 else str

def _to_bool(str):
    str = _trim_to_none(str)
    return False if str is None else bool(strtobool(str))

def _to_int(str):
    str = _trim_to_none(str)
    return 0 if str is None else int(str)

def _to_float(str):
    str = _trim_to_none(str)
    return 0.0 if str is None else float(str)

def _to_frac(str):
    return _to_float(str) / 100.0


class MGAData(object):
    def __init__(self):
        self.datasets = dict()
        self.max_sequence_count = 0
        self.from_sequencing = True
        self.run_id = "210121_A00489_0753_BHYTFKDRXX"

    def add_assignment_from_csv(self, assignment):
        dataset_id = assignment['id']
        dataset = self.datasets.get(dataset_id)

        if dataset is None:
            dataset = MGADataset(dataset_id)
            self.datasets[dataset_id] = dataset

        dataset.add_assignment_from_csv(assignment)
    
    def set_summary_from_csv(self, summary):
        dataset_id = summary['id']
        dataset = self.datasets[dataset_id]
        dataset.set_summary_from_csv(summary)

    def calculate_max_sequences(self):
        self.max_sequence_count = 0
        for dataset in self.datasets.values():
            self.max_sequence_count = max(self.max_sequence_count, dataset.summary.sequences)
        return self.max_sequence_count


class MGADataset(object):
    def __init__(self, id):
        self.id = id
        self.assignments = dict()
        self.unmapped = None
        self.summary = None

    def add_assignment_from_csv(self, assignment):
        genome = assignment['genome']
        if genome == 'unmapped':
            self.unmapped = MGAAssignment(assignment)
        else:
            self.assignments[genome] = MGAAssignment(assignment)
            
    def set_summary_from_csv(self, summary):
        self.summary = MGADatasetSummary(summary)


class MGAAssignment(object):
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
        self.genome_aligned_frac = _to_frac(summary['genome aligned %'])
        self.genome_error_rate = _to_frac(summary['genome error rate'])
        self.control_aligned = _to_int(summary['control aligned'])
        self.control_aligned_frac = _to_frac(summary['control aligned %'])
        self.control_error_rate = _to_frac(summary['control error rate'])
        self.unmapped = _to_int(summary['unmapped'])
        self.unmapped_frac = _to_frac(summary['unmapped %'])
        self.adapter = _to_int(summary['adapter'])
        self.adapter_frac = _to_frac(summary['adapter %'])
