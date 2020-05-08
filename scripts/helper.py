from pathlib import Path
import os
import yaml
import pandas as pd
from snakemake.io import expand

def make_targets_from_runTable(runTable, library_type):
    t = []
    for index, row in runTable.iterrows():
        chip_antibody = row['chip_antibody'].split()[0]
        if chip_antibody == 'none':
            chip_antibody = 'Input'
        e = list([row['Cell_Line'], chip_antibody, library_type, row['Run']])
        p = "/".join(e)
        t.append(p)
    return(t)

def fastp_targets(units):
    """function for creating snakemake targets for executing fastp rule"""
    t = []
    for index, row in units.iterrows():
        t.append(row['batch'] + "/" + row['sample_id'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

def bowtie2_pe_global_input(wildcards, library_type, project, config):
    """function for generating input file names for bowtie2 alignment of PE data"""
    suffix = config['params']['general'][project]['fastq_suffix']
    t = expand("fastp/trimmed/{cell_line}/{chip_antibody}/pe/{run}.{suffix}",
               cell_line = wildcards['cell_line'],
               chip_antibody = wildcards['chip_antibody'],
               run = wildcards['run'],
               suffix = suffix)
    return(t)
