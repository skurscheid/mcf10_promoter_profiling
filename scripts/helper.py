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

def fastp_pe_input(wildcards,runTable):
    """function for generating input file names for fastp processing of PE data"""
    return(runTable.loc[runTable['Run'] == wildcards['Run'], ['fq1_path', 'fq2_path']].values.tolist())
