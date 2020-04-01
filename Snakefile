# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import pandas as pd

#configfile: "config.yaml"
report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

include: "scripts/helper.py"

run="[^.]*"

runTable = pd.read_csv(config['runTable'], sep = ',')
library_type = config['library_type']

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_sra_download:
    input:
        expand("raw/{file}.fastq.gz",
               file = make_targets_from_runTable(runTable, library_type))

rule all_fastp:
    input:
        expand("fastp/trimmed/{file}.fastq.gz",
               file = make_targets_from_runTable(runTable, library_type)),
        expand("fastp/report/{file}.fastp.{suffix}",
               file = make_targets_from_runTable(runTable, library_type),
               suffix = ['html', 'json'])
        
include: "rules/other.smk"
include: "rules/sra_download.smk"
include: "rules/fastq_processing"
