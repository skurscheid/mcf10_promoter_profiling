# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import pandas as pd

configfile: "config.yaml"
report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

include: "scripts/helper.py"

run="[^.]*"

runTable = pd.read_csv(config['params']['general']['runTable'][config['project']]['file'], sep = ',', index_col='row_id')
library_type = config['library_type']
machine = config['machine']

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

rule all_align:
    input:
        expand("samtools/rmdup/{file}.{suffix}",
               file = make_targets_from_runTable(runTable, library_type),
               suffix = ["bam", "bam.bai"])

rule all_align_test:
    input:
        expand("samtools/rmdup/{file}.{suffix}",
               file = make_targets_from_runTable(runTable, library_type)[88],
               suffix = ["bam", "bam.bai"])

rule all_macs2_predictd:
    input:
        expand("macs2/precitd/{file}",
            file = make_targets_from_runTable(runTable, library_type)[88])
    
include: "rules/other.smk"
include: "rules/sra_download.smk"
include: "rules/fastq_processing.smk"
include: "rules/align.smk"
include: "rules/deeptools.smk"
