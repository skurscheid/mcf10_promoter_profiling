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

project = config['project']
runTable = pd.read_csv(config['params']['general'][project]['runTable']['file'], sep = ',', index_col='row_id')
library_type = config['library_type']
machine = config['machine']
selected_columns = config['params']['general'][project]['runTable']['selected_columns']

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

rule all_macs2_predictd:
    input:
        expand("macs2/predictd/{file}_predictd.R",
            file = make_targets_from_runTable(runTable, library_type))

rule all_deeptools_plotCorrelation:
    input:
        expand("deeptools/plotCorrelation/{cell_line}.png",
            cell_line = ['MCF10A', 'MCF7', 'AU565'])

rule all_deeptools_bamCoverage:
    input:
        expand("deeptools/bamCoverage/{file}.{suffix}",
               file = make_targets_from_runTable(runTable, library_type),
               suffix = ["bw"])

rule all_merge_bigwigs:
    input:
        expand("deeptools/bigwigCompare/MCF10A/{chip_antibody}_vs_Input.bw",
               chip_antibody = list(runTable.aggregate_column.unique()))

rule all_deeptools_bigwigCompare:
    input:
        expand("deeptools/bigwigCompare/MCF10A/{chip_antibody}_coverage.bw",
               chip_antibody = list(runTable.aggregate_column.unique()))
        
rule all_deeptools_computeMatrix_referencepoint:
    input:
        expand("deeptools/computeMatrix_referencepoint/{cell_line}/{chip_antibody}_{figure}_matrix.gz",
               cell_line = 'MCF10A',
               chip_antibody = list(runTable.aggregate_column.unique()),
               figure = ['Fig1C', 'Fig1D'])

rule all_deeptools_plotProfile:
    input:
        expand("deeptools/plotProfile/{cell_line}/{figure}/{chip_antibody}.pdf",
               cell_line = 'MCF10A',
               chip_antibody = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me3', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H2BK120ub1', 'H3K23ac', 'H3K4me1', 'H4K8ac'],
               figure = ['Fig1A', 'Fig1B', 'Fig1C', 'Fig1D'])

rule all_deeptools_plotHeatmap:
    input:
        expand("deeptools/plotHeatmap/{cell_line}/{figure}/{chip_antibody}.pdf",
               cell_line = 'MCF10A',
               chip_antibody = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me3', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H2BK120ub1', 'H3K23ac', 'H3K4me1', 'H4K8ac'],
               figure = ['Fig1A', 'Fig1B', 'Fig1C', 'Fig1D'])

# test rules
rule test_deeptools_plotProfile:
    input:
        expand("deeptools/plotProfile/{cell_line}/{figure}/{chip_antibody}.pdf",
               cell_line = 'MCF10A',
               chip_antibody = ['H3K27ac'],
               figure = ['Fig1A', 'Fig1B', 'Fig1C', 'Fig1D'])

rule all_align_test:
    input:
        expand("samtools/rmdup/{file}.{suffix}",
               file = make_targets_from_runTable(runTable, library_type)[88],
               suffix = ["bam", "bam.bai"])

rule all_macs2_predictd_test:
    input:
        expand("macs2/predictd/{file}_predictd.R",
            file = make_targets_from_runTable(runTable, library_type)[88])

rule all_deeptools_plotCorrelation_test:
    input:
        expand("deeptools/plotCorrelation/{cell_line}.png",
            cell_line = 'MCF10A')

rule all_deeptools_bamCoverage_test:
    input:
        expand("deeptools/bamCoverage/{file}.{suffix}",
               file = make_targets_from_runTable(runTable, library_type)[88],
               suffix = ["bw"])

rule all_merge_bigwigs_test:
    input:
        "deeptools/merge_bigwigs/MCF10A/H2BK120ub1_coverage.bw"

# includes of rules
include: "rules/other.smk"
include: "rules/sra_download.smk"
include: "rules/fastq_processing.smk"
include: "rules/align.smk"
include: "rules/deeptools.smk"
