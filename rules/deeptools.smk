__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing HTS data with deepTools
(https://deeptools.readthedocs.io/en/develop/)

For usage, include this in your workflow.
"""

def get_multi_bam_summary_input(runTable, cell_line, input_root, library_type, selected_columns):
    l = []
    sel_rows = runTable[selected_columns[0]] == cell_line
    for index, row in runTable[sel_rows][selected_columns].iterrows():
        l.append('/'.join([input_root, cell_line, row.aggregate_column, library_type, row.Run]) + '.bam')
    return(l)
    
def get_multi_bam_summary_labels(runTable, cell_line):
    l = []
    sel_rows = runTable[selected_columns[0]] == cell_line
    for index, row in runTable[sel_rows][selected_columns].iterrows():
        l.append('_'.join([cell_line, row.aggregate_column, row.Run]))
    return(l)

rule macs2_predictd:
    version:
        1
    conda:
        '../envs/macs2.yaml'
    threads:
        1
    group:
       'deeptools'
    log:
        logfile = 'logs/macs2/predictd/{cell_line}/{chip_antibody}/se/{run}.log'
    params:
        gsize = 'hs'
    input:
        'samtools/rmdup/{cell_line}/{chip_antibody}/se/{run}.bam'
    output:
        directory('macs2/predictd/{cell_line}/{chip_antibody}/se/{run}')
    shell:
        """
            macs2 predictd -i {input}\
                           --gsize {params.gsize}\
                           --outdir {output}\
                           --rfile {wildcards.unit}
        """

rule deeptools_multiBamSummary:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        8
    group:
        "deeptools"
    params:
        labels = get_multi_bam_summary_labels
    log:
        logfile = "logs/deeptools_multiBamSummary/{cell_line}.log"
    input:
        get_multi_bam_summary_input
    output:
        npz = "deeptools/multiBamSummary/{cell_line}.npz"
    shell:
        """
            multiBamSummary bins --bamfiles {input}\
                                 --numberOfProcessors {threads}\
                                 --labels {params.labels}\
                                 --outFileName {output.npz} 2>{log.logfile}
        """
