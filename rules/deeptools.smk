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

def get_multi_bam_summary_input(runTable, cell_line):
    

def get_multi_bam_summary_labels(runTable, cell_line):


rule multi_bam_summary:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        8
    group:
        "deeptools"
    params:
        labels 
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
                                 --labels \
                                 --outFileName {output.npz} 2>{log.logfile}
        """