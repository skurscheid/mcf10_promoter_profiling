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

# input & other functions
def get_multi_bam_summary_input(wildcards):
    selected_columns = config['params']['general']['runTable'][config['project']]['selected_columns']
    library_type = config['library_type']
    cell_line = wildcards['cell_line']
    l = []
    sel_rows = runTable[selected_columns[0]] == cell_line
    for index, row in runTable[sel_rows][selected_columns].iterrows():
        l.append('/'.join(["samtools/rmdup", cell_line, row.aggregate_column, library_type, row.Run]) + '.bam')
    return(l)

def get_multi_bam_summary_labels(wildcards):
    selected_columns = config['params']['general']['runTable'][config['project']]['selected_columns']
    cell_line = wildcards['cell_line']
    l = []
    sel_rows = runTable[selected_columns[0]] == cell_line
    for index, row in runTable[sel_rows][selected_columns].iterrows():
        l.append('_'.join([cell_line, row.aggregate_column, row.Run]))
    return(l)

def get_bigwigCompare_inputs(wildcards):
    selected_columns = config['params']['general']['runTable'][config['project']]['selected_columns']
    library_type = config['library_type']
    cell_line = wildcards['cell_line']
    chip_antibody = wildcards['chip_antibody']
    l = []
    selection = runTable[(runTable['aggregate_column'] == chip_antibody) & (runTable[selected_columns[0]] == cell_line)]
    for index, row in selection.iterrows():
        l.append('/'.join(["deeptools/bamCoverage", cell_line, row.aggregate_column, library_type, row.Run]) + '.bw')
    return(l)

# actual rules
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
        logfile = 'logs/macs2/predictd/{cell_line}/{chip_antibody}/{library_type}/{run}.log'
    params:
        gsize = 'hs'
    input:
        'samtools/rmdup/{cell_line}/{chip_antibody}/{library_type}/{run}.bam'
    output:
        file = ('macs2/predictd/{cell_line}/{chip_antibody}/{library_type}/{run}_predictd.R')
    shell:
        """
            macs2 predictd -i {input}\
                           --gsize {params.gsize}\
                           --rfile {output.file} 2>{log.logfile}
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
        labels = get_multi_bam_summary_labels,
        extendReads = 200
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
                                 --extendReads {params.extendReads}\
                                 --outFileName {output.npz} 2>{log.logfile}
        """

rule deeptools_plotCorrelation:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        1
    group:
        "deeptools"
    params:
    log:
        logfile = "logs/deeptools_plotCorrelation/{cell_line}.log"
    input:
        rules.deeptools_multiBamSummary.output.npz
    output:
        png = "deeptools/plotCorrelation/{cell_line}_heatmap.png"
    shell:
        """
            plotCorrelation -in {input} --whatToPlot heatmap --corMethod pearson -o {output}
        """

rule deeptools_bamCoverage:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        8
    group:
        "deeptools"
    params:
        extendReads = 200,
        binSize = 30,
        smoothLength = 10,
        effectiveGenomeSize = config['params']['deeptools']['genome_size']['GRCh37_hg19_UCSC'],
        normalizeUsing = 'RPKM'
    log:
        logfile = "logs/deeptools_bamCoverage/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    input:
        rules.bam_rmdup.output
    output:
        "deeptools/bamCoverage/{cell_line}/{chip_antibody}/{library_type}/{run}.bw"
    shell:
        """
            bamCoverage -b {input}\
                        --numberOfProcessors {threads}\
                        --effectiveGenomeSize {params.effectiveGenomeSize}\
                        --normalizeUsing {params.normalizeUsing}\
                        --extendReads {params.extendReads}\
                        --binSize {params.binSize}\
                        --smoothLength {params.smoothLength}\
                        --outFileName {output} 2>{log.logfile}
        """

rule merge_bigwigs:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        4
    group:
        "deeptools"
    params:
    log:
        logfile = "logs/merge_bigwigs/{cell_line}/{chip_antibody}_coverage.log"
    input:
        chip = get_bigwigCompare_inputs
    output:
        "deeptools/merge_bigwigs/{cell_line}/{chip_antibody}_coverage.bw"
    script:
        "../scripts/deeptools_wrappers.py"

#rule deeptools_computeMatrix:
#    version:
#        1
#    conda:
#        "../envs/deeptools.yaml"
#    threads:
#        1
#    group:
#        "deeptools"
#    params:
#    log:
#        logfile = "logs/deeptools_computeMatrix/{cell_line}.log"
#    input:
