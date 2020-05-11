__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

rule run_fastp_se:
    conda:
        "../envs/fastp.yaml"
    version:
        "3"
    threads:
        4
    log:
        logfile = "logs/fastp/{cell_line}/{chip_antibody}/se/{run}.log"
    input:
        fastq = "raw/{cell_line}/{chip_antibody}/se/{run}.fastq.gz"
    output:
        trimmed = "fastp/trimmed/{cell_line}/{chip_antibody}/se/{run}.fastq.gz",
        report_html = "fastp/report/{cell_line}/{chip_antibody}/se/{run}.fastp.html",
        report_json = "fastp/report/{cell_line}/{chip_antibody}/se/{run}.fastp.json"
    shell:
        "fastp -i {input[0]} -o {output.trimmed} --html {output.report_html} --json {output.report_json} --thread {threads} 2>{log.logfile}"

rule run_fastp_pe:
    conda:
        "../envs/fastp.yaml"
    version:
        "1"
    threads:
        4
    params:
        fastq_suffix = ['.end1.fastq.gz', '.end2.fastq.gz']
    input:
        fastp_pe_input
    output:
        trimmed1 = "fastp/trimmed/{cell_line}/{chip_antibody}/pe/{run}.end1.fastq.gz",
        trimmed2 = "fastp/trimmed/{cell_line}/{chip_antibody}/pe/{run}.end1.fastq.gz",
        report_html = "fastp/report/{cell_line}/{chip_antibody}/pe/{run}.fastp.html",
        report_json = "fastp/report/{cell_line}/{chip_antibody}/pe/{run}.fastp.json"
    shell:
        """
            fastp -i {input.fq1} -I {input.fq2}\
                  -o {output.trimmed1} -O {output.trimmed2}\
                  --detect_adapter_for_pe\
                  --html {output.report_html} --json {output.report_json} --thread {threads} 
        """


