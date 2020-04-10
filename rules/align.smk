__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for aligning reads with bowtie2
(https://github.com/BenLangmead/bowtie2)

For usage, include this in your workflow.
"""

def get_index(machine, config):
    """ returns path to index"""
    return config["params"]["bowtie2"]["index"][machine]

singularity: "docker://skurscheid/snakemake_baseimage:0.2"

rule bowtie2_se_global:
    """ runs alignment of single-end fastq file, modified parameters specific for HiC data"""
    conda:
        "../envs/alignment.yaml"
    threads:
        12
    group:
        "alignment"
    params:
        index = get_index(machine, config),
        cli_params_global = config['params']['bowtie2']['cli_params_global'],
        samtools_params_global = "-F 4 -bS"
    log:
        logfile = "logs/bowtie2_global/{cell_line}/{chip_antibody}/se/{run}.log"
    input:
        fq = "fastp/trimmed/{cell_line}/{chip_antibody}/se/{run}.fastq.gz"
    output:
        bam = temp("bowtie2/align_global/{cell_line}/{chip_antibody}/se/{run}.bam")
    shell:
        """
            export cli_threads=$(expr {threads} - 2);\
            bowtie2\
                    -x {params.index}\
                    -p $cli_threads\
                    -U {input.fq}\
                    {params.cli_params_global}\
                    --rg-id BMG\
                    --rg SM:{wildcards.run}:{wildcards.cell_line}:{wildcards.chip_antibody}\
                    2>> {log.logfile}\
            | samtools view {params.samtools_params_global} - > {output.bam}
        """

rule bam_quality_filter:
    conda:
        "../envs/alignment.yaml"
    version:
        "1.0"
    group:
        "alignment"
    log:
        logfile = "logs/samtools/quality_filtered/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    params:
        qual = config["params"]["general"]["alignment_quality"]
    input:
        rules.bowtie2_se_global.output
    output:
        temp("samtools/quality_filtered/{cell_line}/{chip_antibody}/{library_type}/{run}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output} 2>{log.logfile}"

rule bam_sort:
    conda:
        "../envs/alignment.yaml"
    version:
        "1.0"
    threads:
        4
    group:
        "alignment"
    log:
        logfile = "logs/samtools/sort/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    input:
        rules.bam_quality_filter.output
    output:
        temp("samtools/sort/{cell_line}/{chip_antibody}/{library_type}/{run}.bam")
    shell:
        "samtools sort -@ {threads} {input} -T {wildcards.run}.sorted -o {output}"

rule bam_mark_duplicates:
    conda:
        "../envs/alignment.yaml"
    version:
        "1.0"
    group:
        "alignment"
    log:
        logfile = "logs/picardTools/MarkDuplicates/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    threads:
        4
    params:
        temp = config["params"]["general"]["temp_dir"]["shiny"]
    input:
        rules.bam_sort.output
    output:
        out= temp("picardTools/MarkDuplicates/{cell_line}/{chip_antibody}/{library_type}/{run}.bam"),
        metrics = "picardTools/MarkDuplicates/{cell_line}/{chip_antibody}/{library_type}/{run}.metrics.txt"
    shell:
        """
            picard MarkDuplicates -XX:ParallelGCThreads={threads} -Xms2g -Xmx8g\
            INPUT={input}\
            OUTPUT={output.out}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output.metrics} 2>{log.logfile}
        """

rule bam_rmdup:
    conda:
        "../envs/alignment.yaml"
    group:
        "alignment"
    log:
        logfile = "logs/samtools/rmdup/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    input:
        rules.bam_mark_duplicates.output.out
    output:
        "samtools/rmdup/{cell_line}/{chip_antibody}/{library_type}/{run}.bam"
    shell:
        "samtools rmdup -s {input} {output} 2>{log.logfile}"

rule bam_index:
    conda:
        "../envs/alignment.yaml"
    group:
        "alignment"
    log:
        logfile = "logs/samtools/index/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    input:
        rules.bam_rmdup.output
    output:
        "samtools/rmdup/{cell_line}/{chip_antibody}/{library_type}/{run}.bam.bai"
    shell:
        "samtools index {input} {output} 2>{log.logfile}"
