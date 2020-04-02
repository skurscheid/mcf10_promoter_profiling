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
        "../envs/bowtie2.yaml"
    threads:
        8
    params:
        index = get_index("shiny", config),
        cli_params_global = config['params']['bowtie2']['cli_params_global'],
        samtools_params_global = "-F 4 -bS"
    log:
        log = "logs/bowtie2_global/{cell_line}/{chip_antibody}/se/{run}.log"
    input:
        fq = "fastp/trimmed/{cell_line}/{chip_antibody}/se/{run}.fastq.gz"
    output:
        bam = "bowtie2/align_global/{cell_line}/{chip_antibody}/se/{run}.bam"
    shell:
        """
            bowtie2\
                    -x {params.index}\
                    -p {threads}\
                    -U {input.fq}\
                    {params.cli_params_global}\
                    --rg-id BMG\
                    --rg SM:{wildcards.run}:{wildcards.cell_line}:{wildcards.chip_antibody}\
                    2>> {log.log}\
            | samtools view {params.samtools_params_global} - > {output.bam}
        """