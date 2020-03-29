rule prefetch:
    version: 
        1
    conda:
        "../envs/sra_tools.yaml"
    threads:
        1
    log:
        "logs/prefetch/{run}.log"
    input:
    output:
        sra_file = temp("sra_download/{cell_line}/{chip_antibody}/{run}.sra")
    shell:
        """
            prefetch {wildcards.run} --output-file {output.sra_file}
        """

rule fastq_dump_se:
    version:
        1
    conda:
        "../envs/sra_tools.yaml"
    threads:
        4
    log:
        "logs/fastq-dump/{run}.log"
    input:
        rules.prefetch.output.sra_file
    output:
        pipe("raw/{cell_line}/{chip_antibody}/{run}.fastq")
    shell:
        """
            fastq-dump {input} {output} 2>{log}
        """

rule pigz_fastq:
    version:
        1
    threads:
        4
    log:
        "logs/pigz_fastq/{run}.log"
    input:
        rules.fastq_dump_se.output
    output:
        "raw/{cell_line}/{chip_antibody}/{run}.fastq.gz"
    shell:
        """
            pigz --stdout {input} > {output} 2>{log}
        """
    