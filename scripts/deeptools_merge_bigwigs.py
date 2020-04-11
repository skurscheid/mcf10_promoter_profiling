def merge_biwgigs_wrapper():
    try:
        len({input.chip}) < 3
    except ValueError:
        print("Rule can only handle maximum of two input files")
    else:
        if len({input.chip}) == 1:
            cmd = 'ln -sr snakemake.input.chip snakemake.output 2> snakemake.log.logfile'
        elif len({input.chip}) == 2:
            cmd = 'bigwigCompare -b1 snakemake.input.chip[0] -b2 snakemake.input.chip[1] --operation mean -o snakemake.output 2>snakemake.log.logfile'
        shell(cmd)