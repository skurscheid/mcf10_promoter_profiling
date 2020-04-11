def merge_biwgigs_wrapper(in_files, out_files, threads, myparams):
    try:
        len(in_files) < 3
    except ValueError:
        print("Rule can only handle maximum of two input files")
    else:
        if len(in_files) == 1:
            cmd = 'ln -sr snakemake.input.chip snakemake.output 2> snakemake.log.logfile'
        elif len(in_files) == 2:
            cmd = 'bigwigCompare -b1 snakemake.input.chip[0] -b2 snakemake.input.chip[1] --operation mean -o snakemake.output 2>snakemake.log.logfile'
        shell(cmd)

merge_bigwigs_wrapper(snakemake.input, snakemake.output, snakemake.threads)