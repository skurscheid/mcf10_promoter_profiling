# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import subprocess as sp

def merge_bigwigs_wrapper(infiles, outfiles, threads, logfile):
        if not isinstance(infiles, list):
            raise TypeError('infiles must be a list')
        if (len(infiles)) > 2:
            raise ValueError("Rule can only handle maximum of two input files")
        elif len(infiles) == 1:
            cmd = ['ln', '-sr', infiles[0], outfiles[0]]
        elif len(infiles) == 2:
            cmd = ['bigwigCompare', '-b1', infiles[0], '-b2', infiles[1], '--operation', 'mean', '-o', outfiles[0], '-p', str(threads)]
        print(cmd)
        sp.run(cmd)

merge_bigwigs_wrapper(snakemake.input, snakemake.output, snakemake.threads, snakemake.log.logfile)
