/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10_promoter_profiling/Snakefile\
    --configfile /home/150/sxk150/mcf10_promoter_profiling/config.yaml\
	--use-conda\
	-d ~/data/PromoterSeqCap/public_data/PRJNA336352\
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10_promoter_profiling/cluster.json\
        --keep-going\
	-pr\
	--config project=PRJNA336352 library_type=se machine=gadi\
	-R `/home/150/sxk150/mcf10_promoter_profiling/scripts/cli_snakemake.sh ${1} --lc`\
	${2}
