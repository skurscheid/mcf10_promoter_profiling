/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10_promoter_profiling/Snakefile ${1}\
	--use-conda\
	--rerun-incomplete \
        --local-cores 48\
        --keep-going\
	 -d /home/150/sxk150/data/PromoterSeqCap/public_data/PRJNA336352 \
	-pr ${2}\
	--config project=PRJNA336352 library_type=se machine=gadi
