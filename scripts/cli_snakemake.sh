/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10_promoter_profiling/Snakefile ${1}\
	--use-conda\
	--rerun-incomplete \
        --local-cores 1\
        --keep-going\
	--cluster-config /home/150/sxk150/mcf10_promoter_profiling/cluster.json\
	 -d /g/data/kv78/PromoterSeqCap/public_data/PRJNA336352 \
	-pr ${2}\
	--config runTable=PRJNA336352_SraRunTable.txt library_type=SE
