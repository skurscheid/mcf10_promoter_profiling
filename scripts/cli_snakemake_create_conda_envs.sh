/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10_promoter_profiling/Snakefile ${1}\
	--use-conda \
	-d /g/data/kv78/PromoterSeqCap/public_data/PRJNA336352 \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10_promoter_profiling/cluster.json\
        --keep-going \
	-pr --create-envs-only\
	--config runTable=PRJNA336352_SraRunTable.txt
