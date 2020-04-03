snakemake -s /home/150/sxk150/mcf10_promoter_profiling/Snakefile ${1}\
	--use-conda\
	--rerun-incomplete \
        --local-cores 48\
        --keep-going\
	 -d /home/150/sxk150/data/PromoterSeqCap/public_data/PRJNA336352 \
	-pr ${2}\
	--config runTable=PRJNA336352_SraRunTable_MFC10A.txt library_type=se
