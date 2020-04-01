snakemake -s /home/sebastian/Development/mcf10_promoter_profiling/Snakefile ${1}\
	--use-conda\
	--rerun-incomplete \
        --local-cores 32\
        --keep-going\
	 -d /home/sebastian/Data/Tremethick/PromoterSeqCap/public_data/PRJNA336352 \
	-pr ${2}\
	--config runTable=PRJNA336352_SraRunTable_MFC10A.txt library_type=se
