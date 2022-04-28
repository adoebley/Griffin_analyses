downloads/DHS_Index_and_Vocabulary_hg38_WM20190703.txt
	sites from: Meuleman et al, Nature, 2020, Index and biological spectrum of human DNase I hypersensitive sites (https://www.nature.com/articles/s41586-020-2559-3)
	Downloaded 2021-05-12
	wget https://zenodo.org/record/3838751/files/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
	gunzip DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz

hg38_raw/
	this folder contains the raw sites split by tissue type
	script: scripts/split_sites_by_tissue.ipynb

top_10000/
	top 10,000 filtered sites
	script: scripts/get_top_10000.ipynb
	#NOTE: summit column was used for position!