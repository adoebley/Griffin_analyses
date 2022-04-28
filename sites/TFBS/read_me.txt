Homo_sapiens_meta_clusters.interval
	raw data downloaded from GTRD

	commands to get content in this folder:
	#download GTRD sites 
	wget gtrd.biouml.org/downloads/19.10/chip-seq/Homo%20sapiens_meta_clusters.interval.gz

	#renamed
	mv Homo\ sapiens_meta_clusters.interval.gz Homo_sapiens_meta_clusters.interval.gz

	The copy in this folder is a symlink from here: /fh/fast/ha_g/user/adoebley/projects/griffin_paper/sites/GTRD/Homo_sapiens_meta_clusters.interval

CIS_BP
	download of the CIS_BP database
	(see CIS_BP/read_me.txt for details)

30000_unfiltered_site_CIS_BP_v2
	Top sites for factors that are also in CIS_BP
	Created with projects/griffin_revisions_1/sites/TFBS/scripts/get_unfiltered_with_10000_sites_and_CIS-BP_v2.ipynb
	Other top site lists not included in the github due to size but they can be regenerated with the provided scripts