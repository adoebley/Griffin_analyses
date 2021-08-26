Homo_sapiens_meta_clusters.interval
	raw data downloaded from GTRD

	commands to get content in this folder:
	#download GTRD sites 
	wget gtrd.biouml.org/downloads/19.10/chip-seq/Homo%20sapiens_meta_clusters.interval.gz

	#renamed
	mv Homo\ sapiens_meta_clusters.interval.gz Homo_sapiens_meta_clusters.interval.gz

	
Unfiltered sites
	This is the data from Homo_sapiens_meta_clusters.interval split into individual TFs
	script: scripts/get_individual_TFs.ipynb

griffin_filter_sites_config
	griffin_filter_sites snakelike
	used to filter the sites by capability
	full outputs not included due to size

CIS_BP
	download of the CIS_BP database
	actual download not included here
	(see CIS_BP_read_me.txt for details of how the download was performed)

10000_filtered_sites_CIS_BP
	Copy of the top 10,000 filtered sites for TFs with cis_bp motif
	script: scripts/get_factors_with_10000_sites_and_CIS-BP.ipynb
	