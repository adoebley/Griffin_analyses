raw_TCGA/brca_peak_Log2Counts_dedup.txt
	ATAC seq data from the Corces et al 2018 Science 
	download commands:
		wget https://atacseq.xenahubs.net/download/brca/brca_peak_Log2Counts_dedup
		mv brca_peak_Log2Counts_dedup brca_peak_Log2Counts_dedup.txt

raw_TCGA/brca_peak_Log2Counts_dedup.json
	The metadata for brca_peak_Log2Counts_dedup.txt
	download command:
		wget https://atacseq.xenahubs.net/download/brca/brca_peak_Log2Counts_dedup.json

raw_TCGA/aav1898_Data_S1.xlsx
	Supplementary data table 1 from Corces et al 2018 Science

raw_TCGA/sample_info.txt
	txt formatted data from aav1898_Data_S1.xlsx
	This is data from the third sheet in the document ('Sequencing Statistics')

raw_TCGA/aav1898_Data_S2.xlsx
	Supplementary data table 2 from Corces et al 2018 Science

raw_TCGA/BRCA_peak_locations.txt
	txt formatted data from aav1898_Data_S2.xlsx
	This is the data from the second sheet in the document ('BRCA-only_PeakCalls')



#pan cancer
raw_TCGA/TCGA_ATAC_peak_Log2Counts_dedup_sample.txt
	commands:
		wget https://tcgaatacseq.s3.us-east-1.amazonaws.com/latest/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz
		gunzip TCGA_ATAC_peak_Log2Counts_dedup_sample.gz
		mv TCGA_ATAC_peak_Log2Counts_dedup_sample TCGA_ATAC_peak_Log2Counts_dedup_sample.txt

raw_TCGA/TCGA_ATAC_peak_Log2Counts_dedup_sample.json
	commands:
		wget https://tcgaatacseq.s3.us-east-1.amazonaws.com/latest/TCGA_ATAC_peak_Log2Counts_dedup_sample.json

raw_TCGA/pan_cancer_peak_locations
	txt formatted data from aav1898_Data_S2.xlsx
	This is the data from the first sheet in the document ('Pan-Cancer_PeakCalls')


griffin_filter_sites_config/
	config used to filter the site lists
		pan_cancer_peak_locations.txt 
		BRCA_peak_locations.txt

ER_differential
	script and outputs for identifying ER status differential atac sites

