This data is from:
Satpathy AT, Granja JM, Yost KE, Qi Y et al. Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nat Biotechnol 2019 Aug;37(8):925-936. PMID: 31375813


sites_hg19/GSE129785_scATAC-Hematopoiesis-All.peaks.txt
	Merged hematopoietic ATAC peaks from single cell data
	
	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785

	download commands:
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129785/suppl/GSE129785%5FscATAC%2DHematopoiesis%2DAll%2Epeaks%2Etxt%2Egz

		gunzip GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz 
		
		mv GSE129785_scATAC-Hematopoiesis-All.peaks.txt GSE129785_scATAC-Hematopoiesis-All.peaks_hg19.txt

sites_hg19/GSE129785_scATAC-Hematopoiesis-All.peaks_hg19.bed
	Convereted the raw download into bed format using scripts/convert_heme_to_bed.ipynb

sites_hg38/GSE129785_scATAC-Hematopoiesis-All.peaks_hg38.bed
	primary output of scripts/liftover_heme.sh

sites_hg38/GSE129785_scATAC-Hematopoiesis-All.peaks_hg38.unmapped.bed
	Another output of scripts/liftover_heme.sh
	This file contains the regions that couldn't be lifted over