#this analysis was run using the DGE tool at https://xenabrowser.net/
(http://xena.ucsc.edu/ then click Launch Xena)

Tutorial is here: https://ucsc-xena.gitbook.io/project/overview-of-features/differential-gene-expression

select study: 
	TCGA TARGET GTEX

	select variable 1:
		main category
	select variable 2: (not used)
		study

	under variable 1 click the 3 dots and select 'differential expression'
		Subgroup1:
			GTEX Blood (337 samples)
		Subgroup2:
			TCGA_Breast_Invasive_Carcinoma (1099 samples)

		submit analysis
		download results: 
			DEG_results_GTEX_Blood vs. TCGA_Breast_Invasive_Carcinoma.csv



	under variable 1 click the 3 dots and select 'differential expression'
		Subgroup1:
			GTEX Blood (337 samples)
		Subgroup2:
			GTEX Breast (179 samples)

		submit analysis
		download results:
			DEG_results_GTEX_Blood vs. GTEX_Breast.csv

select study: 
	GTEX

	select variable 1:
		Body site detail
	select variable 2: (not used)
		_primary_site

	under variable 1 click the 3 dots and select 'differential expression'
		Subgroup 1:
			Whole Blood (456 samples)
		Subgroup 2: 
			Breast - Mammary Tissue (221 samples)

		submit analysis
		download results: 
			DEG_results_Whole_Blood vs. Breast_Mammary_Tissue.csv

Select dataset:
	TCGA Breast Cancer (BRCA)

	select variable 1:
		ER_Status_nature2012

	select variable 2: (not used)
		sample type

	under variable 1 click the 3 dots and select 'differential expression'
		Subgroup 1:
			Positive
		Subgroup 2:
			Negative	

		submit analysis
		download results: 	
			DEG_results_Positive vs. Negative.csv
