1mb_delfi_HD
	This folder contains scripts and results for creating a PoN for hg38 from delfi data

	First, I created readCounter wigs using scripts/run_readCounter_1mb.sh

	These wigs are in readCounter/

	scripts/1mb_HD_wigs.txt contains a list of the readCounter wigs (this list is used as input for the next step)

	Then I used scripts/run_create_PoN.sh to run the command to create a panel of normals

	The results of this are in results

	This version only uses 215 HD samples (excludes the higher coverage and replicate samples)
