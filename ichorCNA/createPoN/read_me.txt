1mb_delfi_HD
	This folder contains scripts and results for creating a PoN for hg38 from delfi data

	First, I created readCounter wigs using scripts/run_readCounter_1mb.sh

	These wigs are in readCounter/

	scripts/1mb_HD_wigs.txt contains a list of the readCounter wigs (this list is used as input for the next step)

	I made a copy of createPanelOfNormals.R (scripts/mod_createPanelOfNormals.R) so that I could correct the parameters on one of the functions
		This script is copied from the ichor version: ichorCNA_2020-11-20

	Then I used scripts/run_create_PoN.sh to run the command to create a panel of normals

	The results of this are in results
