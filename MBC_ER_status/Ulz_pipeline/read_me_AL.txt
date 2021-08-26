downloads
	this is an edited copy of the Ulz pipeline (https://github.com/PeterUlz/TranscriptionFactorProfiling)
	(original unedited download is in: /fh/fast/ha_g/user/adoebley/downloads/TranscriptionFactorProfiling)
	original download was around 10-22-2019
	modified things as needed to work on our server
	I tried to mark all edits with ‘AL edit or AL comment’ 

	Did not copy over Ref/GTRD or Ref/GTRD_slop to this location because they are large and not used for this analysis

Software
	This folder contains symlinks of software
	ln -s /app/software/BEDTools/2.29.2-GCC-8.3.0/bin/bedtools

Scripts
	This is a symlink of the Scripts directory in the downloads folder
	ln -s downloads/Scripts/ .

Ref
	This is a symlink of the Ref directory in the downloads folder
	ln -s downloads/Ref/ .


write_scripts
	contains scripts (made by Katharine or me) to write the slurm scripts that run the analysis

config
	contains the sample yamls (this analysis only works on hg19)

sh_scripts
	created by the write_scripts_new.py script
	contains scripts to run the analysis

merge
	contains the merged output
	
#PROCEDURE FOR RUNNING THE ANALYSIS

Write scripts
	python write_scripts/write_scripts_new.py -i config/MBC_0.1X_0.1TFx_hg19.yaml -n MBC_ULP

	OUTPUT (run these commands in subsequent steps):
		Commands to run scripts:
		launchSlurm -p campus-new -b sh_scripts/trim_align_MBC_ULP.sh
		python write_scripts/reformat_log2_file.py -i <path/to/results/ichorCNA>
		launchSlurm -p campus-new -b sh_scripts/run_tf_analysis_MBC_ULP.sh
		launchSlurm -p campus-new -b sh_scripts/run_scoring_MBC_ULP.sh

Trim align MBC ULP
	PATH="$PATH:/fh/fast/ha_g/user/adoebley/projects/griffin_paper/Ulz_nucleosome_profiling/downloads"
	launchSlurm -p campus-new -b sh_scripts/trim_align_MBC_ULP.sh
	sbatch --array=1-181 trim_align_MBC_ULP.slurm
	sacct -j 27191212

run ichorCNA
	Ran ichorCNA (50mb bins, hg19) here: 
	/fh/fast/ha_g/user/adoebley/projects/griffin_paper/ichorCNA/ichorCNA_2020-11-20/scripts/MBC_ULP_50kb_for_Ulz_pipeline
	The snakemake failed out because there isn't a panel of normals for 50mb bins (I think I had previously figured out that this was the cause of the failure), however all the relevant outputs were created before that point
	I modified the ichorCNA snakefile to ignore the missing files

copy over ichor results and reformat them
	python write_scripts/reformat_log2_file.py -i /fh/fast/ha_g/user/adoebley/projects/griffin_paper/ichorCNA/ichorCNA_2020-11-20/scripts/MBC_ULP_50kb_for_Ulz_pipeline/results/ichorCNA/

Run analysis
	launchSlurm -p campus-new -b sh_scripts/run_tf_analysis_MBC_ULP.sh
	sbatch --array=1-181 run_tf_analysis_MBC_ULP.slurm
	sacct -j 27408032

Run scoring 
	launchSlurm -p campus-new -b sh_scripts/run_scoring_MBC_ULP.sh
	sbatch --array=1-181 run_scoring_MBC_ULP.slurm
	This script is hard coded to score all the different analysis types at the same time. It generates some empty files and error messages but continues to the actual scoring pipeline. 
	sacct -j 28215669


#RUNNING ADDTIONAL SAMPLES (lower tumor fraction and samples with new labels)

write scripts
	python write_scripts/write_scripts_new.py -i config/additional_MBC_hg19.yaml -n more_MBC_ULP

	outputs:
		Commands to run scripts:
		launchSlurm -p campus-new -b sh_scripts/trim_align_more_MBC_ULP.sh
		python write_scripts/reformat_log2_file.py -i <path/to/results/ichorCNA>
		launchSlurm -p campus-new -b sh_scripts/run_tf_analysis_more_MBC_ULP.sh
		launchSlurm -p campus-new -b sh_scripts/run_scoring_more_MBC_ULP.sh

Trim align MBC ULP
	PATH="$PATH:/fh/fast/ha_g/user/adoebley/projects/griffin_paper/Ulz_nucleosome_profiling/downloads"
	launchSlurm -p campus-new -b sh_scripts/trim_align_more_MBC_ULP.sh
	sbatch --array=1-149 trim_align_more_MBC_ULP.slurm
	sacct -j 30211410

run ichorCNA
	Ran ichorCNA (50mb bins, hg19) here: 
	/fh/fast/ha_g/user/adoebley/projects/griffin_paper/ichorCNA/ichorCNA_2020-11-20/scripts/MBC_ULP_50kb_for_Ulz_pipeline

copy over ichor results and reformat them
	python write_scripts/reformat_log2_file.py -i /fh/fast/ha_g/user/adoebley/projects/griffin_paper/ichorCNA/ichorCNA_2020-11-20/scripts/MBC_ULP_50kb_for_Ulz_pipeline/results/ichorCNA/

Run analysis
	launchSlurm -p campus-new -b sh_scripts/run_tf_analysis_more_MBC_ULP.sh
	sbatch --array=1-149 run_tf_analysis_more_MBC_ULP.slurm
	sacct -j 30272863

Run scoring
	launchSlurm -p campus-new -b sh_scripts/run_scoring_more_MBC_ULP.sh
	sbatch --array=1-149 run_scoring_more_MBC_ULP.slurm
	sacct -j 30288010