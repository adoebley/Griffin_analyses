#griffin_nucleosome_profiling.snakefile
#Anna-Lisa Doebley
#Template made 2021-04-12
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1 #this loads deeptools
PATH="$PATH:/fh/fast/ha_g/user/adoebley/projects/griffin_paper/scripts/"


#command to run snakemake (remove -np at end when done validating):
snakemake -s griffin_nucleosome_profiling.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

#command to run snakemake on restart (remove -np at end when done validating):
#need to add '-q restart-new' according to scicomp. This should be fixed at some point.
snakemake -s griffin_nucleosome_profiling.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p restart-new -q restart-new --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName} --requeue" -j 40 -np
"""

configfile: "config/samples.GC.yaml" #this is generated by griffin_GC_correction
configfile: "config/sites.yaml"
configfile: "config/config.yaml" 
configfile: "config/cluster_slurm.yaml" 

rule all:
	input: #commented out files are produced at the same time as other files
		expand("{results_dir}/coverage/all_sites/{samples}.all_sites.coverage.txt",results_dir=config['results_dir'],samples=config['samples']),
		#expand("{results_dir}/plots/{site_files}.nucleosome_center.pdf", results_dir = config['results_dir'], site_files = config['site_files'])

rule calc_cov:
	input:
		bam = lambda wildcards: config["samples"][wildcards.samples]['bam'],
		GC_bias = lambda wildcards: config["samples"][wildcards.samples]['GC_bias']
	output:
		cov_file = config['results_dir']+"/coverage/all_sites/{samples}.all_sites.coverage.txt"
	params:
		sample_name = "{samples}",
		background_normalization = config['background_normalization'],

		results_dir=config['results_dir'],
		sites_yaml = config['sites_yaml'],
		reference_genome = config['reference_genome'],

		chrom_col=config['chrom_column'],
		chroms = config['chroms'],
		norm_window = config['norm_window'],
		plot_window = config['plot_window'],
		fragment_length = config['fragment_length'],
		
		step=config['step'],
		size_range=config['size_range'],
		map_q=config['map_quality'],
		strand_col=config['strand_column'],

		individual=config['individual'],
		smoothing = config['smoothing'],

		number_of_sites=config['number_of_sites'],
		sort_by=config['sort_by'],
		ascending=config['ascending'],

		cpus = config['calc_cov']['ncpus']
		
	shell:
		"time griffin_calc_coverage.py --sample_name {params.sample_name} --bam {input.bam} --GC_bias {input.GC_bias} \
		--background_normalization {params.background_normalization} \
		--sites_yaml {params.sites_yaml} --reference_genome {params.reference_genome} --results_dir {params.results_dir} \
		--chrom_column {params.chrom_col} --chroms {params.chroms} --norm_window {params.norm_window} \
		--plot_window {params.plot_window} --fragment_length {params.fragment_length} \
		--step {params.step} --size_range {params.size_range} --map_quality {params.map_q} \
		--strand_column {params.strand_col} --individual {params.individual} --smoothing {params.smoothing} \
		--num_sites {params.number_of_sites} --sort_by {params.sort_by} --ascending {params.ascending} \
		--cpu {params.cpus} "

# rule generate_plots:
# 	input:
# 		cov_files = expand("{results_dir}/coverage/all_sites/{samples}.all_sites.coverage.txt", results_dir = config['results_dir'], samples = config['samples'])
# 	output:
# 		plots = expand("{results_dir}/plots/{site_files}.nucleosome_center.pdf", results_dir = config['results_dir'], site_files = config['site_files'])
# 	params:
# 		plot_values=config['plot_window_values'],
# 		fragment_length=config['fragment_length'],
# 		step=config['step'],
# 		individual=config['individual'],
# 		results_dir=config['results_dir']
# 	shell:
# 		'simple_plotter.py --in_files {input.cov_files} --plot_window {params.plot_values} \
# 		--fragment_length {params.fragment_length} --step {params.step} \
# 		--individual {params.individual} --out_dir {params.results_dir}'

