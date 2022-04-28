#get_GC_and_mappability_content.snakefile
#Anna-Lisa Doebley
#Template made 2021-12-12
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1 

#command to run snakemake (remove -np at end when done validating):
snakemake -s get_GC_and_mappability_content.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

"""

configfile: "config/sites.yaml"
configfile: "config/config.fred_hutch.yaml" 
configfile: "config/cluster_slurm.yaml" 

rule all:
	input: #commented out files are produced at the same time as other files
		expand("{results_dir}/{site_lists}.smoothed_GC_content.tsv",results_dir=config['results_dir'],site_lists=config['site_lists']),
		expand("{results_dir}/{site_lists}.mappability.tsv",results_dir=config['results_dir'],site_lists=config['site_lists'])


rule fetch_vals:
	input:
		site_list = lambda wildcards: config["site_lists"][wildcards.site_lists]
	output:
		GC = config['results_dir']+"/{site_lists}.smoothed_GC_content.tsv",
		mappability = config['results_dir']+"/{site_lists}.mappability.tsv",
		tmp_dir = temp(directory(config['tmp_dir']+"/{site_lists}"))
	params:
		GC_and_map_script = config['GC_and_map_script'],
		TF_name = "{site_lists}",

		results_dir=config['results_dir'],

		reference_genome = config['reference_genome'],
		mappability_bw = config['mappability_bw'],
		chrom_sizes_path = config['chrom_sizes_path'],

		chrom_column=config['chrom_column'],
		position_column=config['position_column'],
		chroms = config['chroms'],

		save_window = config['save_window'],

		encode_exclude = config['encode_exclude'],
		centromeres = config['centromeres'],
		gaps = config['gaps'],
		patches = config['patches'],
		alternative_haplotypes = config['alternative_haplotypes'],

	shell:
		"time {params.GC_and_map_script} \
		--TF_name {params.TF_name} \
		--site_list {input.site_list} \
		--tmp_dir {output.tmp_dir} \
		--results_dir {params.results_dir} \
		--reference_genome {params.reference_genome} \
		--mappability_bw {params.mappability_bw} \
		--chrom_sizes_path {params.chrom_sizes_path} \
		--chrom_column {params.chrom_column} \
		--position_column {params.position_column} \
		--chroms {params.chroms} \
		--save_window {params.save_window} \
		--exclude_paths {params.encode_exclude} {params.centromeres} {params.gaps} {params.patches} {params.alternative_haplotypes} "
