k100.Umap.MultiTrackMappability.hg38.bw
	wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw
	mv k100.Umap.MultiTrackMappability.bw k100.Umap.MultiTrackMappability.hg38.bw

	Hg38 has Umap single-read and multi-read tracks for different read lengths
	Info here: https://bismap.hoffmanlab.org/
	"multi-read mappability, the probability that a randomly selected read of length k in a given region is uniquely mappable"

k100.Umap.MultiTrackMappability.hg38.bedGraph
	#converted to bed graph
	bigWigToBedGraph k100.Umap.MultiTrackMappability.hg38.bw k100.Umap.MultiTrackMappability.hg38.bedGraph

#################
encode_unified_GRCh38_exclusion_list.bed
	wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
	gunzip ENCFF356LFX.bed.gz
	mv ENCFF356LFX.bed encode_unified_GRCh38_exclusion_list.bed
	Details here: https://www.encodeproject.org/files/ENCFF356LFX/

#################
https://genome.ucsc.edu/cgi-bin/hgTables

hg38_centromeres.bed
hg38_gaps.bed
hg38_fix_patches.bed
hg39_alternative_haplotypes.bed
	Downloaded from UCSC table browser
	See screenshots for settings


#################
k100_minus_exclusion_lists.mappable_regions.hg38.bed
	#a filter containing all fully mapable, non-excluded regions for the GC bias correction
	#created with:
	scripts/get_non_excluded_mappable_regions.ipynb
