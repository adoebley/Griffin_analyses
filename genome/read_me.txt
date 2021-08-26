##################
Hg38 has Umap single-read and multi-read tracks for different read lengths

commands to download tracks and convert to bed:

k50.Umap.MultiTrackMappability.hg38.bw
	wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw
	mv k50.Umap.MultiTrackMappability.bw k50.Umap.MultiTrackMappability.hg38.bw

k50.Umap.MultiTrackMappability.hg38.bedGraph
	#converted to bed
	bigWigToBedGraph k50.Umap.MultiTrackMappability.hg38.bw k50.Umap.MultiTrackMappability.hg38.bedGraph

##################
RepeatMasker.hg38.bed
	This is the repeatMasker track for hg38 downloaded from the table browser
	http://genome.ucsc.edu/cgi-bin/hgTables

	The file was downloaded according to the directions here (except hg38 rather than hg19):
	https://deeptools.readthedocs.io/en/latest/content/tools/computeGCBias.html

	The settings used are in download_repeat_masker.png screenshot (same as the settings described in the link above)
	No modifiations were made to the settings on the second download page so I did not screenshot these settings

#################
repeat_masker.mapable.k50.Umap.hg38
	#a filter containing all fully mapable, non repetitive regions
	#created with:
	scripts/create_repeat_mask_filter.ipynb
