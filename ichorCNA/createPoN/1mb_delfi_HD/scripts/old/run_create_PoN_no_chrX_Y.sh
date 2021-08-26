Rscript ../../mod_createPanelOfNormals.R \
	--gcWig /fh/fast/ha_g/user/adoebley/ichorCNA/ichorCNA_2020-11-20/inst/extdata/gc_hg38_1000kb.wig \
	--mapWig /fh/fast/ha_g/user/adoebley/ichorCNA/ichorCNA_2020-11-20/inst/extdata/map_hg38_1000kb.wig \
	--filelist 1mb_delfi_HD_wigs.txt --outfile ../results/HD_delfi_PoN_hg38_1Mb \
	--centromere /fh/fast/ha_g/user/adoebley/ichorCNA/ichorCNA_2020-11-20/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
	--chrs "paste0('chr', c(1:22, \"X\",\"Y\"))" \
	--genomeStyle UCSC \
	--chrNormalize "paste0('chr', c(1:22))"