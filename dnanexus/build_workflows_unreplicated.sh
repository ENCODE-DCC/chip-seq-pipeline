#!/bin/bash

#######
## deploy workflows for unreplicated experiments
##GRCh38
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unreplicated (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--simplicate_experiment

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unreplicated (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--simplicate_experiment
