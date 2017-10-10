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

chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unary Control Unreplicated (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--simplicate_experiment \
--unary_control

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unary Control Unreplicated (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--simplicate_experiment \
--unary_control

##hg19
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unreplicated (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--simplicate_experiment

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unreplicated (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz" \
--simplicate_experiment

chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unary Control Unreplicated (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--simplicate_experiment \
--unary_control

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unary Control Unreplicated (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz" \
--simplicate_experiment \
--unary_control

## no reference pre-loaded
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unreplicated (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--use_existing_folders \
--simplicate_experiment

chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unary Control Unreplicated (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--unary_control \
--use_existing_folders \
--simplicate_experiment

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unreplicated (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--use_existing_folders \
--simplicate_experiment

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unary Control Unreplicated (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--unary_control \
--use_existing_folders \
--simplicate_experiment

##mm10
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unreplicated (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/blacklists/mm10.blacklist.bed.gz" \
--simplicate_experiment


chip_workflow.py \
--target tf \
--unary_control \
--name "ENCODE TF ChIP-seq Unary Control Unreplicated (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/blacklists/mm10.blacklist.bed.gz" \
--simplicate_experiment


chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unreplicated (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--simplicate_experiment

chip_workflow.py \
--target histone \
--unary_control \
--name "ENCODE Histone ChIP-seq Unary Control Unreplicated (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--simplicate_experiment

