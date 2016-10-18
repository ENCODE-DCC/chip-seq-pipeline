#!/bin/bash

#######
## deploy workflows to the ENCODE Universal Pipelines project
## no reference pre-loaded
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--use_existing_folders

chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unary Control (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--unary_control \
--use_existing_folders

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--use_existing_folders

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq Unary Control (specify reference)" \
--outf "/ChIP-seq/workflows/" \
--unary_control \
--use_existing_folders

##hg19
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

chip_workflow.py \
--target tf \
--unary_control \
--name "ENCODE TF ChIP-seq Unary Control (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz"

chip_workflow.py \
--target histone \
--unary_control \
--name "ENCODE Histone ChIP-seq Unary Control (hg19)" \
--outf "/ChIP-seq/workflows/hg19/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/ChIP-seq/male.hg19.tar.gz"

##GRCh38
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--blacklist "ENCODE Reference Files:/GRCh38/blacklists/GRCh38.blacklist.bed.gz"

chip_workflow.py \
--target tf \
--unary_control \
--name "ENCODE TF ChIP-seq Unary Control (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--blacklist "ENCODE Reference Files:/GRCh38/blacklists/GRCh38.blacklist.bed.gz"

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz"

chip_workflow.py \
--target histone \
--unary_control \
--name "ENCODE Histone ChIP-seq Unary Control (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz"

##mm10
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--blacklist "ENCODE Reference Files:/mm10/blacklists/mm10.blacklist.bed.gz"

chip_workflow.py \
--target tf \
--unary_control \
--name "ENCODE TF ChIP-seq Unary Control (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--blacklist "ENCODE Reference Files:/mm10/blacklists/mm10.blacklist.bed.gz"

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz"

chip_workflow.py \
--target histone \
--unary_control \
--name "ENCODE Histone ChIP-seq Unary Control (mm10)" \
--outf "/ChIP-seq/workflows/mm10/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz"
