#!/bin/bash

#######
## deploy workflows for mapping only (no peak-calling)
##GRCh38
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Mapping (GRCh38)" \
--outf "/ChIP-seq/workflows/GRCh38/" \
--use_existing_folders \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--maponly
