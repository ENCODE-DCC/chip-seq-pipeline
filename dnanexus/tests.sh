#!/usr/bin/env bash

## DCC-spinner tests
# Map ECSR000EEB SE TF
map_only ENCSR000EEB \
--debug \
--assembly GRCh38 \
--outf /test_runs/ENCSR000EEB-map_only-$(date +"%Y%m%d%H%M") \
--key www \
--yes

# Call peaks on ECSR000EEB SE TF
call_chip_from_tas ENCSR000EEB \
--debug \
--assembly GRCh38 \
--inf "ENCODE - ChIP Production:/mapping_GRCh38/" \
--outf /test_runs/ENCSR000EEB-call_chip_from_tas-$(date +"%Y%m%d%H%M") \
--key www \
--yes

## full templates
chip_workflow.py --target tf --debug --title TF_template --outp "E3 ChIP-seq" --use_existing_folders
chip_workflow.py --target histone --debug --title histone_template --outp "E3 ChIP-seq" --use_existing_folders

## Complete experiments

# ECSR000EEB SE TF
chip_workflow.py \
--debug \
--target tf \
--title ENCSR000EEB-hs-MAFK \
--outf /test_runs/ENCSR000EEB-$(date +"%Y%m%d%H%M") \
--rep1 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep1-ENCFF000XUL.fastq.gz" \
--rep2 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep2-ENCFF000XUK.fastq.gz" \
--ctl1 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/ctl1-ENCFF000XTF.fastq.gz" \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--genomesize hs \
--yes

# ENCSR000EEB SE TF from TA
chip_workflow.py \
--debug \
--target tf \
--nomap \
--title ENCSR000EEB-hs-MAFK-nomap \
--outf /test_runs/ENCSR000EEB-nomap-$(date +"%Y%m%d%H%M")  \
--rep1pe false --rep2pe false \
--rep1 "ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR000EEB/rep1/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--rep2 "ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR000EEB/rep2/ENCFF000XUK.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--ctl1 "ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR000EEN/rep1/ENCFF000XTF.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--genomesize hs \
--yes

## chr extracts

# ECSR000EEB SE TF chr21
chip_workflow.py \
--debug \
--target tf \
--title ENCSR000EEB-hs-MAFK-chr21 \
--outf /test_runs/ENCSR000EEB-chr21-$(date +"%Y%m%d%H%M") \
--rep1 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep1-chr21.fq.gz" \
--rep2 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep2-chr21.fq.gz" \
--ctl1 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/ctl1-chr21.fq.gz" \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--genomesize hs \
--yes

# ENCSR000EEB SE TF from TA chr21
chip_workflow.py \
--debug \
--target tf \
--title ENCSR000EEB-hs-MAFK-nomap-chr21 \
--nomap \
--outf /test_runs/ENCSR000EEB-nomap-chr21-$(date +"%Y%m%d%H%M") \
--rep1pe false --rep2pe false \
--rep1 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep1-chr21.tagAlign.gz" \
--rep2 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep2-chr21.tagAlign.gz" \
--ctl1 "E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/ctl1-chr21.tagAlign.gz" \
--reference "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--chrom_sizes "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GRCh38_EBV.chrom.sizes" \
--blacklist "ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--genomesize hs \
--yes

## histones

# ENCSR678FIT chr19  from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-chr19-ta-IDR2 --outf /ENCSR678FIT-chr19-ta-IDR2-$(date +"%Y%m%d%H%M") --idr --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--rep2 /test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl1 /test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl2 /test_data/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT  from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-ta-IDR2 --outf /ENCSR678FIT-ta-IDR2-$(date +"%Y%m%d%H%M") --idr --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep1/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep2/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT chr19 overlap only from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-chr19-ta-OL --outf /ENCSR678FIT-chr19-ta-OL-$(date +"%Y%m%d%H%M") --nomap \
--rep1pe false --rep2pe false \
--rep1 /test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--rep2 /test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl1 /test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl2 /test_data/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT overlap only from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-ta-OL --outf /ENCSR678FIT-ta-OL-$(date +"%Y%m%d%H%M") --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep1/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep2/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR


## single applets

# encode_macs2 simplicate
dx run \
--input "rep1_ta=/test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "ctl1_ta=/test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "rep1_xcor=/test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "rep1_paired_end=false" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--input "genomesize=mm" \
--input "narrowpeak_as=ENCODE Reference Files:/narrowPeak.as" \
--input "gappedpeak_as=ENCODE Reference Files:/gappedPeak.as" \
--input "broadpeak_as=ENCODE Reference Files:/broadPeak.as" \
--input "prefix=test" \
--verbose \
--destination /encode_macs2_test/test_simplicate_$(date +"%Y%m%d%H%M") \
--name encode_macs2_test \
--delay-workspace-destruction \
--priority high \
--yes \
applets/encode_macs2

# encode_macs2 replicate
dx run \
--input "rep1_ta=/test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "rep2_ta=/test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "ctl1_ta=/test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "ctl2_ta=/test_data/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "rep1_xcor=/test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "rep2_xcor=/test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "rep1_paired_end=false" \
--input "rep2_paired_end=false" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--input "genomesize=mm" \
--input "narrowpeak_as=ENCODE Reference Files:/narrowPeak.as" \
--input "gappedpeak_as=ENCODE Reference Files:/gappedPeak.as" \
--input "broadpeak_as=ENCODE Reference Files:/broadPeak.as" \
--input "prefix=test" \
--verbose \
--destination /encode_macs2_test/test_replicate_$(date +"%Y%m%d%H%M") \
--name encode_macs2_test \
--delay-workspace-destruction \
--priority high \
--yes \
applets/encode_macs2

# encode_spp simplicate
dx run \
--input "rep1_ta=/test_data/ENCFF000XUL-chr21.tagAlign.gz" \
--input "ctl1_ta=/test_data/ENCFF000XTF-chr21.tagAlign.gz" \
--input "rep1_xcor=/test_data/ENCFF000XUL-chr21.tagAlign.sample.15.SE.tagAlign.gz.cc.qc" \

--input "rep1_paired_end=false" \
--input "chrom_sizes=ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "idr_peaks=true" \
--verbose \
--destination /encode_spp_test/test_simplicate_$(date +"%Y%m%d%H%M") \
--name encode_spp_test \
--delay-workspace-destruction \
--priority high \
--yes \
applets/encode_spp

# encode_spp replicate
dx run \
--input "rep1_ta=/test_data/ENCFF000XUL-chr21.tagAlign.gz" \
--input "rep2_ta=/test_data/ENCFF000XUK-chr21.tagAlign.gz" \
--input "ctl1_ta=/test_data/ENCFF000XTF-chr21.tagAlign.gz" \
--input "ctl2_ta=/test_data/ENCFF000XTF-chr21.tagAlign.gz" \
--input "rep1_xcor=/test_data/ENCFF000XUL-chr21.tagAlign.sample.15.SE.tagAlign.gz.cc.qc" \
--input "rep2_xcor=/test_data/ENCFF000XUK.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "rep1_paired_end=false" \
--input "rep2_paired_end=false" \
--input "chrom_sizes=ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "idr_peaks=true" \
--verbose \
--destination /encode_spp_test/test_replicate_$(date +"%Y%m%d%H%M") \
--name encode_spp_test \
--delay-workspace-destruction \
--priority high \
--yes \
applets/encode_spp

# idr2 replicate
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL-chr21.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK-chr21.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL-chr21-ENCFF000XUK-chr21_pooled.regionPeak.gz" \
--verbose \
--destination /IDRv2_test_replicate/test_$(date +"%Y%m%d%H%M") \
--name IDRv2_test_replicate \
--delay-workspace-destruction \
--priority high \
--yes \
applets/idr2

# idr2 simplicate
dx run \
--input "rep1_peaks=/test_data/simplicate/ENCFF000XUL-chr21_PR1.regionPeak.gz" \
--input "rep2_peaks=/test_data/simplicate/ENCFF000XUL-chr21_PR2.regionPeak.gz" \
--input "pooled_peaks=/test_data/simplicate/ENCFF000XUL-chr21.regionPeak.gz" \
--verbose \
--destination /IDRv2_test_simplicate/test_$(date +"%Y%m%d%H%M") \
--name IDRv2_test_simplicate \
--delay-workspace-destruction \
--priority high \
--yes \
applets/idr2

# encode_idr replicate
dx run \
--input "experiment=ENCSR936XTK" \
--input "r1pr_peaks=file-F09ZPy80xZ2gB29z41B97421" \
--input "rep1_ta=/test_data/ENCSR936XTK/bams/rep2/ENCFF960TNPENCFF640CBP.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz" \
--input "rep1_xcor=/test_data/ENCSR936XTK/bams/rep2/ENCFF960TNPENCFF640CBP.raw.srt.filt.srt.nodup.filt.nodup.sample.15.MATE1.tagAlign.gz.cc.qc" \
--input "rep1_signal=/test_data/ENCSR936XTK/peaks/encode_macs2/r1.pvalue_signal.bw" \
--input "r2pr_peaks=/test_data/ENCSR936XTK/peaks/idr2/ENCFF246DIPvENCFF246DIP.IDRv2.IDR0.05.narrowPeak.gz" \
--input "rep2_ta=/test_data/ENCSR936XTK/bams/rep3/ENCFF246DIPENCFF616WSS.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz" \
--input "rep2_xcor=/test_data/ENCSR936XTK/bams/rep3/ENCFF246DIPENCFF616WSS.raw.srt.filt.srt.nodup.filt.nodup.sample.15.MATE1.tagAlign.gz.cc.qc" \
--input "rep2_signal=/test_data/ENCSR936XTK/peaks/encode_macs2/r2.pvalue_signal.bw" \
--input "reps_peaks=/test_data/ENCSR936XTK/peaks/idr2/ENCFF960TNPvENCFF246DIP.IDRv2.IDR0.05.narrowPeak.gz" \
--input "pooled_signal=/test_data/ENCSR936XTK/peaks/encode_macs2/pool.pvalue_signal.bw" \
--input "pooledpr_peaks=file-F09ZQGQ0vPXkpF9yx59g566V" \
--input "paired_end=true" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "blacklist=ENCODE Reference Files:/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--verbose \
--destination /test_results/encode_idr/$(date +"%Y%m%d%H%M") \
--name encode_idr_test \
--delay-workspace-destruction \
--priority high \
--yes \
applets/encode_idr

# encode_idr simplicate
dx run \
--input "experiment=ENCSR218GSN" \
--input "r1pr_peaks=/test_data/ENCSR218GSN/idr2/R1PR1.fixcovR1PR2.fixco.IDRv2.IDR0.05.narrowPeak.gz" \
--input "rep1_ta=/test_data/ENCSR218GSN/bams/rep1/ENCFF664IIEENCFF292UBY.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz" \
--input "rep1_xcor=/test_data/ENCSR218GSN/bams/rep1/ENCFF664IIEENCFF292UBY.raw.srt.filt.srt.nodup.filt.nodup.sample.15.MATE1.tagAlign.gz.cc.qc" \
--input "rep1_signal=/test_data/ENCSR218GSN/encode_macs2/r1.pvalue_signal.bw" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "blacklist=ENCODE Reference Files:/GRCh38/blacklists/GRCh38.blacklist.bed.gz" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--verbose \
--destination /test_results/encode_idr/$(date +"%Y%m%d%H%M") \
--name encode_idr_test_singlicate \
--delay-workspace-destruction \
--priority high \
--yes \
applets/encode_idr

# Run ENCODE_map on 36 bp SE chr21 extract
dx run \
--input "reads1=ENCODE Uniform Processing Pipelines:/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R1-ENCFF000XTT.chr21.fq.gz" \
--input "reference_tar=ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--verbose \
--destination /encode_map_test/test_$(date +"%Y%m%d%H%M") \
--name encode_map_test \
--delay-workspace-destruction \
--priority high \
--yes \
encode_map


# Run ENCODE_map on 36 bp SE chr21 extract, crop to 25 bp
dx run \
--input "reads1=ENCODE Uniform Processing Pipelines:/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R1-ENCFF000XTT.chr21.fq.gz" \
--input "reference_tar=ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--input "crop_length=25" \
--verbose \
--destination /encode_map_test/test_$(date +"%Y%m%d%H%M") \
--name encode_map_test \
--delay-workspace-destruction \
--priority high \
--yes \
encode_map


# Run ENCODE_map on 100 bp PE
dx run \
--input "reads1=/test_data/TF_PE/ENCFF109UIV.fastq.gz" \
--input "reads2=/test_data/TF_PE/ENCFF748SHJ.fastq.gz" \
--input "reference_tar=ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--verbose \
--destination /encode_map_test/test_$(date +"%Y%m%d%H%M") \
--name encode_map_test \
--delay-workspace-destruction \
--priority high \
--yes \
encode_map


# Run ENCODE_map on 100 bp PE crop to 36
dx run \
--input "reads1=/test_data/TF_PE/ENCFF109UIV.fastq.gz" \
--input "reads2=/test_data/TF_PE/ENCFF748SHJ.fastq.gz" \
--input "reference_tar=ENCODE Uniform Processing Pipelines:/Reference Files/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--input "crop_length=36" \
--verbose \
--destination /encode_map_test/test_$(date +"%Y%m%d%H%M") \
--name encode_map_test \
--delay-workspace-destruction \
--priority high \
--yes \
encode_map


# Run IDRv1 applet with ENCSR000EEB chr21 rep1 vs rep2
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL-chr21.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK-chr21.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL-chr21-ENCFF000XUK-chr21_pooled.regionPeak.gz" \
--verbose \
--destination /IDRv1_test \
--name IDRv1_test \
--delay-workspace-destruction \
--priority high \
--yes \
idr

# Run IDRv1 applet with ENCSR000EEB rep1 vs rep2
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.fixcoord.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK.raw.srt.filt.nodup.srt.SE.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE-ENCFF000XUK.raw.srt.filt.nodup.srt.SE_pooled.regionPeak.gz" \
--verbose \
--destination /IDRv1_test \
--name IDRv1_test \
--delay-workspace-destruction \
--priority high \
--yes \
idr

# Run IDRv2 applet with ENCSR000EEB chr21 rep1 vs rep2
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL-chr21.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK-chr21.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL-chr21-ENCFF000XUK-chr21_pooled.regionPeak.gz" \
--verbose \
--destination /IDRv2_test \
--name IDRv2_test \
--delay-workspace-destruction \
--priority high \
--yes \
idr2

# Run IDRv2 applet with ENCSR000EEB rep1 vs rep2
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.fixcoord.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK.raw.srt.filt.nodup.srt.SE.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE-ENCFF000XUK.raw.srt.filt.nodup.srt.SE_pooled.regionPeak.gz" \
--verbose \
--destination /IDRv2_test \
--name IDRv2_test \
--delay-workspace-destruction \
--priority high \
--yes \
idr2

# Run SPP with ENCSR000EEB chr21 rep1 vs input
dx run \
--input "control=E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/ctl1-chr21.tagAlign.gz" \
--input "experiment=E3 ChIP-seq:/test_data/TF/ENCSR000EEB-hs-MAFK/rep1-chr21.tagAlign.gz" \
--input "xcor_scores_input=E3 ChIP-seq:/reference_analyses/hg19/ENCSR000EEB-chr21-201706221826/encode_map/rep1-chr21.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--verbose \
--destination /spp_test \
--name spp_test \
--delay-workspace-destruction \
--priority high \
--yes \
/ChIP-seq/applets/spp

# Run SPP with ENCSR000EEB rep1 vs input - bedtobigbed fails with code 255 with coodinates in scientific notation
dx run \
--input "control=/test_data/ENCFF000XTF.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "experiment=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "xcor_scores_input=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--verbose \
--destination /spp_test \
--name spp_test \
--delay-workspace-destruction \
--priority high \
--yes \
spp

# Test SPP with full GRCh38 vs GRCh38-minimalized tagAligns
dx run \
--input "control=ENCODE - ChIP Production:/test_TF_38/ENCFF972NVS.raw.srt.filt.nodup.srt.SE.minimal.tagAlign.gz" \
--input "experiment=ENCODE - ChIP Production:/test_TF_38/ENCFF002AXG.raw.srt.filt.nodup.srt.SE.minimal.tagAlign.gz" \
--input "xcor_scores_input=ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR464DKE/rep1/ENCFF002AXG.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input 'spp_version=1.14' \
--verbose \
--destination /spp_test \
--name spp_test_minimal \
--delay-workspace-destruction \
--priority high \
--yes \
/applets/spp

dx run \
--input "control=ENCODE - ChIP Production:/test_TF_38/ENCFF972NVS.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "experiment=ENCODE - ChIP Production:/test_TF_38/ENCFF002AXG.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "xcor_scores_input=ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR464DKE/rep1/ENCFF002AXG.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input 'spp_version=1.14' \
--verbose \
--destination /spp_test \
--name spp_test_full \
--delay-workspace-destruction \
--priority high \
--yes \
/applets/spp

dx run \
--input "control=ENCODE - ChIP Production:/test_TF_38/c1_22.tagAlign.gz" \
--input "experiment=ENCODE - ChIP Production:/test_TF_38/r1_22.tagAlign.gz" \
--input "xcor_scores_input=ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR464DKE/rep1/ENCFF002AXG.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input 'spp_version=1.14' \
--verbose \
--destination /spp_test \
--name spp_test_22 \
--delay-workspace-destruction \
--priority high \
--yes \
/applets/spp

dx run \
--input "control=ENCODE - ChIP Production:/test_TF_38/c1_22_Un.tagAlign.gz" \
--input "experiment=ENCODE - ChIP Production:/test_TF_38/r1_22_Un.tagAlign.gz" \
--input "xcor_scores_input=ENCODE - ChIP Production:/mapping_GRCh38/bams/ENCSR464DKE/rep1/ENCFF002AXG.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--verbose \
--destination /spp_test \
--name spp_test_22_Un \
--delay-workspace-destruction \
--priority high \
--yes \
/applets/spp

dx run \
--input "control=/test_TF_38/c1_22_all.tagAlign.gz" \
--input "experiment=/test_TF_38/r1_22_all.tagAlign.gz" \
--input "xcor_scores_input=/mapping_GRCh38/bams/ENCSR464DKE/rep1/ENCFF002AXG.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "bigbed=true" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--verbose \
--destination /spp_test \
--name spp_test_22_all \
--delay-workspace-destruction \
--priority high \
--yes \
--instance-type mem3_ssd1_x16 \
--debug-on All \
--ssh \
/applets/spp


# Run MACS2 with ENCSR678FIT chr19
dx run \
--input "experiment=/test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "xcor_scores_input=/test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "control=/test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--input "narrowpeak_as=ENCODE Reference Files:/narrowPeak.as" \
--input "gappedpeak_as=ENCODE Reference Files:/gappedPeak.as" \
--input "broadpeak_as=ENCODE Reference Files:/broadPeak.as" \
--input "genomesize=mm" \
--verbose \
--destination /macs2_test \
--name macs2_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
macs2

# Run MACS2 with ENCSR678FIT (narrow mark)
dx run \
--input "experiment=mm10_mapping:/e115/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "xcor_scores_input=mm10_mapping:/e115/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "control=mm10_mapping:/e115/bams/ENCSR817FFF/rep1/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--input "narrowpeak_as=ENCODE Reference Files:/narrowPeak.as" \
--input "gappedpeak_as=ENCODE Reference Files:/gappedPeak.as" \
--input "broadpeak_as=ENCODE Reference Files:/broadPeak.as" \
--input "genomesize=mm" \
--verbose \
--destination /macs2_test \
--name macs2_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
macs2

# Run MACS2 with ENCSR311TLE (broad mark)
dx run \
--input "experiment=/mm10_mapping/e115_50a/bams/ENCSR311TLE/rep1/ENCFF002BVQ-ENCFF002BVU_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "xcor_scores_input=/mm10_mapping/e115_50a/bams/ENCSR311TLE/rep1/ENCFF002BVQ-ENCFF002BVU_pooled.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc" \
--input "control=/mm10_mapping/e115_controls/bams/ENCSR081UQY/rep1/ENCFF001ZUV.raw.srt.filt.nodup.srt.SE.tagAlign.gz" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--input "narrowpeak_as=ENCODE Reference Files:/narrowPeak.as" \
--input "gappedpeak_as=ENCODE Reference Files:/gappedPeak.as" \
--input "broadpeak_as=ENCODE Reference Files:/broadPeak.as" \
--input "genomesize=mm" \
--verbose \
--destination /macs2_test_broad \
--name macs2_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
macs2

# Run pool for ENCSR936XTK rep1 and rep2 tas
dx run \
--input "inputs=/test_data/ENCSR936XTK/bams/rep2/ENCFF960TNPENCFF640CBP.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz" \
--input "inputs=/test_data/ENCSR936XTK/bams/rep3/ENCFF246DIPENCFF616WSS.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz" \
--verbose \
--destination /test_results/pool/$(date +"%Y%m%d%H%M") \
--name pool_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
/applets/pool

# Run xcor_only for ENCSR936XTK rep1_2 pool
dx run \
--input "input_tagAlign=/test_data/ENCSR936XTK/ENCFF960TNPENCFF640CBP.raw.srt.filt.srt.nodup.PE2SE-ENCFF246DIPENCFF616WSS.raw.srt.filt.srt.nodup.PE2SE_pooled.tagAlign.gz" \
--input "paired_end=true" \
--verbose \
--destination /test_results/xcor_only/$(date +"%Y%m%d%H%M") \
--name xcor_only_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
/applets/xcor_only

# Run peak overlap for ENCSR678FIT chr19
dx run \
--input "rep1_peaks=/test_data/ENCSR494LJG/peaks/encode_macs2/r1.narrowPeak.gz" \
--input "rep2_peaks=/test_data/ENCSR494LJG/peaks/encode_macs2/r2.narrowPeak.gz" \
--input "pooled_peaks=/test_data/ENCSR494LJG/peaks/encode_macs2/pool.narrowPeak.gz" \
--input "pooledpr1_peaks=/test_data/ENCSR494LJG/peaks/encode_macs2/ppr1.narrowPeak.gz" \
--input "pooledpr2_peaks=/test_data/ENCSR494LJG/peaks/encode_macs2/ppr2.narrowPeak.gz" \
--input "chrom_sizes=ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "peak_type=narrowPeak" \
--input "pool_ta=/test_data/ENCSR494LJG/ENCFF480CWX-ENCFF074KIQ-ENCFF002BBF_pooled-crop.raw.srt.filt.nodup.srt.SE-ENCFF002AUU.raw.srt.filt.nodup.srt.SE_pooled.tagAlign.gz" \
--input "pool_xcor=/test_data/ENCSR678FIT-chr19_rep1_2_pooled.tagAlign.sample.15.SE.tagAlign.gz.cc.qc" \
--verbose \
--destination /test_results/overlap_peaks/$(date +"%Y%m%d%H%M") \
--name overlap_peaks_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
/applets/overlap_peaks


dx run \
--input "paired_end=false" \
--input "input_bam=/test_data/ENCFF000XUL.chr21.raw.srt.bam" \
--verbose \
--destination /test_output/ \
--name filter_qc_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
/applets/filter_qc


###### new tests for ChIP Demo

# Build full TF pipeline on ENCSR464DKE chr21 extracts:
chip_workflow.py \
--target tf \
--chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "ENCSR464DKE-hCTCF-chr21-$(date +"%Y%m%d%H%M")" \
--title "ENCSR464DKE-hCTCF-chr21-$(date +"%Y%m%d%H%M")" \
--rep1 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/R1-ENCFF921SED.chr21.fq.gz" \
--rep2 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/R2-ENCFF812KOM.chr21.fq.gz" \
--ctl1 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/C1-ENCFF690VPV.chr21.fq.gz" \
--ctl2 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/C2-ENCFF357TLV.chr21.fq.gz" \
--yes

# Build full TF pipeline on ENCSR286PCG chr21 extracts:
chip_workflow.py \
--target tf \
--chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "ENCSR286PCG-hZBED1-chr21-$(date +"%Y%m%d%H%M")" \
--title "ENCSR286PCG-hZBED1-chr21-$(date +"%Y%m%d%H%M")" \
--rep1 "/ChIP-seq/test_data/ENCSR286PCG-hZBED1/R1-ENCFF016MFU.chr21.fq.gz" \
--rep2 "/ChIP-seq/test_data/ENCSR286PCG-hZBED1/R2-ENCFF986OUP.chr21.fq.gz" \
--ctl1 "/ChIP-seq/test_data/ENCSR286PCG-hZBED1/C1-ENCFF048VYQ.chr21.fq.gz" \
--ctl2 "/ChIP-seq/test_data/ENCSR286PCG-hZBED1/C2-ENCFF839YOM.chr21.fq.gz" \
--yes

# Build full TF pipeline on ENCSR000EEB chr21 extracts:
chip_workflow.py \
--target tf \
--chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/hg19/ChIP-seq/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "ENCSR000EEB-hMAFK-chr21-$(date +"%Y%m%d%H%M")" \
--title "ENCSR000EEB-hMAFK-chr21-$(date +"%Y%m%d%H%M")" \
--rep1 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R1-ENCFF000XTT.chr21.fq.gz" \
--rep2 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R2-ENCFF000XTU.chr21.fq.gz" \
--ctl1 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/C1-ENCFF000XSJ.chr21.fq.gz" \
--yes

# Build full TF pipeline on ENCSR000EEB chr21 extracts as simplicate:
chip_workflow.py \
--target tf \
--chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/hg19/ChIP-seq/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "ENCSR000EEB-hMAFK-chr21-simplicate-$(date +"%Y%m%d%H%M")" \
--title "ENCSR000EEB-hMAFK-chr21-simplicate-$(date +"%Y%m%d%H%M")" \
--rep1 "/test_data/ENCFF000XUL.chr21.fq.gz" \
--ctl1 "/test_data/ENCFF000XTF.chr21.fq.gz" \
--yes


# Build full histone pipeline on ENCSR087PLZ chr19 extracts:
chip_workflow.py \
--target histone \
--chrom_sizes "ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--outf "ENCSR087PLZ-mH3K9ac-chr19-$(date +"%Y%m%d%H%M")" \
--title "ENCSR087PLZ-mH3K9ac-chr19-$(date +"%Y%m%d%H%M")" \
--rep1 "/test_data/histone_demo/ENCSR087PLZ-rep1-mm-chr19.bam.fq.gz" \
--ctl1 "/test_data/histone_demo/ENCSR087PLZ-ctl1-mm-chr19.bam.fq.gz" \
--rep2 "/test_data/histone_demo/ENCSR087PLZ-rep2-mm-chr19.bam.fq.gz" \
--ctl2 "/test_data/histone_demo/ENCSR087PLZ-ctl2-mm-chr19.bam.fq.gz" \
--yes

# Build full histone pipeline on ENCSR087PLZ chr19 extracts as simplicate:
chip_workflow.py \
--target histone \
--chrom_sizes "ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--outf "ENCSR087PLZ-mH3K9ac-chr19-simplicate-$(date +"%Y%m%d%H%M")" \
--title "ENCSR087PLZ-mH3K9ac-chr19-simpliate-$(date +"%Y%m%d%H%M")" \
--rep1 "/test_data/histone_demo/ENCSR087PLZ-rep1-mm-chr19.bam.fq.gz" \
--ctl1 "/test_data/histone_demo/ENCSR087PLZ-ctl1-mm-chr19.bam.fq.gz" \
--yes

# Build full histone pipeline on ENCSR238SGC:
chip_workflow.py \
--target histone \
--chrom_sizes "ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--outf "ENCSR238SGC-H3K4me1-$(date +"%Y%m%d%H%M")" \
--title "ENCSR238SGC-H3K4me1-$(date +"%Y%m%d%H%M")" \
--rep1 "/test_data/histone_chip/ENCFF833BLU.fastq.gz" \
--rep2 "/test_data/histone_chip/ENCFF646LXU.fastq.gz" \
--ctl1 "/test_data/histone_chip/ENCFF524CAC.fastq.gz" \
--ctl2 "/test_data/histone_chip/ENCFF163AJI.fastq.gz" \
--yes

# Build full histone pipeline on ENCSR238SGC chr19 extracts:
chip_workflow.py \
--target histone \
--chrom_sizes "ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--outf "ENCSR238SGC-H3K4me1-chr19-$(date +"%Y%m%d%H%M")" \
--title "ENCSR238SGC-H3K4me1-chr19-$(date +"%Y%m%d%H%M")" \
--rep1 "/test_data/histone_chip/ENCSR238SGC-H3K4me1/ENCFF833BLU-chr19.fq.gz" \
--rep2 "/test_data/histone_chip/ENCSR238SGC-H3K4me1/ENCFF646LXU-chr19.fq.gz" \
--ctl1 "/test_data/histone_chip/ENCSR238SGC-H3K4me1/ENCFF524CAC-chr19.fq.gz" \
--ctl2 "/test_data/histone_chip/ENCSR238SGC-H3K4me1/ENCFF163AJI-chr19.fq.gz" \
--yes

# Build full histone pipeline on ENCSR238SGC chr19 using rep1 as simplicate:
chip_workflow.py \
--target histone \
--chrom_sizes "ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz" \
--outf "ENCSR238SGC-H3K4me1-simplicate-chr19-$(date +"%Y%m%d%H%M")" \
--title "ENCSR238SGC-H3K4me1-simplicate-chr19-$(date +"%Y%m%d%H%M")" \
--rep1 "/test_data/histone_chip/ENCSR238SGC-H3K4me1/ENCFF833BLU-chr19.fq.gz" \
--ctl1 "/test_data/histone_chip/ENCSR238SGC-H3K4me1/ENCFF524CAC-chr19.fq.gz" \
--yes
