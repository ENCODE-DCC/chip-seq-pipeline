#!/usr/bin/env bash

# full IDR template
chip_workflow.py --target tf --debug --title IDR_Template --outp "E3 ChIP-seq"


## Complete experiments

# ECSR000EEB SE IDR
chip_workflow.py --target tf --debug --title ENCSR000EEB-fullIDR --outf /ENCSR000EEB-fullIDR-$(date +"%Y%m%d%H%M") --idr --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ECSR000EEB SE IDRv2
chip_workflow.py --target tf --debug --title ENCSR000EEB-fullIDR --outf /ENCSR000EEB-fullIDR-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ENCSR000EEB IDRv2 from TA
chip_workflow.py --target tf --debug --nomap --yes \
--title ENCSR000EEB-fullIDRnomap --outf /ENCSR000EEB-fullIDRnomap-$(date +"%Y%m%d%H%M")  \
--rep1pe false --rep2pe false \
--idr --idrversion 2 \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR000BUA SE IDRv2
chip_workflow.py --target tf --debug --yes \
--title ENCSR000BUA-fullIDR --outf /ENCSR000BUA-fullIDRv2-$(date +"%Y%m%d%H%M") \
--idr --idrversion 2 \
--rep1 /ENCSR000BUA/rep1/ENCFF000RBI.fastq.gz \
--rep2 /ENCSR000BUA/rep2/ENCFF000RBC.fastq.gz \
--ctl1 /ENCSR000BUA/ctl1/ENCFF000RCK.fastq.gz \
--ctl2 /ENCSR000BUA/ctl2/ENCFF000RCP.fastq.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR000BUA SE IDRv2 from TA
chip_workflow.py --target tf --debug --nomap --yes \
--title ENCSR000BUA-fullIDR --outf /ENCSR000BUA-fullIDRv2-$(date +"%Y%m%d%H%M") \
--rep1pe false --rep2pe false \
--idr --idrversion 2 \
--rep1 /ENCSR000BUA/rep1/ENCFF000RBI.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /ENCSR000BUA/rep2/ENCFF000RBC.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /ENCSR000BUA/ctl1/ENCFF000RCK.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /ENCSR000BUA/ctl2/ENCFF000RCP.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR000DMT SE IDRv2
chip_workflow.py --target tf --debug --title ENCSR000DMT-fullIDR --outf /ENCSR000DMT-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR000DMT/rep1/ENCFF000SBI.fastq.gz \
--rep2 /ENCSR000DMT/rep2/ENCFF000SBK.fastq.gz \
--ctl1 /ENCSR000DMT/ctl/ENCFF000SAZ.fastq.gz \
--ctl2 /ENCSR000DMT/ctl/ENCFF000SAZ.fastq.gz

# ENCSR000DMT SE IDRv2 from TA
chip_workflow.py --target tf --debug --nomap --yes \
--title ENCSR000DMT-fullIDR --outf /ENCSR000DMT-fullIDRv2-$(date +"%Y%m%d%H%M") \
--idr  --idrversion 2 \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000DMT/rep1/ENCFF000SBI.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /ENCSR000DMT/rep2/ENCFF000SBK.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /ENCSR000DMT/ctl/ENCFF000SAZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /ENCSR000DMT/ctl/ENCFF000SAZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR769ZTN PE IDRv2
chip_workflow.py --target tf --debug --title ENCSR769ZTN-fullIDR --outf /ENCSR769ZTN-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR769ZTN/rep1/ENCFF002ELM.fastq.gz /ENCSR769ZTN/rep1/ENCFF002ELL.fastq.gz \
--rep2 /ENCSR769ZTN/rep2/ENCFF002ELK.fastq.gz /ENCSR769ZTN/rep2/ENCFF002ELJ.fastq.gz \
--ctl1 /ENCSR769ZTN/ctl1/ENCFF002EFT.fastq.gz /ENCSR769ZTN/ctl1/ENCFF002EFQ.fastq.gz \
--ctl2 /ENCSR769ZTN/ctl2/ENCFF002EFU.fastq.gz /ENCSR769ZTN/ctl2/ENCFF002EFS.fastq.gz

# ENCSR769ZTN PE IDRv2 from TA
chip_workflow.py --target tf --debug --nomap --yes \
--title ENCSR769ZTN-fullIDR --outf /ENCSR769ZTN-fullIDRv2-$(date +"%Y%m%d%H%M") \
--idr --idrversion 2 \
--rep1pe true --rep2pe true \
--rep1 /ENCSR769ZTN/rep1/ENCFF002ELMENCFF002ELL.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--rep2 /ENCSR769ZTN/rep2/ENCFF002ELKENCFF002ELJ.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--ctl1 /ENCSR769ZTN/ctl1/ENCFF002EFTENCFF002EFQ.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--ctl2 /ENCSR769ZTN/ctl2/ENCFF002EFUENCFF002EFS.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR795HTY PE IDRv2 map and call peaks
chip_workflow.py --target tf --debug --yes \
--title ENCSR795HTY-map-idr2-$(date +"%Y%m%d%H%M") --outf /ENCSR795HTY-map-idr2-$(date +"%Y%m%d%H%M") \
--idr --idrversion 2 \
--rep1 /ENCSR795HTY/rep1/ENCFF240UHP.fastq.gz /ENCSR795HTY/rep1/ENCFF555NNB.fastq.gz \
--rep2 /ENCSR795HTY/rep2/ENCFF859GVE.fastq.gz /ENCSR795HTY/rep2/ENCFF974VFA.fastq.gz \
--ctl1 /ENCSR795HTY/ctl1/ENCFF445UEI.fastq.gz /ENCSR795HTY/ctl1/ENCFF141YIN.fastq.gz \
--ctl2 /ENCSR795HTY/ctl2/ENCFF210VLR.fastq.gz /ENCSR795HTY/ctl2/ENCFF057GZE.fastq.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR795HTY PE IDRv2 from TA
chip_workflow.py --target tf --debug --nomap --yes \
--title ENCSR795HTY-nomap-idr2 --outf /ENCSR795HTY-nomap-idr2-$(date +"%Y%m%d%H%M") \
--idr --idrversion 2 \
--rep1pe true --rep2pe true \
--rep1 file-BZ05pb80zp5PFPKxXxXV5ZV8 \
--rep2 file-BZ0KG380YBv2BV14PxXqYfBJ \
--ctl1 file-BZ0pPy00ByzPFPKxXxXV8BJ4 \
--ctl2 file-BZ1Bp0j0xByP4350J52qPb7Z \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"


## chr1 extracts

# ENCSR000EEB chr1 IDR from TA
chip_workflow.py --target tf --debug --title ENCSR000EEB-fullIDRtachr1 --outf /ENCSR000EEB-fullIDRtachr1 --idr --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr1.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr1.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr1.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr1.tagAlign.gz


## chr21 extracts

# ECSR000EEB chr21 SE IDRv2
chip_workflow.py --target tf --debug --title ENCSR000EEBchr21-fullIDR --outf /ENCSR000EEBchr21-fullIDR-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /test_data/ENCFF000XUL.chr21.fq.gz \
--rep2 /test_data/ENCFF000XUK.chr21.fq.gz \
--ctl1 /test_data/ENCFF000XTF.chr21.fq.gz \
--ctl2 /test_data/ENCFF000XTF.chr21.fq.gz

# ENCSR000EEB chr21 IDR from TA
chip_workflow.py --target tf --debug --title ENCSR000EEB-fullIDRtachr21 --outf /ENCSR000EEB-fullIDRtachr21-$(date +"%Y%m%d%H%M") --idr --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr21.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz

# ENCSR000EEB chr21 IDRv2 from TA
chip_workflow.py --target tf --debug --yes --nomap \
--title ENCSR000EEB-fullIDRtachr21 --outf /ENCSR000EEB-fullIDRtachr21-$(date +"%Y%m%d%H%M") \
--idr --idrversion 2 \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr21.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--genomesize hs --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--reference "ENCODE Reference Files:/hg19/hg19_XY.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# ENCSR000EEB chr21 no idr from TA
chip_workflow.py --target tf --debug --title ENCSR000EEB-fullIDRtachr21 --outf /ENCSR000EEB-fullIDRtachr21-$(date +"%Y%m%d%H%M") --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr21.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz

## histones

# ENCSR678FIT chr19 IDRv2 from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-chr19-ta-IDR2 --outf /ENCSR678FIT-chr19-ta-IDR2-$(date +"%Y%m%d%H%M") --idr --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--rep2 /test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl1 /test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl2 /test_data/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT IDRv2 from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-ta-IDR2 --outf /ENCSR678FIT-ta-IDR2-$(date +"%Y%m%d%H%M") --idr --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep1/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep2/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT chr19 overlap only from TA
chip_workflow.py --target histone --debug --title ENCSR678FIT-chr19-ta-OL --outf /ENCSR678FIT-chr19-ta-OL-$(date +"%Y%m%d%H%M") --nomap --yes \
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
--input "control=/ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz" \
--input "experiment=/ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz" \
--input "xcor_scores_input=/test_data/ENCFF000XUL-chr21.tagAlign.sample.15.SE.tagAlign.gz.cc.qc" \
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

# Run peak overlap for ENCSR678FIT chr19
dx run \
--input "rep1_peaks=/ENCSR678FIT-chr19-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.narrowPeak.gz" \
--input "rep2_peaks=/ENCSR678FIT-chr19-ta-IDR2-201505041839/encode_macs2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.narrowPeak.gz" \
--input "pooled_peaks=/ENCSR678FIT-chr19-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19-ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19_pooled.tagAlign.narrowPeak.gz" \
--input "pooledpr1_peaks=/ENCSR678FIT-chr19-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.SE.pr1-ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.SE.pr1_pooled.tagAlign.narrowPeak.gz" \
--input "pooledpr2_peaks=/ENCSR678FIT-chr19-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz.SE.pr2-ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz.SE.pr2_pooled.tagAlign.narrowPeak.gz" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "peak_type=narrowPeak" \
--verbose \
--destination /overlap_peaks_test_np \
--name overlap_peaks_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
overlap_peaks

# Run peak overlap for ENCSR678FIT
dx run \
--input "rep1_peaks=/ENCSR678FIT-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.narrowPeak.gz" \
--input "rep2_peaks=/ENCSR678FIT-ta-IDR2-201505041839/encode_macs2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.narrowPeak.gz" \
--input "pooled_peaks=/ENCSR678FIT-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE-ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE_pooled.tagAlign.narrowPeak.gz" \
--input "pooledpr1_peaks=/ENCSR678FIT-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.SE.pr1-ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.SE.pr1_pooled.tagAlign.narrowPeak.gz" \
--input "pooledpr2_peaks=/ENCSR678FIT-ta-IDR2-201505041839/encode_macs2/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz.SE.pr2-ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz.SE.pr2_pooled.tagAlign.narrowPeak.gz" \
--input "chrom_sizes=ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--input "as_file=ENCODE Reference Files:/narrowPeak.as" \
--input "peak_type=narrowPeak" \
--verbose \
--destination /overlap_peaks_test_np \
--name overlap_peaks_test \
--delay-workspace-destruction \
--priority high \
--yes \
--watch \
overlap_peaks


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
--reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "ENCSR000EEB-hMAFK-chr21-$(date +"%Y%m%d%H%M")" \
--title "ENCSR000EEB-hMAFK-chr21-$(date +"%Y%m%d%H%M")" \
--rep1 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R1-ENCFF000XTT.chr21.fq.gz" \
--rep2 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R2-ENCFF000XTU.chr21.fq.gz" \
--ctl1 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/C1-ENCFF000XSJ.chr21.fq.gz" \
--yes

# Build full histone pipeline on ENCSR087PLZ chr21 extracts:
chip_workflow.py \
--target histone \
--chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/male.mm10.tar.gz" \
--outf "ENCSR087PLZ-mH3K9ac-chr19-$(date +"%Y%m%d%H%M")" \
--title "ENCSR087PLZ-mH3K9ac-chr19-$(date +"%Y%m%d%H%M")" \
--rep1 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/R1-ENCFF560GLI.chr19.fq.gz" \
--rep2 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/R2-ENCFF891NNX.chr19.fq.gz" \
--ctl1 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/C1-ENCFF069WCH.chr19.fq.gz" \
--ctl2 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/C2-ENCFF101KOM.chr19.fq.gz" \
--yes

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


#######
## deploy workflows to the ENCODE Universal Pipelines project
chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq (no reference)" \
--outf "/ChIP-seq/"

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (no reference)" \
--outf "/ChIP-seq/"

chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq (hg19)" \
--chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "/ChIP-seq/"

chip_workflow.py \
--target tf \
--name "ENCODE TF ChIP-seq Unary Control (hg19)" \
--chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \
--blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
--outf "/ChIP-seq/" \
--unary_control

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (mm10)" \
--chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \
--genomesize mm \
--reference "ENCODE Reference Files:/mm10/male.mm10.tar.gz" \
--outf "/ChIP-seq/"

chip_workflow.py \
--target histone \
--name "ENCODE Histone ChIP-seq (GRCh38)" \
--chrom_sizes "ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \
--genomesize hs \
--reference "ENCODE Reference Files:/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz" \
--outf "/ChIP-seq/"
