#!/usr/bin/env bash

# full IDR template
tf_workflow.py --debug --name IDR_Template --outp "E3 ChIP-seq"


## Complete experiments

# ECSR000EEB SE full IDR
tf_workflow.py --debug --name ENCSR000EEB-fullIDR --outf /ENCSR000EEB-fullIDR-$(date +"%Y%m%d%H%M") --idr --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ECSR000EEB SE full IDRv2
tf_workflow.py --debug --name ENCSR000EEB-fullIDR --outf /ENCSR000EEB-fullIDR-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ENCSR000EEB SE full IDR only (no naive peaks)
tf_workflow.py --debug --name ENCSR000EEB-fullIDRonly --outf /ENCSR000EEB-fullIDRonly-$(date +"%Y%m%d%H%M") --idronly --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ENCSR000EEB IDRv2 only from TA
tf_workflow.py --debug --name ENCSR000EEB-fullIDRnomap --outf /ENCSR000EEB-fullIDRnomap-$(date +"%Y%m%d%H%M") --idronly --nomap --yes \
--rep1pe false --rep2pe false \
--idrversion 2 \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.raw.srt.filt.nodup.srt.SE.tagAlign.gz

# ENCSR000BUA SE full IDRv2
tf_workflow.py --debug --name ENCSR000BUA-fullIDR --outf /ENCSR000BUA-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR000BUA/rep1/ENCFF000RBI.fastq.gz \
--rep2 /ENCSR000BUA/rep2/ENCFF000RBC.fastq.gz \
--ctl1 /ENCSR000BUA/ctl1/ENCFF000RCK.fastq.gz \
--ctl2 /ENCSR000BUA/ctl2/ENCFF000RCP.fastq.gz

# ENCSR000BUA SE full IDRv2 nomap
tf_workflow.py --debug --name ENCSR000BUA-fullIDR --outf /ENCSR000BUA-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --nomap --yes \
--rep1pe false --rep2pe false \
--idrversion 2 \
--rep1 /ENCSR000BUA/rep1/ENCFF000RBI.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /ENCSR000BUA/rep2/ENCFF000RBC.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /ENCSR000BUA/ctl1/ENCFF000RCK.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /ENCSR000BUA/ctl2/ENCFF000RCP.raw.srt.filt.nodup.srt.SE.tagAlign.gz

# ENCSR000DMT SE full IDRv2
tf_workflow.py --debug --name ENCSR000DMT-fullIDR --outf /ENCSR000DMT-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR000DMT/rep1/ENCFF000SBI.fastq.gz \
--rep2 /ENCSR000DMT/rep2/ENCFF000SBK.fastq.gz \
--ctl1 /ENCSR000DMT/ctl/ENCFF000SAZ.fastq.gz \
--ctl2 /ENCSR000DMT/ctl/ENCFF000SAZ.fastq.gz

# ENCSR000DMT SE full IDRv2 nomap
tf_workflow.py --debug --name ENCSR000DMT-fullIDR --outf /ENCSR000DMT-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --nomap --yes \
--idrversion 2 \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000DMT/rep1/ENCFF000SBI.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /ENCSR000DMT/rep2/ENCFF000SBK.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /ENCSR000DMT/ctl/ENCFF000SAZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /ENCSR000DMT/ctl/ENCFF000SAZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz

# ENCSR769ZTN PE full IDRv2
tf_workflow.py --debug --name ENCSR769ZTN-fullIDR --outf /ENCSR769ZTN-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /ENCSR769ZTN/rep1/ENCFF002ELM.fastq.gz /ENCSR769ZTN/rep1/ENCFF002ELL.fastq.gz \
--rep2 /ENCSR769ZTN/rep2/ENCFF002ELK.fastq.gz /ENCSR769ZTN/rep2/ENCFF002ELJ.fastq.gz \
--ctl1 /ENCSR769ZTN/ctl1/ENCFF002EFT.fastq.gz /ENCSR769ZTN/ctl1/ENCFF002EFQ.fastq.gz \
--ctl2 /ENCSR769ZTN/ctl2/ENCFF002EFU.fastq.gz /ENCSR769ZTN/ctl2/ENCFF002EFS.fastq.gz

# ENCSR769ZTN PE full IDRv2 nomap
tf_workflow.py --debug --name ENCSR769ZTN-fullIDR --outf /ENCSR769ZTN-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --nomap --yes \
--rep1pe true --rep2pe true \
--idrversion 2 \
--rep1 /ENCSR769ZTN/rep1/ENCFF002ELMENCFF002ELL.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--rep2 /ENCSR769ZTN/rep2/ENCFF002ELKENCFF002ELJ.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--ctl1 /ENCSR769ZTN/ctl1/ENCFF002EFTENCFF002EFQ.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz \
--ctl2 /ENCSR769ZTN/ctl2/ENCFF002EFUENCFF002EFS.raw.srt.filt.srt.nodup.PE2SE.tagAlign.gz

# ENCSR795HTY PE full IDRv2
tf_workflow.py --debug --name ENCSR795HTY-fullIDR --outf /ENCSR795HTY-fullIDR --idr --yes \
--idrversion 2 \
--rep1 /ENCSR795HTY/rep1/ENCFF240UHP.fastq.gz /ENCSR795HTY/rep1/ENCFF555NNB.fastq.gz \
--rep2 /ENCSR795HTY/rep2/ENCFF859GVE.fastq.gz /ENCSR795HTY/rep2/ENCFF974VFA.fastq.gz \
--ctl1 /ENCSR795HTY/ctl1/ENCFF445UEI.fastq.gz /ENCSR795HTY/ctl1/ENCFF141YIN.fastq.gz \
--ctl2 /ENCSR795HTY/ctl2/ENCFF210VLR.fastq.gz /ENCSR795HTY/ctl2/ENCFF057GZE.fastq.gz

# ENCSR795HTY PE full IDRv2 nomap
tf_workflow.py --debug --name ENCSR795HTY-fullIDR --outf /ENCSR795HTY-fullIDRv2-$(date +"%Y%m%d%H%M") --idr --nomap --yes \
--idrversion 2 \
--rep1pe true --rep2pe true \
--rep1 file-BZ05pb80zp5PFPKxXxXV5ZV8 \
--rep2 file-BZ0KG380YBv2BV14PxXqYfBJ \
--ctl1 file-BZ0pPy00ByzPFPKxXxXV8BJ4 \
--ctl2 file-BZ1Bp0j0xByP4350J52qPb7Z


## chr1 extracts

# ENCSR000EEB chr1 IDR only from TA
tf_workflow.py --debug --name ENCSR000EEB-fullIDRtachr1 --outf /ENCSR000EEB-fullIDRtachr1 --idronly --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr1.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr1.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr1.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr1.tagAlign.gz


## chr21 extracts

# ECSR000EEB chr21 SE full IDRv2
tf_workflow.py --debug --name ENCSR000EEBchr21-fullIDR --outf /ENCSR000EEBchr21-fullIDR-$(date +"%Y%m%d%H%M") --idr --yes \
--idrversion 2 \
--rep1 /test_data/ENCFF000XUL.chr21.fq.gz \
--rep2 /test_data/ENCFF000XUK.chr21.fq.gz \
--ctl1 /test_data/ENCFF000XTF.chr21.fq.gz \
--ctl2 /test_data/ENCFF000XTF.chr21.fq.gz

# ENCSR000EEB chr21 IDR only from TA
tf_workflow.py --debug --name ENCSR000EEB-fullIDRtachr21 --outf /ENCSR000EEB-fullIDRtachr21-$(date +"%Y%m%d%H%M") --idronly --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr21.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz

# ENCSR000EEB chr21 IDRv2 only from TA
tf_workflow.py --debug --name ENCSR000EEB-fullIDRtachr21 --outf /ENCSR000EEB-fullIDRtachr21-$(date +"%Y%m%d%H%M") --idronly --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr21.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz

## histones

# ENCSR678FIT chr19 IDRv2 from TA
histone_workflow.py --debug --name ENCSR678FIT-chr19-ta-IDR2 --outf /ENCSR678FIT-chr19-ta-IDR2-$(date +"%Y%m%d%H%M") --idr --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--rep2 /test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl1 /test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl2 /test_data/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT IDRv2 from TA
histone_workflow.py --debug --name ENCSR678FIT-ta-IDR2 --outf /ENCSR678FIT-ta-IDR2-$(date +"%Y%m%d%H%M") --idr --idrversion 2 --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep1/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep2/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT chr19 overlap only from TA
histone_workflow.py --debug --name ENCSR678FIT-chr19-ta-OL --outf /ENCSR678FIT-chr19-ta-OL-$(date +"%Y%m%d%H%M") --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /test_data/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--rep2 /test_data/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl1 /test_data/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--ctl2 /test_data/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.chr19.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR678FIT overlap only from TA
histone_workflow.py --debug --name ENCSR678FIT-ta-OL --outf /ENCSR678FIT-ta-OL-$(date +"%Y%m%d%H%M") --nomap --yes \
--rep1pe false --rep2pe false \
--rep1 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep1/ENCFF926URZ.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--rep2 /mm10_mapping/e165_experiments/bams/ENCSR678FIT/rep2/ENCFF593LFI-ENCFF919IQP_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl1 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep1/ENCFF002BXW.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--ctl2 /mm10_mapping/e165_controls/bams/ENCSR817FFF/rep2/ENCFF482JNU-ENCFF249RWY-ENCFF118FEW_pooled.raw.srt.filt.nodup.srt.SE.tagAlign.gz \
--genomesize mm --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes"

# ENCSR


## single applets

# Run IDRv1 applet with ENCSR000EEB chr21 rep1 vs rep2
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL-chr21.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK-chr21.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL-chr21-ENCFF000XUK-chr21_pooled.regionPeak.gz" \
--input "idr_version=1" \
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
--input "idr_version=2" \
--verbose \
--destination /IDRv2_test \
--name IDRv2_test \
--delay-workspace-destruction \
--priority high \
--yes \
idr

# Run IDRv2 applet with ENCSR000EEB rep1 vs rep2
dx run \
--input "rep1_peaks=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE.fixcoord.regionPeak.gz" \
--input "rep2_peaks=/test_data/ENCFF000XUK.raw.srt.filt.nodup.srt.SE.regionPeak.gz" \
--input "pooled_peaks=/test_data/ENCFF000XUL.raw.srt.filt.nodup.srt.SE-ENCFF000XUK.raw.srt.filt.nodup.srt.SE_pooled.regionPeak.gz" \
--input "idr_version=2" \
--verbose \
--destination /IDRv2_test \
--name IDRv2_test \
--delay-workspace-destruction \
--priority high \
--yes \
idr

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
