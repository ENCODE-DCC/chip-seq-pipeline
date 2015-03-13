#!/usr/bin/env bash

# full IDR template
tf_workflow.py --debug --name IDR_Template --outp "E3 ChIP-seq"

# ECSR000EEB SE full IDR
tf_workflow.py --debug --name ENCSR000EEB-fullIDR --outf /ENCSR000EEB-fullIDR --idr --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ENCSR000EEB SE full IDR only (no naive peaks)
tf_workflow.py --debug --name ENCSR000EEB-fullIDRonly --outf /ENCSR000EEB-fullIDRonly --idronly --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL.fastq.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK.fastq.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF.fastq.gz

# ENCSR000EEB chr1 IDR only from TA
tf_workflow.py --debug --name ENCSR000EEB-fullIDRtachr1 --outf /ENCSR000EEB-fullIDRtachr1 --idronly --nomap --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr1.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr1.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr1.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr1.tagAlign.gz

# ENCSR000EEB chr21 IDR only from TA
tf_workflow.py --debug --name ENCSR000EEB-fullIDRtachr21  --outf /ENCSR000EEB-fullIDRtachr21 --idronly --nomap --yes \
--rep1 /ENCSR000EEB/rep1/ENCFF000XUL-chr21.tagAlign.gz \
--rep2 /ENCSR000EEB/rep2/ENCFF000XUK-chr21.tagAlign.gz \
--ctl1 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz \
--ctl2 /ENCSR000EEB/ctl/ENCFF000XTF-chr21.tagAlign.gz

# ENCSR000BUA SE full IDR
tf_workflow.py --debug --name ENCSR000BUA-fullIDR --outf /ENCSR000BUA-fullIDR --idr --yes \
--rep1 /ENCSR000BUA/rep1/ENCFF000RBI.fastq.gz \
--rep2 /ENCSR000BUA/rep2/ENCFF000RBC.fastq.gz \
--ctl1 /ENCSR000BUA/ctl1/ENCFF000RCK.fastq.gz \
--ctl2 /ENCSR000BUA/ctl2/ENCFF000RCP.fastq.gz

# ENCSR000DMT SE full IDR
tf_workflow.py --debug --name ENCSR000DMT-fullIDR --outf /ENCSR000DMT-fullIDR --idr --yes \
--rep1 /ENCSR000DMT/rep1/ENCFF000SBI.fastq.gz \
--rep2 /ENCSR000DMT/rep2/ENCFF000SBK.fastq.gz \
--ctl1 /ENCSR000DMT/ctl/ENCFF000SAZ.fastq.gz \
--ctl2 /ENCSR000DMT/ctl/ENCFF000SAZ.fastq.gz

# ENCSR769ZTN PE full IDR
tf_workflow.py --debug --name ENCSR769ZTN-fullIDR --outf /ENCSR769ZTN-fullIDR --idr --yes \
--rep1 /ENCSR769ZTN/rep1/ENCFF002ELM.fastq.gz /ENCSR769ZTN/rep1/ENCFF002ELL.fastq.gz \
--rep2 /ENCSR769ZTN/rep2/ENCFF002ELK.fastq.gz /ENCSR769ZTN/rep2/ENCFF002ELJ.fastq.gz \
--ctl1 /ENCSR769ZTN/ctl1/ENCFF002EFT.fastq.gz /ENCSR769ZTN/ctl1/ENCFF002EFQ.fastq.gz \
--ctl2 /ENCSR769ZTN/ctl2/ENCFF002EFU.fastq.gz /ENCSR769ZTN/ctl2/ENCFF002EFS.fastq.gz

# ENCSR795HTY PE full IDR
tf_workflow.py --debug --name ENCSR795HTY-fullIDR --outf /ENCSR795HTY-fullIDR --idr --yes \
--rep1 /ENCSR795HTY/rep1/ENCFF240UHP.fastq.gz /ENCSR795HTY/rep1/ENCFF555NNB.fastq.gz \
--rep2 /ENCSR795HTY/rep2/ENCFF859GVE.fastq.gz /ENCSR795HTY/rep2/ENCFF974VFA.fastq.gz \
--ctl1 /ENCSR795HTY/ctl1/ENCFF445UEI.fastq.gz /ENCSR795HTY/ctl1/ENCFF141YIN.fastq.gz \
--ctl2 /ENCSR795HTY/ctl2/ENCFF210VLR.fastq.gz /ENCSR795HTY/ctl2/ENCFF057GZE.fastq.gz
