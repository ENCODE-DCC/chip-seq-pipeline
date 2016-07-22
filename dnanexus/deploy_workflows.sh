chip_workflow.py \
	--name  "ENCODE TF ChIP-seq (hg19)" \
	--title "TF ChIP-seq" \
	--target tf \
	--genomesize hs \
	--reference     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.tar.gz" \
	--chrom_sizes   "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
	--narrowpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/narrowPeak.as" \
	--gappedpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/gappedPeak.as" \
	--broadpeak_as  "ENCODE Uniform Processing Pipelines:/Reference Files/as files/broadPeak.as" \
	--outp          "ENCODE Uniform Processing Pipelines" \
	--outf			"/ChIP-seq/test_workflows/" \
	--applets       "ENCODE Uniform Processing Pipelines" \
	--blacklist     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

chip_workflow.py \
	--name  "ENCODE TF ChIP-seq Unary Control (hg19)" \
	--title "TF ChIP-seq" \
	--target tf \
	--genomesize hs \
	--unary_control \
	--reference     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.tar.gz" \
	--chrom_sizes   "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
	--narrowpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/narrowPeak.as" \
	--gappedpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/gappedPeak.as" \
	--broadpeak_as  "ENCODE Uniform Processing Pipelines:/Reference Files/as files/broadPeak.as" \
	--outp          "ENCODE Uniform Processing Pipelines" \
	--outf			"/ChIP-seq/test_workflows/" \
	--applets       "ENCODE Uniform Processing Pipelines" \
	--blacklist     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# chip_workflow.py \
# 	--name  "Peaks-only TF ChIP-seq (hg19)" \
# 	--title "TF ChIP-seq Peaks" \
# 	--nomap \
# 	--target tf \
# 	--genomesize hs \
# 	--reference     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.tar.gz" \
# 	--chrom_sizes   "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
# 	--narrowpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/narrowPeak.as" \
# 	--gappedpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/gappedPeak.as" \
# 	--broadpeak_as  "ENCODE Uniform Processing Pipelines:/Reference Files/as files/broadPeak.as" \
# 	--outp          "ENCODE Uniform Processing Pipelines" \
# 	--outf			"/ChIP-seq/test_workflows/" \
# 	--applets       "ENCODE Uniform Processing Pipelines" \
# 	--blacklist     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

# chip_workflow.py \
# 	--name  "Peaks-only TF ChIP-seq Unary Control (hg19)" \
# 	--title "TF ChIP-seq Peaks" \
# 	--nomap \
# 	--target tf \
# 	--genomesize hs \
# 	--unary_control \
# 	--reference     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.tar.gz" \
# 	--chrom_sizes   "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/male.hg19.chrom.sizes" \
# 	--narrowpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/narrowPeak.as" \
# 	--gappedpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/gappedPeak.as" \
# 	--broadpeak_as  "ENCODE Uniform Processing Pipelines:/Reference Files/as files/broadPeak.as" \
# 	--outp          "ENCODE Uniform Processing Pipelines" \
# 	--outf			"/ChIP-seq/test_workflows/" \
# 	--applets       "ENCODE Uniform Processing Pipelines" \
# 	--blacklist     "ENCODE Uniform Processing Pipelines:/Reference Files/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

chip_workflow.py \
	--name  "ENCODE TF ChIP-seq (no reference)" \
	--title "TF ChIP-seq" \
	--target tf \
	--narrowpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/narrowPeak.as" \
	--gappedpeak_as "ENCODE Uniform Processing Pipelines:/Reference Files/as files/gappedPeak.as" \
	--broadpeak_as  "ENCODE Uniform Processing Pipelines:/Reference Files/as files/broadPeak.as" \
	--outp          "ENCODE Uniform Processing Pipelines" \
	--outf			"/ChIP-seq/test_workflows/" \
	--applets       "ENCODE Uniform Processing Pipelines"
