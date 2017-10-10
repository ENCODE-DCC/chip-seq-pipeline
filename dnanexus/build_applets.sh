#!/usr/bin/env bash

DEFAULT_FOLDER='/ChIP-seq/'

project=$1
folder=${2:-$DEFAULT_FOLDER}

PRODUCTION_APPLETS=('bam2tagalign' 'encode_map' 'filter_qc' 'xcor' 'xcor_only' 'spp' 'pool' 'pseudoreplicator' 'encode_spp' 'encode_macs2' 'macs2' 'idr2' 'encode_idr' 'overlap_peaks')
ACCESSORY_APPLETS=('input_shield' 'accession_analysis' 'shell')
dx mkdir -p "$project:$folder/applets/"

for applet in ${PRODUCTION_APPLETS[@]}; do
	dest="$project:$folder/applets/$applet"
	echo $dest
	# cp common.py $appl/resources/home/dnanexus/common.py
	dx build --archive --destination "$dest" "$applet/"
done

for applet in ${ACCESSORY_APPLETS[@]}; do
	dest="$project:$folder/applets/$appl"
	echo $dest
	# cp common.py $appl/resources/home/dnanexus/common.py
	dx build --archive --destination "$dest" "$applet/"
done

