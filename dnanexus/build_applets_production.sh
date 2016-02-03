#!/usr/bin/env bash

DEFAULT_PROJECT='ENCODE - ChIP Production'
DEFAULT_FOLDER='applets/'

project=${1:-$DEFAULT_PROJECT}
folder=${2:-$DEFAULT_FOLDER}

APPLETS=('input_shield' 'accession_analysis' 'encode_bwa' 'filter_qc' 'xcor' 'xcor_only' 'spp' 'pool' 'pseudoreplicator' 'encode_spp' 'encode_macs2' 'macs2' 'idr' 'idr2' 'encode_idr' 'overlap_peaks')

for appl in ${APPLETS[@]}; do
		dest="$project:$folder$appl"
	echo $dest
	dx build --archive --destination "$dest" "$appl/"
done

