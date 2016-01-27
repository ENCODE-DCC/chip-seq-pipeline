#!/usr/bin/env bash

DEFAULT_PROJECT='ENCODE - ChIP Production'

project=${1:-$DEFAULT_PROJECT}
APPLETS=('input_shield' 'accession_analysis' 'encode_bwa' 'filter_qc' 'xcor' 'xcor_only' 'spp' 'pool' 'pseudoreplicator' 'encode_spp' 'encode_macs2' 'macs2' 'idr' 'idr2' 'encode_idr' 'overlap_peaks')

for appl in ${APPLETS[@]}; do
	dest="$project:/applets/$appl"
	echo $dest
	dx build --archive --destination "$dest" "$appl/"
done

