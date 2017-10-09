#!/usr/bin/env bash

DEFAULT_FOLDER='/ChIP-seq/'

project=$1
folder=${2:-$DEFAULT_FOLDER}

ASSETS=('awscli_asset' 'bedtools_asset' 'bioconductor_asset' 'common_asset' 'macs2_asset' 'spp_asset')
dx mkdir -p "$project:$folder/assets/"

for asset in ${ASSETS[@]}; do
	dest="$project:$folder/assets/$asset"
	echo $dest
	dx build_asset --destination "$project:$folder/assets/$asset" "$asset/"
done

