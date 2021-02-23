#!/bin/bash

# Get the input file
file=$1

# Check if there is one input
if [ $# -ne 1 ]; then
	echo "Please provide the input file with the Google Cloud URLs..."
	echo "USAGE: gcloud_download INPUT_FILE"
	exit 1
fi

# Download the files using gsutil
gsutil -m cp -r -n -I  . < $file
