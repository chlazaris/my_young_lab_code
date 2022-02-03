#!/bin/bash

# Check for the number of args
#if [ $# -ne 1 ]; then
#        echo "Please provide the correct number of arguments..."
#        echo "USAGE: plot_heatmap.sh MATRIX_FILE COLOR_SCHEME [default: Reds] [MAX-VALUES]"
#        exit 1
#fi

# Get the matrix file as input
matrix=$1
color_scheme=${2:-'Reds'}
values=$3

# Get the name for the output heatmap file
outfile=${matrix%.gz}.pdf

if [ $color_scheme == 'Reds' ]; then
	# Plot a heatmap
	echo plotHeatmap -m $matrix --dpi 300 --colorMap 'Reds' \
		--missingDataColor '#FFFFFF' --averageTypeSummaryPlot mean \
		--sortUsingSamples 1 \
		--whatToShow 'heatmap and colorbar' -x "Distance" -z "Peaks" -o $outfile
else
	# Specify the colors for the heatmaps
	colors=$(perl -pe 'chomp if eof' $color_scheme | sed 's/^/\"white,/' | sed 's/$/\"/' | tr '\n' ' ')

	# Specify the max values for the heatmaps
	max_values=$(perl -pe 'chomp if eof' $values | tr '\n' ' ') 

	# Plot a heatmap
	echo plotHeatmap -m $matrix --dpi 300 \
		--missingDataColor "'#FFFFFF'" --averageTypeSummaryPlot mean \
		--sortUsingSamples 1 \
		--whatToShow \"'heatmap and colorbar'\" -x "Distance" -z "Peaks" --zMax $max_values --colorList $colors \
		-o $outfile > command
	chmod a+x command
	./command
	# Get rid of the command file
	rm -rf command
fi
