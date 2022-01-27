#!/bin/bash

# Check for the number of args
#if [ $# -ne 1 ]; then
#        echo "Please provide the correct number of arguments..."
#        echo "USAGE: plot_heatmap.sh MATRIX_FILE"
#        exit 1
#fi

# Get the matrix file as input
matrix=$1
color_scheme=${2:-'Reds'}
#values=$3

# Get the name for the output heatmap file
outfile=${matrix%.gz}.pdf

if [ $color_scheme == 'Reds' ]; then
	# Plot a heatmap
	plotHeatmap -m $matrix --dpi 300 --colorMap 'Reds' \
		--missingDataColor '#FFFFFF' --averageTypeSummaryPlot mean \
		--whatToShow 'heatmap and colorbar' -x "Distance" -z "Peaks" -o $outfile
else
	# Specify the colors for the heatmaps
	colors=$(perl -pe 'chomp if eof' $color_scheme | sed 's/^/\"white,/' | sed 's/$/\"/' | tr '\n' ' ')
	echo $colors

	# Specify the max values for the heatmaps
	max_values=$(perl -pe 'chomp if eof' $values | tr '\n' ' ')
   
	missing_color='#FFFFFF'
	echo $max_values

	# Plot a heatmap
	plotHeatmap -m ${matrix} --dpi 300 --colorList ${colors}\
	    --missingDataColor '#FFFFFF' --averageTypeSummaryPlot mean\
	    --zMax ${max_values} --whatToShow 'heatmap and colorbar'\
	    -x "Distance" -z "Peaks" -o $outfile
fi
