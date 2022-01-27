#!/bin/bash
date=$(date '+%Y%m%d')

# Find the top directories in terms of space
du -a /lab/solexa_young/lazaris/ | sort -nr | head -n500 > space_usage/$date_top500_dirs.txt
