#!/bin/bash

# Use the key and download data in the current directory
rsync -av --progress -e "ssh -i /lab/solexa_young/lazaris/keys/AleData.pem" AleData@3.91.148.26:~/20210819* .
