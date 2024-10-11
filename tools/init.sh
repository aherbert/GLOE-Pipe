#!/bin/sh

# Initialise the directory containing the GLOE-Pipe source for use
# as a project directory:
# <project dir>/GLOEPipe

if [ ! -d GLOEPipe ]
then
    echo "Missing GLOEPipe source directory"
    exit 1
fi
ln -s GLOEPipe/pipelines/GLOEseq/bpipe.config . 
ln -s GLOEPipe/pipelines/GLOEseq/gloeseq.pipeline.groovy . 
ln -s GLOEPipe/modules/GLOEseq/essential.vars.groovy .

# Update the project directory
DIR=`pwd`
sed -i "s@ESSENTIAL_PROJECT=.*@ESSENTIAL_PROJECT=\"${DIR}\"@" essential.vars.groovy

if [ ! -f targets.txt ]
then
	head -n 1 GLOEPipe/pipelines/GLOEseq/targets.txt > targets.txt
fi

echo "Configure targets.txt for the input fastq data files"
