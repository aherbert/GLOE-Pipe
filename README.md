![IMB-logo](resources/IMB_logo.png)

# GLOE-Pipe #

The Genome-wide Ligation of 3'OH Ends sequencing (GLOE-seq) is a sensitive method to monitor double-strand breaks (DSBs) and single-strand breaks (SSBs) genome-wide at base-pair resolution.

## Prerequisites ##
#### Programs required ####
- Trimmomatic
- FastQC
- Bowtie 2
- Samtools
- BED tools
- UCSC utilities
- MACS2
- ChIPSeeker

#### Files required ####
- targets.txt (sample names)
- raw reads (.fastq.gz) or mapped data (.bam)

## Preparations to run ##

### Get GLOE-Pipe to your project ###
This can be done cloning the most recent development version from the GitLab repository

    git clone https://gitlab.com/GPetrosino/GLOE-Pipe.git <project_dir>/GLOEPipe
    cd <project_dir>
    ln -s GLOEPipe/pipelines/GLOEseq/* . 
    ln -s GLOEPipe/modules/GLOEseq/essential.vars.groovy .  

### Customise GLOE-Pipe to your needs ###

Adjust the information found in the following files:

- *gloeseq.pipeline.groovy* describes the pipeline steps and the location of the respective modules
- *essential.vars.groovy* specifies the main project variables like project dir and reference genome 
- *targets.txt* contains the filenames and sample comparisons
- *tool.location.groovy* and *bpipe.config* specify the paths and resource allocation for the tools

## Run GLOE-Pipe ##

Copy the input FastQ files into the <project_dir>/rawdata folder.

Using GNU Screen (for persistence) load the bpipe module customised for the Slurm job manager, e.g.

    screen
    module load bpipe/0.9.9.3.slurm

Start running GLOE-Pipe (direct mode as default)

    bpipe run gloeseq.pipeline.groovy rawdata/*.fastq.gz

## Notes ##
To use indirect mode the main pipeline needs to be modified in the gloeseq.pipeline.groovy file.

The REs folder contains bed files including break sites of restriction endonucleases used in the GLOE-seq paper.
The break sites overlapping transposons, telomeres and mtDNA on yeast genome are excluded.