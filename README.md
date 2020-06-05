![IMB-logo](resources/IMB_logo.png)

# GLOE-Pipe

The **G**enome-wide **L**igation of 3'-**O**H **E**nds sequencing (GLOE-Seq) is a sensitive method to capture DNA strand breaks genomewide at nucleotide resolution.

To facilitate the analysis of the raw GLOE-Seq data, we have developed an easy-to-use, modular and versatile computational pipeline, called **GLOE-Pipe**.

Given that the sequences produced by GLOE-Seq are the reverse-complement of the original captured DNA fragment, GLOE-Pipe assigns the 5'-end of a read to the strand opposite to that on which it mapped. 

Two output modes are available:
- in the **direct** mode, each read represents one unit of signal that is positioned exactly at the 3'-terminal nucleotide of the original captured fragment. This mode would normally be used to visualise strand breaks
- the **indirect** mode is essentially identical, but it positions the signal one nucleotide immediately upstream of the 5'-end of the original captured fragment, which corresponds to the nucleotide immediately 3' of the break. The indirect mode would therefore normally be used to visualise the position of modified or damaged bases (under the premise that the enzymatic treatment used to detect this modification or lesion generates a nick 5’ to the affected nucleotide)

A flowchart of GLOE-Pipe, detailing each one of its steps, can be found [here](https://app.diagrams.net/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=GLOEseq_pipeline.html#Uhttps%3A%2F%2Fdrive.google.com%2Fuc%3Fid%3D12Ke7Tz_CBC-Hke5WHT6FVtbFMwhMCQew%26export%3Ddownload).

## System Requirements
GLOE-Pipe is designed to run in a high-performance cluster (HPC) environment with a Linux distribution (*e.g.* Debian 9), Bpipe (version 0.9.9.3) as the domain specific language (DSL), and Lmod (version 6.6) as the module system.

## Prerequisites

### Software

| Software | Version used | Link |
| --- | --- | --- |
| `FastQC` | 0.11.5 | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| `Trimmomatic` | 0.36 | http://www.usadellab.org/cms/?page=trimmomatic |
| `Bowtie2` | 2.3.4 | http://bowtie-bio.sourceforge.net/bowtie2/index.shtml |
| `Samtools` | 1.5 | http://samtools.sourceforge.net/ |
| `BEDTools` | 2.25.0 | https://bedtools.readthedocs.io/en/latest/ |
| `bedGraphToBigWig` | 365 | https://github.com/ENCODE-DCC/kentUtils |
| `MACS2` | 2.1.1 | https://github.com/taoliu/MACS |
| `deepTools` | 3.1.0 | https://deeptools.readthedocs.io/en/develop |
| `R` | 3.5.1 | https://www.r-project.org/ |
| `ChIPseeker` | 1.8.0 | https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html |
| `rtracklayer` | 1.42.2 | https://bioconductor.org/packages/release/bioc/html/rtracklayer.html |
| `regioneR` | 1.14.0 | https://bioconductor.org/packages/release/bioc/html/regioneR.html |
| `Cairo` | 1.15-10 | https://cran.r-project.org/web/packages/Cairo/index.html |
| `openxlsx` | 4.1.0 | https://cran.r-project.org/web/packages/openxlsx/index.html |

### Files
- `targets.txt`, which includes sample names
- raw reads (`*.fastq.gz`)

## Before you run GLOE-Pipe

### How to get GLOE-Pipe into your project
This can be done by cloning the most recent development version of GLOE-Pipe from this Git repository

    git clone https://gitlab.com/GPetrosino/GLOE-Pipe.git <project_dir>/GLOEPipe
    cd <project_dir>
    ln -s GLOEPipe/pipelines/GLOEseq/* . 
    ln -s GLOEPipe/modules/GLOEseq/essential.vars.groovy .

### Customise GLOE-Pipe to your needs

Adjust the information found in the following files:

- `gloeseq.pipeline.groovy`, which describes the pipeline steps and the location of the respective modules
- `essential.vars.groovy`, which specifies the main project variables, such as the path to the project directory and the reference genome 
- `targets.txt`, which contains file names and sample comparisons
- `tool.location.groovy` and `bpipe.config`, which specify the paths and resource allocation for each module

## Run GLOE-Pipe

Copy the input FastQ files into `<project_dir>/rawdata folder`.

Using GNU `screen` (for persistence), load the bpipe module customised for the Slurm job manager

    screen
    module load bpipe/0.9.9.3.slurm

Start running GLOE-Pipe (NB: **the indirect mode is the default**)

    bpipe run gloeseq.pipeline.groovy rawdata/*.fastq.gz

## Notes
- To use the direct mode, the main pipeline needs to be changed in the `gloeseq.pipeline.groovy` file
- The `REs` folder contains bed files that include information on the break sites produced by the restriction endonucleases used in the GLOE-seq article
- For the yeast genome, the break sites that overlap with transposons, telomeres and mtDNA are excluded
