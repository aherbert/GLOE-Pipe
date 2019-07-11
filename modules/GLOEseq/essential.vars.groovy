ESSENTIAL_PROJECT=".../project"   // full path to the project directory
ESSENTIAL_BOWTIE_REF=".../genome" // full path to the reference index files for Bowtie2
ESSENTIAL_BOWTIE_GENOME=".../genome.fa" // full path to the reference genome FASTA file
ESSENTIAL_BOWTIE_GENOME_INDEX="../saccer3.chrom.sizes" // chromosome sizes file of the reference genome
ESSENTIAL_SAMPLE_PREFIX=""  // sample prefix to be trimmed in the results and reports
ESSENTIAL_BSGENOME="BSgenome.Scerevisiae.UCSC.sacCer3" // Bioconductor genome reference used by some modules
ESSENTIAL_TXDB="TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" // needed for breaks annotation
ESSENTIAL_ANNODB="org.Sc.sgd.db"                    // needed for breaks annotation
ESSENTIAL_FRAGLEN=200 // mean length of library inserts and also minimum peak size called by MACS2
ESSENTIAL_MACS2_GSIZE="1.20E+07" // mappable genome size for MACS2 
ESSENTIAL_THREADS=8 // number of threads for parallel tasks
ESSENTIAL_DUP="all" // relevant for MACS2, it instructs it to keep duplicate reads


ESSENTIAL_RE1=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/BsrDI.bed" // Expected BsrDI breaks in the yeast genome
ESSENTIAL_RE2=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/Nb_BsrDI.bed" // Expected NbBsrDI breaks in the yeast genome
ESSENTIAL_RE3=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/NotI.bed" // Expected NotI breaks in the yeast genome


//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
TRIMREADS=PROJECT + "/rawdata_trimmed"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
TRACKS_NORMALIZED=PROJECT + "/tracks_normalized"
BG=PROJECT + "/tracks/bg"
SS=PROJECT + "/tracks/strandspecific"
RFD=PROJECT + "/tracks/strandspecific/rfd"
SS_NORMALIZED=PROJECT + "/tracks_normalized/strandspecific"
BED=PROJECT + "/bed"
READS=PROJECT + "/bed/reads"
SITES=PROJECT + "/bed/sites"
COUNTS=PROJECT + "/bed/counts"
TPM=PROJECT + "/bed/tpm"