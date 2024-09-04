ESSENTIAL_PROJECT="..."   // full path to the project directory
ESSENTIAL_TRIMMOMATIC="/data/genomes/TruSeq3-SE.fa" // full path to the adapter fasta file
ESSENTIAL_MIN_READLENGTH=36
ESSENTIAL_BOWTIE_REF="/data/genomes/hg38/hg38" // full path to the reference index files for Bowtie2
ESSENTIAL_BOWTIE_GENOME_INDEX="/data/genomes/hg38/hg38.chrom.sizes" // chromosome sizes file of the reference genome
ESSENTIAL_PAIRED="no" // single-end mode ("no") or paired-end mode ("yes") according to the sequencing mode.
ESSENTIAL_TXDB="TxDb.Hsapiens.UCSC.hg38.knownGene" // needed for breaks annotation
ESSENTIAL_ANNODB="org.Hs.eg.db" // needed for breaks annotation
ESSENTIAL_TSS=200 // TSS region parameter from -200 to +200
ESSENTIAL_FRAGLEN=200 // mean length of library inserts and also minimum peak size called by MACS2
ESSENTIAL_MACS2_GSIZE="2.94E+09" // mappable genome size for MACS2 
ESSENTIAL_THREADS=4 // number of threads for parallel tasks
ESSENTIAL_DUP="all" // relevant for MACS2, it instructs it to keep duplicate reads
ESSENTIAL_BEDCOVERAGE="--smoothLength 1 --binSize 1 --normalizeUsing BPM"  // deepTools options for making normalised bigWig tracks
                                                                           // If you want to exclude chromsomes for normalisation e.g. chrM add 
                                                                           // the following parameter --ignoreForNormalization chrM
ESSENTIAL_RFD_BINSIZE=1000 // Size of the bins, in bases, for the output of the RFD bigwig file

ESSENTIAL_EB1=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/BsrDI.bed" // Expected BsrDI breaks in the yeast genome
ESSENTIAL_EB2=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/Nb_BsrDI.bed" // Expected NbBsrDI breaks in the yeast genome
ESSENTIAL_EB3=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/NotI.bed" // Expected NotI breaks in the yeast genome

ESSENTIAL_RE_3end=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/NotI_3end_saccer3.bed" // Expected 3'end breaks in the genome
ESSENTIAL_RE_5end=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/NotI_5end_saccer3.bed" // Expected 5'end breaks in the genome
ESSENTIAL_QPCR=ESSENTIAL_PROJECT + "/GLOEPipe/tools/REs/qpcr_perc_cut_efficiency.txt" // qPCR cutting efficiency results

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MERGEDREADS=PROJECT + "/rawdata_merged"
UMIREADS=PROJECT + "/rawdata_umi"
TRIMREADS=PROJECT + "/rawdata_trimmed"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
TRACKS_NORMALIZED=PROJECT + "/tracks_normalized"
TRACKS_RE_NORMALIZED=PROJECT + "/tracks_re_normalized"
BG=PROJECT + "/tracks/bg"
SS=PROJECT + "/tracks/strandspecific"
RFD=PROJECT + "/tracks/strandspecific/rfd"
SS_NORMALIZED=PROJECT + "/tracks_normalized/strandspecific"
ALPHA=RESULTS + "/alpha"
SS_RE_NORMALIZED=PROJECT + "/tracks_re_normalized/strandspecific"
BED=PROJECT + "/bed"
READS=PROJECT + "/bed/reads"
SITES=PROJECT + "/bed/sites"
COUNTS=PROJECT + "/bed/counts"
TPM=PROJECT + "/bed/tpm"

//Necessary for Fastqscreen
ESSENTIAL_BOWTIE_Human="/data/FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38"
ESSENTIAL_BOWTIE_Yeast="/data/FastQ_Screen_Genomes/Yeast/Saccharomyces_cerevisiae.R64-1-1"
ESSENTIAL_BOWTIE_Bovine="/data/FastQ_Screen_Genomes/Bovine/ARS-UCD2.0"
ESSENTIAL_BOWTIE_Bacteria="/data/FastQ_Screen_Genomes/E_coli/Ecoli"
ESSENTIAL_BOWTIE_Mycoplasma="/data/FastQ_Screen_Genomes/Mycoplasma/Mycoplasma"
ESSENTIAL_BOWTIE_PhiX="/data/FastQ_Screen_Genomes/PhiX/phi_plus_SNPs"
ESSENTIAL_BOWTIE_ERCC="/data/FastQ_Screen_Genomes/ERCC/ERCC92"
ESSENTIAL_BOWTIE_rRNA="/data/FastQ_Screen_Genomes/rRNA/GRCm38_rRNA"
