ESSENTIAL_PROJECT=".../project"   // full path to the project directory
ESSENTIAL_TRIMMOMATIC=".../TruSeq3-SE.fa" // full path to the adapter fasta file
ESSENTIAL_MIN_READLENGTH=36
ESSENTIAL_BOWTIE_REF=".../genome" // full path to the reference index files for Bowtie2
ESSENTIAL_BOWTIE_GENOME_INDEX="../saccer3.chrom.sizes" // chromosome sizes file of the reference genome
ESSENTIAL_PAIRED="no" // single-end mode ("no") or paired-end mode ("yes") according to the sequencing mode.
ESSENTIAL_TXDB="TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" // needed for breaks annotation
ESSENTIAL_ANNODB="org.Sc.sgd.db" // needed for breaks annotation
ESSENTIAL_TSS=200 // TSS region parameter from -200 to +200
ESSENTIAL_FRAGLEN=200 // mean length of library inserts and also minimum peak size called by MACS2
ESSENTIAL_MACS2_GSIZE="1.20E+07" // mappable genome size for MACS2 
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
ESSENTIAL_BOWTIE_Human="/fsimb/common/genomes/homo_sapiens/gencode/release-25_GRCh38.p7/full/index/bowtie2/GRCh38.p7.genome"
ESSENTIAL_BOWTIE_Yeast="/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/bowtie2/2.3.2/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel"
ESSENTIAL_BOWTIE_Bovine="/fsimb/common/genomes/bos_taurus/ensembl/3.1/full/index/bowtie2/2.2.9/UMD3.1"
ESSENTIAL_BOWTIE_Bacteria="/fsimb/common/genomes/Escherichia_coli/ensembl/full/index/bowtie2/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.31.dna.genome"
ESSENTIAL_BOWTIE_Mycoplasma="/fsimb/common/genomes/contaminants/fastqscreen_references/mycoplasma/v1/index/bowtie2/2.3.4.3/mycoplasma_all_ref"
ESSENTIAL_BOWTIE_PhiX="/fsimb/common/genomes/phix/19930428/NCBI/index/bowtie2/2.3.4.3/ncbi_phix"
ESSENTIAL_BOWTIE_ERCC="/fsimb/common/genomes/ERCC/index/bowtie2/2.3.4.3/ERCC92"
ESSENTIAL_BOWTIE_rRNA="/fsimb/common/genomes/contaminants/fastqscreen_references/rrna/v1/index/bowtie2/2.3.4.3/hs_mm_ce_dm_rn_dr_xt_rRNA"