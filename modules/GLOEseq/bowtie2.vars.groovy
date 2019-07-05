// bowtie parameters with suggested typical defaults
BOWTIE_THREADS=" -p" + Integer.toString(ESSENTIAL_THREADS) // threads to use
BOWTIE_SAMTOOLS_THREADS="-@" + Integer.toString(ESSENTIAL_THREADS)
BOWTIE_REF=ESSENTIAL_BOWTIE_REF // prefix of the bowtie reference genome
