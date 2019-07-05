//vars for task macs2 from catalog ChIPseq, version 1
MACS2_TARGETS="targets.txt" // targets file describing the samples
MACS2_MFOLD="-m 5 50"	// range of enrichment ratio (default: 5,50)
MACS2_GSIZE="-g " + ESSENTIAL_MACS2_GSIZE // the mappable genome size
MACS2_BWIDTH="--bw " + Integer.toString(ESSENTIAL_FRAGLEN)	  // bandwidth use for model building
MACS2_INPUT=BED // where the bed files are stored
MACS2_FORMAT="--format BED"
MACS2_EXTRA="--extsize 1 --nomodel --shift 0 --keep-dup " + ESSENTIAL_DUP		// other parms sent to macs2

