bed2bz2 = {
	doc title: "bed2bz2",
		desc:  "Compress Bed file",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"

	output.dir=BED
	def PBZIP2_FLAG = " "  + PBZIP2_THREADS + " "
	
		exec """
               pbzip2 $PBZIP2_FLAG ${BED}/*.bed

		""","bed2bz2"
}