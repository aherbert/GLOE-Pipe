FastQC = {
	doc title: "FastQC",
		desc:  "Quality control of input file",
		constraints: "Only supports compressed FASTQ files",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino"
    
	output.dir   = FASTQC_OUTDIR
	def FASTQC_FLAGS = "--extract --quiet"
	
	transform(".fastq.gz") to ("_fastqc.zip") {
		exec """
			module load fastqc/${FASTQC_VERSION} &&

			fastqc $FASTQC_FLAGS -o $output.dir $input
		""","FastQC"
	}

	forward input
}
