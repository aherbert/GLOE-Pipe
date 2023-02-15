FastQC = {
	doc title: "FastQC",
		desc:  "Quality control of input file",
		constraints: "Only supports compressed FASTQ files",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino, Anke Busch"
        
        var subdir : ""
        output.dir = FASTQC_OUTDIR + "/$subdir"     
	def FASTQC_FLAGS = "--extract --quiet" + FASTQC_THREADS
	
	transform("*.fastq.gz") to ("_fastqc.zip") {
		exec """
			module load fastqc/${FASTQC_VERSION} &&

			fastqc $FASTQC_FLAGS -o $output.dir $inputs
		""","FastQC"
	}

	forward input
}
