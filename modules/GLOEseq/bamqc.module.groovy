BamQC = {
	doc title: "BamQC",
		desc:  "Quality control of bam file",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"
    
	output.dir   = BamQC_OUTDIR
	def BamQC_FLAGS = "--extract --quiet"
	
	transform(".bam") to ("_bamqc.zip") {
		exec """
			module load BamQC/${BamQC_VERSION} &&

			bamqc $BamQC_FLAGS -o $output.dir $input
		""","BamQC"
	}

	forward input
}
