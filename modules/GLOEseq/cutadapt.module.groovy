Cutadapt = {
	doc title: "Cutadapt module for single-end data",
		desc:  "A flexible read trimming tool for Illumina NGS data",
		bpipe_version: "tested with bpipe 0.9.9.8.slurm",
		author: "Giuseppe Petrosino, Anke Busch"

	output.dir   = TRIMREADS

	// create the log folder if it doesn't exists
        def TRIM_LOGS = new File(LOGS + "/trimAdapter")
        if (!TRIM_LOGS.exists()) {
                TRIM_LOGS.mkdirs()
        }

	// extract the base of the input file name (w/o the directory, w/o fastq.gz)
        def INPUTBASE = input
        int path_index = INPUTBASE.lastIndexOf("/")
        INPUTBASE = INPUTBASE.substring(path_index+1)
        INPUTBASE = (INPUTBASE =~ /.fastq.gz/).replaceFirst("") 


	// Cutadapt options
	def CUTADAPT_FLAGS = "-a file:" + CUTADAPT_ADAPTER + 
                             " --overlap=" + CUTADAPT_OVERLAP + 
                             " --minimum-length=" + CUTADAPT_MINREADLEN +
                             " --error-rate=" + CUTADAPT_ERROR_RATE + 
                             " --nextseq-trim=" + CUTADAPT_NEXTSEQTRIM +
                             " " + CUTADAPT_EXTRA

	
	transform(".fastq.gz") to ("_trimmed.fastq.gz") {
		exec """

                        module load cutadapt/${CUTADAPT_VERSION} &&

                        if [ -n "\$SLURM_JOBID" ]; then
                                export TMPDIR=/jobdir/\${SLURM_JOBID};
                        fi &&

                        cutadapt $CUTADAPT_FLAGS --output=$output $input 1> ${TRIM_LOGS}/${INPUTBASE}.cutadapt.log 

		""","Cutadapt"
	}

}
