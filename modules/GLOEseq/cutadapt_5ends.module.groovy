Cutadapt_5ends = {
	doc title: "Cutadapt module for paired-end data",
		desc:  "A flexible read trimming tool for Illumina NGS data",
		bpipe_version: "tested with bpipe 0.9.9.8.slurm",
		author: "Giuseppe Petrosino, Anke Busch"

	output.dir   = TRIMREADS

	// create the log folder if it doesn't exists
    def TRIM_LOGS = new File(LOGS + "/trimAdapter")
    if (!TRIM_LOGS.exists()) {
            TRIM_LOGS.mkdirs()
    }


    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_PRUNED = SAMPLENAME_BASE.replace(".R1.m.end5", "") // delete .R1 in combined log file of pe design


	// Cutadapt options
	def CUTADAPT_FLAGS = "-a file:" + CUTADAPT_ADAPTER + 
                         " -A file:" + CUTADAPT_ADAPTER + 
                             " --overlap=" + CUTADAPT_OVERLAP + 
                             " --minimum-length=" + CUTADAPT_MINREADLEN +
                             " --error-rate=" + CUTADAPT_ERROR_RATE + 
                             " --nextseq-trim=" + CUTADAPT_NEXTSEQTRIM +
                             " " + CUTADAPT_EXTRA

	transform("*.fastq.gz") to (".trimmed.fastq.gz") {
		exec """

            module load cutadapt/${CUTADAPT_VERSION} &&

                        if [ -n "\$SLURM_JOBID" ]; then
                                export TMPDIR=/jobdir/\${SLURM_JOBID};
                        fi &&

                        cutadapt $CUTADAPT_THREADS $CUTADAPT_FLAGS --output=\${TMPDIR}/${SAMPLENAME_BASE_PRUNED}.R1.m.end5.trimmed.fastq.gz --paired-output \${TMPDIR}/${SAMPLENAME_BASE_PRUNED}.R2.m.end5.trimmed.fastq.gz $input1 $input2 1> ${TRIM_LOGS}/${SAMPLENAME_BASE_PRUNED}.end5.cutadapt.log &&
                        mv -t $output.dir \${TMPDIR}/${SAMPLENAME_BASE_PRUNED}*.trimmed.fastq.gz 

		""","Cutadapt_5ends"
	}

}