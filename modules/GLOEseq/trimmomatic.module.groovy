Trimmomatic = {
	doc title: "Trimmomatic",
		desc:  "A flexible read trimming tool for Illumina NGS data",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"

	output.dir   = TRIMREADS
	def Trimmomatic_FLAGS = "SE -phred33 " +
                            Trimmomatic_THREADS + " "
                            
        def Trimmomatic_OPTIONS = "ILLUMINACLIP:/fsimb/common/tools/trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"
	
	transform(".fastq.gz") to ("_trimmed.fastq.gz") {
		exec """
			module load trimmomatic/${TRIMMOMATIC_VERSION} &&

                        java -jar ${TOOL_TRIMMOMATIC}/trimmomatic-0.36.jar $Trimmomatic_FLAGS $input $output $Trimmomatic_OPTIONS	
		""","Trimmomatic"
	}

}
