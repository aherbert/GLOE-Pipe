Trimmomatic = {
	doc title: "Trimmomatic",
		desc:  "A flexible read trimming tool for Illumina NGS data",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino"

	output.dir   = TRIMREADS
	def Trimmomatic_FLAG1 = "SE -phred33 " +
                            Trimmomatic_THREADS 
	def Trimmomatic_FLAG2 = "ILLUMINACLIP:" + Trimmomatic_ADAPTER + ":2:30:10 SLIDINGWINDOW:4:15" + Trimmomatic_MINREADLEN
	
	transform(".fastq.gz") to ("_trimmed.fastq.gz") {
		exec """
			module load trimmomatic/${TRIMMOMATIC_VERSION} &&

                        java -jar ${TOOL_TRIMMOMATIC}/trimmomatic-0.36.jar $Trimmomatic_FLAG1 $input $output $Trimmomatic_FLAG2	
		""","Trimmomatic"
	}

}
