SingleReads = {
	doc title: "SingleReads",
                desc:  "Select only R1 mapped in correct orientation and within insert size. https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/",
                constraints: "Samtools multithreaded version expected (>=1.2)",
                bpipe_version: "tested with bpipe 0.9.9.5.slurm",
                author: "Giuseppe Petrosino"

	output.dir=MAPPED

    def SAMTOOLS_FLAGS = SAMTOOLS_THREADS

	transform(".bam") to (".sr.bam") {
          exec """
            module load samtools/${SAMTOOLS_VERSION} &&
	
            samtools view $SAMTOOLS_FLAGS -hbf 99 $input > ${output.prefix}_bowtie2_q30_99.bam &&
            samtools view $SAMTOOLS_FLAGS -hbf 83 $input > ${output.prefix}_bowtie2_q30_83.bam &&
            samtools merge $SAMTOOLS_FLAGS $output ${output.prefix}_bowtie2_q30_99.bam ${output.prefix}_bowtie2_q30_83.bam &&
            rm ${output.prefix}_bowtie2_q30_99.bam ${output.prefix}_bowtie2_q30_83.bam

          ""","SingleReads"
	}
}
