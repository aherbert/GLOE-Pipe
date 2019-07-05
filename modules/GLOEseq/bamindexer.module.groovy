//rule for task BAMindexer from catalog NGS, version 1
//desc: Samtools index a bam file
BAMindexer = {
	doc title: "BAMindexer",
		desc:  "Call samtools to index a bam file",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"

	output.dir = MAPPED

	transform(".bam") to(".bam.bai") {
		exec """
			module load samtools/${SAMTOOLS_VERSION} &&
		        samtools flagstat $input > ${output.prefix}.flagstat &&
			samtools index $input
		""","BAMindexer"
	}

	forward input
}


