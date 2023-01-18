MarkDups = {
	doc title: "MarkDups",
                desc:  "Call picard tools to mark without removing duplicated reads from a bam file",
                constraints: "Picard tools",
                bpipe_version: "tested with bpipe 0.9.9.3.slurm",
                author: "Giuseppe Petrosino"

	output.dir=MAPPED
	def JAVA_FLAGS  = "-Xmx" + MARKDUPS_MAXMEM + "m"
	def MARKDUPS_FLAGS  = "REMOVE_DUPLICATES=FALSE ASSUME_SORTED=TRUE"

	transform(".bam") to (".dupmarked.bam") {
		exec """
			module load jdk/${JAVA_VERSION} &&
			module load picard/${PICARD_VERSION} &&
			
			java $JAVA_FLAGS -jar ${TOOL_PICARD}/picard.jar MarkDuplicates $MARKDUPS_FLAGS INPUT=$input OUTPUT=$output METRICS_FILE=${input.prefix}_dupmetrics.tsv
		""","MarkDups"
	}
}
