bam2bedI5ends = {
	doc title: "bam2bedI5ends",
		desc:  "Convert BAM file to bed and the perl code identifies the first base before the 5' of each read on the same strand (indirect mode)",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"

	output.dir=BED
	
	def PBZIP2_FLAG = " "  + PBZIP2_THREADS + " "

	transform(".bam") to (".bed") {
		exec """

			module load kentUtils/${KENTUTILS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load bedtools/${BEDTOOLS_VERSION} &&

			bamToBed -i $input > ${output.prefix}.reads.bed &&
			perl ${TOOL_GLOEseq}/indirect_mode_5ends.pl ${output.prefix}.reads.bed | sort -k1,1 -k2,2n -k 6 - |  awk '(\$2 >= 0)' > ${output.prefix}.sites.bed &&
			awk  '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" "0" "\\t" \$6}' ${output.prefix}.sites.bed > $output.bed &&
			grep "+" $output.bed > ${output.prefix}.for.bed &&
			grep "-" $output.bed > ${output.prefix}.rev.bed &&

			rm ${output.prefix}.reads.bed &&
			rm ${output.prefix}.sites.bed

		""","bam2bedI5ends"
	}
}
