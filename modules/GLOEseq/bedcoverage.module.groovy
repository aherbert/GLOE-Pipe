bedcoverage = {
	doc title: "bedcoverage",
		desc:  "Convert Bed file to bigWig",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"
		
	output.dir=TRACKS_NORMALIZED
	BEDCOVERAGE_FLAGS = BEDCOVERAGE_CORES + " " + BEDCOVERAGE_OTHER
	
	transform(".bed") to (".normalized.bw") {
		exec """
                if [ ! -d ${SS_NORMALIZED} ]; then
                       mkdir -p ${SS_NORMALIZED};
                fi &&

                module load bedtools/${BEDTOOLS_VERSION} &&
                module load samtools/${SAMTOOLS_VERSION} &&
                module load deepTools/${DEEPTOOLS_VERSION} &&

                bedToBam -i $input -g $ESSENTIAL_BOWTIE_GENOME_INDEX > ${input.prefix}.bam  &&
                samtools index ${input.prefix}.bam &&

                bamCoverage $BEDCOVERAGE_FLAGS --bam ${input.prefix}.bam -o ${output.prefix}.bw &&
            
                bamCoverage $BEDCOVERAGE_FLAGS --bam ${input.prefix}.bam -o ${output.prefix}.fwd.bw --samFlagExclude 16 &&

                bamCoverage $BEDCOVERAGE_FLAGS --bam ${input.prefix}.bam -o ${output.prefix}.rev.bw --samFlagInclude 16 &&

                rm ${input.prefix}.bam &&
                rm ${input.prefix}.bam.bai &&
                mv ${output.prefix}.fwd.bw ${SS_NORMALIZED}/ &&
                mv ${output.prefix}.rev.bw ${SS_NORMALIZED}/


		""","bedcoverage"
	}
	forward input
}