bedcoverage3ends = {
	doc title: "bedcoverage3ends",
		desc:  "Convert Bed file to RE normalized bigWig",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino"
		
	output.dir=TRACKS_RE_NORMALIZED
	BEDCOVERAGE_FLAGS = BEDCOVERAGE_CORES + " " + BEDCOVERAGE_OTHER

    breaks_signal_3END_FLAGS = BREAKS_3END + " " + QPCR_CUT_EFF + " " + BREAKS_OUTDIR_3END + " "
	
	transform(".bed") to (".RE.normalized.bw") {
		exec """
                if [ ! -d ${ALPHA} ]; then
                       mkdir -p ${RESULTS} &&
                       mkdir -p ${ALPHA};
                fi &&

                if [ ! -d ${ALPHA}/3end ]; then
                       mkdir -p ${ALPHA}/3end;
                fi &&

                if [ ! -d ${SS_RE_NORMALIZED} ]; then
                       mkdir -p ${SS_RE_NORMALIZED};
                fi &&

                module load bedtools/${BEDTOOLS_VERSION} &&
                module load samtools/${SAMTOOLS_VERSION} &&
                module load deepTools/${DEEPTOOLS_VERSION} &&
                module load R/${R_VERSION} &&

                Rscript ${TOOL_BREAKS_DETECTED}/Breaks_signal_3ends.R Bed=$input $breaks_signal_3END_FLAGS &&

                bedToBam -i $input -g $ESSENTIAL_BOWTIE_GENOME_INDEX > ${input.prefix}.bam  &&
                samtools index ${input.prefix}.bam &&

                ALPHA=\$( cat ${ALPHA}/3end/\$(basename ${input.prefix}).alpha.txt) &&
                SCALE=\$(echo "1/\$ALPHA" | bc -l) &&
                
                bamCoverage $BEDCOVERAGE_FLAGS --scaleFactor \${SCALE} --bam ${input.prefix}.bam -o ${output.prefix}.bw &&
            
                bamCoverage $BEDCOVERAGE_FLAGS --scaleFactor \${SCALE} --bam ${input.prefix}.bam -o ${output.prefix}.fwd.bw --samFlagExclude 16 &&

                bamCoverage $BEDCOVERAGE_FLAGS --scaleFactor \${SCALE} --bam ${input.prefix}.bam -o ${output.prefix}.rev.bw --samFlagInclude 16 &&

                rm ${input.prefix}.bam &&
                rm ${input.prefix}.bam.bai &&
                mv ${output.prefix}.fwd.bw ${SS_RE_NORMALIZED}/ &&
                mv ${output.prefix}.rev.bw ${SS_RE_NORMALIZED}/


		""","bedcoverage3ends"
	}
	forward input
}