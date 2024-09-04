bed2bw = {
	doc title: "bed2bw",
		desc:  "Convert Bed file to bigWig",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino"

	output.dir=TRACKS
	
	transform(".bed") to (".bw") {
		exec """
			         if [ ! -d ${SS} ]; then
				             mkdir -p ${SS};
			         fi &&
			         
			         if [ ! -d ${TMP} ]; then
				             mkdir -p ${TMP};
			         fi &&

			module load bedtools/${BEDTOOLS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load kentUtils/${KENTUTILS_VERSION} &&


			CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
			samtools idxstats ${MAPPED}/\$(basename ${input.prefix}).bam | cut -f1-2 | grep -v '^\\*' > ${CHRSIZES} &&
			sort -k1,1 -k2,2n $input -o ${output.prefix}.sorted.bed &&
			genomeCoverageBed -bg -i ${output.prefix}.sorted.bed -trackline -g ${CHRSIZES} > ${output.prefix}.bedgraph &&
			genomeCoverageBed -bg -i ${output.prefix}.sorted.bed -trackline -g ${CHRSIZES} -strand + > ${output.prefix}.fwd.bedgraph &&
			genomeCoverageBed -bg -i ${output.prefix}.sorted.bed -trackline -g ${CHRSIZES} -strand - > ${output.prefix}.rev.bedgraph &&
			bedGraphToBigWig ${output.prefix}.bedgraph ${CHRSIZES} $output &&
			bedGraphToBigWig ${output.prefix}.fwd.bedgraph ${CHRSIZES} ${output.prefix}.fwd.bw &&
			bedGraphToBigWig ${output.prefix}.rev.bedgraph ${CHRSIZES} ${output.prefix}.rev.bw &&


			rm ${output.prefix}.sorted.bed &&
			rm ${output.prefix}.bedgraph &&
			rm ${output.prefix}.fwd.bedgraph &&
			rm ${output.prefix}.rev.bedgraph &&
			mv ${output.prefix}.fwd.bw ${SS}/ &&
			mv ${output.prefix}.rev.bw ${SS}/ &&
			rm ${CHRSIZES}

		""","bed2bw"
	}

            // to generate other intermediate files

            //       if [ ! -d ${BG} ]; then
            //           	mkdir -p ${BG};      
            //       fi &&

            // mv ${output.prefix}.bedgraph ${BG}/ &&

}
