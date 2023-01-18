rfd = {
	doc title: "rfd",
		desc:  "Calculate the proportions rightward- and leftward- moving forks within each window (1kb) [RFD=(C-W)/(C+W)] and generate RFD track",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino"
		
		output.dir=RFD
		RFD_FLAGS = RFD_CORES
	
	transform(".bw") to (".1kb.rfd.bw") {
		exec """
                if [ ! -d ${RFD} ]; then
                       mkdir -p ${RFD};
                fi &&
                
                module load deepTools/${DEEPTOOLS_VERSION} &&

                FWD=${SS}/\$(basename ${input.prefix}).fwd.bw &&
                REV=${SS}/\$(basename ${input.prefix}).rev.bw &&

                bigwigCompare $RFD_FLAGS -b1 $REV -b2 $FWD -bs 1000 --operation subtract -of bigwig -o ${output.prefix}.rev-fwd.bw &&
                
                bigwigCompare $RFD_FLAGS -b1 $REV -b2 $FWD -bs 1000 --operation add -of bigwig -o ${output.prefix}.rev+fwd.bw &&
                  
                bigwigCompare $RFD_FLAGS -b1 ${output.prefix}.rev-fwd.bw -b2 ${output.prefix}.rev+fwd.bw -bs 1000 --pseudocount 0 --operation ratio -of bigwig -o ${output.prefix}.bw &&
 
                rm ${output.prefix}.rev-fwd.bw ${output.prefix}.rev+fwd.bw
      

		""","rfd"
	}
}
