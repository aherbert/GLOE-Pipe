MultiQC = {
	doc title: "MultiQC",
		desc:  "MultiQC is a reporting tool that parses summary statistics from results and log files generated by other bioinformatics tools",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.9.5.slurm",
		author: "Giuseppe Petrosino"
    
	output.dir   = MultiQC_OUTDIR
	def MultiQC_FLAGS = "--ignore .bpipe"
	
    produce("multiqc_report.html") {
		exec """

                module load MultiQC/${MultiQC_VERSION} &&
                multiqc $QC $LOGS $MultiQC_FLAGS -o $output.dir &&

                if [ -e "${QC}/fastqc/trimmed" ]; then
                       mkdir -p ${output.dir}/raw;
            	       multiqc $QC/fastqc/raw $MultiQC_FLAGS -o ${output.dir}/raw;
                fi;  
  
		""","MultiQC"
	}
}
