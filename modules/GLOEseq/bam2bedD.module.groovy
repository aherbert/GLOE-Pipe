bam2bedD = {
	doc title: "bam2bedD",
		desc:  "Convert BAM file to bed and the perl code identifies the first base  each read on the opposite strand (direct mode)",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.9.3.slurm",
		author: "Giuseppe Petrosino"

	output.dir=BED
	
	def PBZIP2_FLAG = " "  + PBZIP2_THREADS + " "

	transform(".bam") to (".bed") {
		exec """
			if [ ! -d ${READS} ]; then
				 mkdir -p ${READS};
			fi &&

                        if [ ! -d ${SITES} ]; then
	                         mkdir -p ${SITES};      
                        fi &&

                        if [ ! -d ${COUNTS} ]; then
				mkdir -p ${COUNTS};      
                        fi &&

                        if [ ! -d ${TPM} ]; then
	                        mkdir -p ${TPM};      
                        fi &&

			module load kentUtils/${KENTUTILS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load bedtools/${BEDTOOLS_VERSION} &&

			bamToBed -i $input > ${output.prefix}.reads.bed &&
			perl ${TOOL_GLOEseq}/direct_mode.pl ${output.prefix}.reads.bed | sort -k1,1 -k2,2n -k 6 - |  awk '(\$2 >= 0)' > ${output.prefix}.sites.bed &&
			awk  '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" "0" "\\t" \$6}' ${output.prefix}.sites.bed > $output.bed &&
			grep "+" $output.bed > ${output.prefix}.for.bed &&
			grep "-" $output.bed > ${output.prefix}.rev.bed &&
			
			bedtools groupby -i ${output.prefix}.sites.bed -g 1,3,6 -c 1 -full -o count > ${output.prefix}.counts.bed &&
			export TOTALFPE=`wc -l ${output.prefix}.sites.bed | awk '{print \$1}'` &&
			cat ${output.prefix}.counts.bed | perl -alne '\$F[3]=(\$F[6]/\$ENV{TOTALFPE})*1e6; print join "\\t", @F;' > ${output.prefix}.tpm.bed &&
			
			pbzip2 $PBZIP2_FLAG ${output.prefix}.reads.bed &&
			pbzip2 $PBZIP2_FLAG ${output.prefix}.sites.bed &&
			pbzip2 $PBZIP2_FLAG ${output.prefix}.counts.bed &&
			pbzip2 $PBZIP2_FLAG ${output.prefix}.tpm.bed &&

			mv ${output.prefix}.reads.bed.bz2 ${READS}/ &&
			mv ${output.prefix}.sites.bed.bz2 ${SITES}/ &&                   
			mv ${output.prefix}.counts.bed.bz2 ${COUNTS}/ &&
			mv ${output.prefix}.tpm.bed.bz2 ${TPM}/ 


		""","bam2bedD"
	}
}