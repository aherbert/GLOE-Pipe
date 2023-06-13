Barcode3ends = {
    doc title: "Barcode 3'ends",
        desc:  "extracts reads based to the index for the 3'ends in the R1 fastq",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Giuseppe Petrosino"

    output.dir = UMIREADS

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_PRUNED = SAMPLENAME_BASE.replace(".R1.merged", "") // delete .R1 in combined log file of pe design


	transform("*.fastq.gz") to (".end3.fastq.gz") {
		exec """

             if [ ! -d ${UMIREADS} ]; then
		             mkdir -p ${UMIREADS};
	         fi &&

		     if [ ! -d ${TMP} ]; then
		             mkdir -p ${TMP};
		     fi &&


		module load umitools/${UMITOOLS_VERSION} &&

        umi_tools extract --extract-method=regex --bc-pattern='^(?P<discard_1>[AG][CT][AG][CT][AG][CT]){s<=1}(?P<umi_1>.{6})' -I $input1 --stdout ${UMIREADS}/${SAMPLENAME_BASE_PRUNED}.R1.merged.end3.fastq.gz --read2-in $input2 --read2-out=${UMIREADS}/${SAMPLENAME_BASE_PRUNED}.R2.merged.end3.fastq.gz


        ""","Barcode3ends"
    }
}
