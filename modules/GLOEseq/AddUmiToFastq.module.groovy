AddUmiToFastq = {
    doc title: "Adds UMI to the fastq header",
        desc:  "adds UMI of the second read in GLOE-Seq samples to the fastq header using umitools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Giuseppe Petrosino"

    output.dir = UMIREADS

    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /.R1.fastq.gz/).replaceFirst("")


	produce(OUTPUTFILE + ".umibarcode.fastq.gz") {
		exec """

                                 if [ ! -d ${UMIREADS} ]; then
				             mkdir -p ${UMIREADS};
			         fi &&

			         if [ ! -d ${TMP} ]; then
				             mkdir -p ${TMP};
			         fi &&


				module load umitools/${UMITOOLS_VERSION} &&


                umi_tools extract --bc-pattern=NNNNNNNNNN --stdin $input2 --stdout ${TMP}/\$(basename ${input2.prefix}).barcode.fastq.gz --read2-in $input1 --read2-out=$output &&

                rm ${TMP}/\$(basename ${input2.prefix}).barcode.fastq.gz 
        ""","AddUmiToFastq"
    }
}

