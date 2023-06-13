Merge_fastq = {
    doc title: "Merge Fastq files",
        desc:  "adds index to the R1 fastq",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Giuseppe Petrosino"

    output.dir = MERGEDREADS


    // extract the base of the input file name (w/o the directory, w/o fastq.gz)
       def OUTPUTFILE = input1
       int path_index = OUTPUTFILE.lastIndexOf("/")
       OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
       OUTPUTFILE = (OUTPUTFILE =~ /.R1.fastq.gz/).replaceFirst("")


    produce(OUTPUTFILE + "*.merged.fastq.gz") {
		exec """

                        if [ -n "\$SLURM_JOBID" ]; then
                                export TMPDIR=/jobdir/\${SLURM_JOBID};
                        fi &&



               paste -d "" <(zcat $input3) <(zcat $input1 | sed 's/^+.*//' | sed 's/^@.*//') | gzip > ${MERGEDREADS}/${OUTPUTFILE}.R1.merged.fastq.gz &&
               cp $input2 ${MERGEDREADS}/${OUTPUTFILE}.R2.merged.fastq.gz

        ""","Merge_fastq"
    }

}
