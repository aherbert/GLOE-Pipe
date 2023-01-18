bowtie2_pe = {
    doc title: "Bowtie2 PE alignment",
        desc:  "Align paired end reads and select only R1 mapped in correct orientation and within insert size. https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/",
        constraints: "Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.9.5.slurm",
        author: "Giuseppe Petrosino"

    output.dir = MAPPED


    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /.fastq.gz/).replaceFirst("")

    def BOWTIE2_FLAGS = "-q "  +
                       BOWTIE2_QUALS    + " " + 
                       BOWTIE2_MM_SEED  + " " + 
                       BOWTIE2_INSERT   + " " + 
                       BOWTIE2_THREADS  + " " + 
                       BOWTIE2_EXTRA


    def SAMTOOLS_VIEW_FLAG1 = "-bhSu"
    def SAMTOOLS_SORT_FLAGS = "-O bam " + BOWTIE_SAMTOOLS_THREADS
    def SAMTOOLS_VIEW_FLAG2  = "-q 30 -bhu"

    produce(OUTPUTFILE + ".bam") {
        exec """
            module load bowtie2/${BOWTIE2_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi                                       &&  

            bowtie2 $BOWTIE2_FLAGS $BOWTIE2_REF -1 $input1 -2 $input2 | samtools view $SAMTOOLS_VIEW_FLAG1 - | samtools sort $SAMTOOLS_SORT_FLAGS -T $TMPDIR/\$(basename $output.prefix)_bowtie2 - > ${output.prefix}_bowtie2.bam &&
            samtools view $SAMTOOLS_VIEW_FLAG2 ${output.prefix}_bowtie2.bam | samtools sort $BOWTIE_SAMTOOLS_THREADS -T $TMPDIR/\$(basename $output.prefix)_bowtie2_q30 -o ${output.prefix}_bowtie2_q30.bam &&
            samtools index ${output.prefix}_bowtie2_q30.bam &&
            samtools view -hbf 99 ${output.prefix}_bowtie2_q30.bam > ${output.prefix}_bowtie2_q30_99.bam &&
            samtools view -hbf 83 ${output.prefix}_bowtie2_q30.bam > ${output.prefix}_bowtie2_q30_83.bam &&
            samtools merge $output ${output.prefix}_bowtie2_q30_99.bam ${output.prefix}_bowtie2_q30_83.bam &&
            rm ${output.prefix}_bowtie2.bam ${output.prefix}_bowtie2_q30_99.bam ${output.prefix}_bowtie2_q30_83.bam
 
        ""","bowtie2_pe"
    }
}
