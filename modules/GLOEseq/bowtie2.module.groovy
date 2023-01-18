bowtie2 = {
    doc title: "Bowtie2 SE alignment",
        desc:  "Align single end reads",
        constraints: "Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.9.5.slurm",
        author: "Giuseppe Petrosino"

    output.dir = MAPPED

    def BOWTIE_FLAGS = " "  +
                       BOWTIE_THREADS  + " "


    def SAMTOOLS_VIEW_FLAG1 = "-bhSu"
    def SAMTOOLS_SORT_FLAGS = "-O bam " + BOWTIE_SAMTOOLS_THREADS
    def SAMTOOLS_VIEW_FLAG2  = "-q 30 -bhu"

    transform(".fastq.gz") to (".bam") {
        exec """
            module load bowtie2/${BOWTIE2_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi                                       &&  
            bowtie2 $BOWTIE_FLAGS -x $BOWTIE_REF $input | samtools view $SAMTOOLS_VIEW_FLAG1 - | samtools sort $SAMTOOLS_SORT_FLAGS -T $TMPDIR/\$(basename $output.prefix)_bowtie2_sorted - > ${output.prefix}_bowtie2_sorted.bam &&
            samtools view $SAMTOOLS_VIEW_FLAG2 ${output.prefix}_bowtie2_sorted.bam | samtools sort $BOWTIE_SAMTOOLS_THREADS -T $TMPDIR/\$(basename $output.prefix) -o $output &&
            rm ${output.prefix}_bowtie2_sorted.bam
        ""","bowtie2"
    }
}
