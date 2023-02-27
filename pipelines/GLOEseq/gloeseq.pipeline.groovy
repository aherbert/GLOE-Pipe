MODULE_FOLDER="./GLOEPipe/modules/"

load MODULE_FOLDER + "GLOEseq/essential.vars.groovy"
load MODULE_FOLDER + "GLOEseq/tool.locations.groovy"
load MODULE_FOLDER + "GLOEseq/tool.versions.groovy"

load MODULE_FOLDER + "GLOEseq/trimmomatic.vars.groovy"
load MODULE_FOLDER + "GLOEseq/trimmomatic.module.groovy"

load MODULE_FOLDER + "GLOEseq/cutadapt.vars.groovy"
load MODULE_FOLDER + "GLOEseq/cutadapt.module.groovy"

load MODULE_FOLDER + "GLOEseq/cutadapt_pe.vars.groovy"
load MODULE_FOLDER + "GLOEseq/cutadapt_pe.module.groovy"

load MODULE_FOLDER + "GLOEseq/fastqc.vars.groovy"
load MODULE_FOLDER + "GLOEseq/fastqc.module.groovy"

load MODULE_FOLDER + "GLOEseq/bowtie2.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bowtie2.module.groovy"

load MODULE_FOLDER + "GLOEseq/bowtie2_pe.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bowtie2_pe.module.groovy"

load MODULE_FOLDER + "GLOEseq/bamindexer.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bamindexer.module.groovy"

load MODULE_FOLDER + "GLOEseq/singlereads.module.groovy"
load MODULE_FOLDER + "GLOEseq/singlereads.vars.groovy"

load MODULE_FOLDER + "GLOEseq/bamqc.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bamqc.module.groovy"

load MODULE_FOLDER + "GLOEseq/markdups.vars.groovy"
load MODULE_FOLDER + "GLOEseq/markdups.module.groovy"

load MODULE_FOLDER + "GLOEseq/bam2bedD.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bam2bedD.module.groovy"

load MODULE_FOLDER + "GLOEseq/bam2bedI.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bam2bedI.module.groovy"

load MODULE_FOLDER + "GLOEseq/bedcoverage.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bedcoverage.module.groovy"

load MODULE_FOLDER + "GLOEseq/bed2bw.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bed2bw.module.groovy"

load MODULE_FOLDER + "GLOEseq/rfd.vars.groovy"
load MODULE_FOLDER + "GLOEseq/rfd.module.groovy"

load MODULE_FOLDER + "GLOEseq/macs2.vars.groovy"
load MODULE_FOLDER + "GLOEseq/macs2.module.groovy"

load MODULE_FOLDER + "GLOEseq/breaks_annotation.vars.groovy"                                                                                                                    
load MODULE_FOLDER + "GLOEseq/breaks_annotation.module.groovy" 

load MODULE_FOLDER + "GLOEseq/breaks_detected.vars.groovy"
load MODULE_FOLDER + "GLOEseq/breaks_detected.module.groovy"

load MODULE_FOLDER + "GLOEseq/bed2bz2.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bed2bz2.module.groovy"

load MODULE_FOLDER + "GLOEseq/multiqc.vars.groovy"
load MODULE_FOLDER + "GLOEseq/multiqc.module.groovy"

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"


//MAIN PIPELINE TASK - SINGLE END READS with Trimmomatic (INDIRECT mode - default)
run {
	"%.fastq.gz" * [ FastQC.using(subdir="raw"), Trimmomatic + [ FastQC.using(subdir:"trimmed"), bowtie2 + BAMindexer + BamQC + bam2bedI + [ bedcoverage, bed2bw + rfd ] + macs2 ] ] + [ breaks_annotation, breaks_detected ] + collectBpipeLogs + MultiQC
}

//MAIN PIPELINE TASK - SINGLE END READS with Cutadapt (INDIRECT mode - default) 
//run {
//
//	"%.fastq.gz" * [ 
//              FastQC.using(subdir="raw") , 
//              Cutadapt + 
//              [ 
//                  FastQC.using(subdir:"trimmed"), 
//                  bowtie2 + BAMindexer + BamQC + bam2bedI + [ 
//                                                              bedcoverage, 
//                                                              bed2bw + rfd,
//                                                            ] 
//              ] ] + collectBpipeLogs + MultiQC
//    }


//MAIN PIPELINE TASK - PAIRED END READS With Cutadapt (INDIRECT mode - default) 
//run {
// 
//         "%.fastq.gz" * [
//              FastQC.using(subdir="raw")
//                         ] +
//          "%.R*.fastq.gz" * [
//               Cutadapt_pe + 
//               [ FastQC.using(subdir:"trimmed"), bowtie2_pe + BAMindexer + MarkDups + BAMindexer + BamQC + SingleReads + BAMindexer + bam2bedI + [ bedcoverage, bed2bw + rfd ] ]
//                            ] + collectBpipeLogs + MultiQC
//     }


//MAIN PIPELINE TASK (DIRECT mode - optional) 
//run {
//	    "%.fastq.gz" * 
//	   [ FastQC.using(subdir="raw"), Trimmomatic + [ FastQC.using(subdir:"trimmed"), bowtie2 + BAMindexer + bam2bedD + [ bedcoverage, bed2bw + rfd ] + macs2 ] ]  + 
//	   [ breaks_annotation, breaks_detected ] + collectBpipeLogs + MultiQC
//}



//MAIN PIPELINE TASK (DIRECT mode - optional) USED FOR THE SAMPLES WITH UMIs 
//run {
//           "%.R*.fastq.gz" * [ AddUmiToFastq ] + "%.fastq.gz" *
//           [ FastQC.using(subdir="raw"), Trimmomatic + [ FastQC.using(subdir:"trimmed"), bowtie2 + BAMindexer + bam2bedD + [ bedcoverage, bed2bw ] + umidedup + BAMindexer + bam2bedD + [ bedcoverage, bed2bw ] ] ] + collectBpipeLogs + MultiQC
//}
