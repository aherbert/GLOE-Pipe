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

load MODULE_FOLDER + "GLOEseq/cutadapt_5ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/cutadapt_5ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/cutadapt_3ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/cutadapt_3ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/fastqc.vars.groovy"
load MODULE_FOLDER + "GLOEseq/fastqc.module.groovy"

load MODULE_FOLDER + "GLOEseq/fastqscreen.vars.groovy"
load MODULE_FOLDER + "GLOEseq/fastqscreen.module.groovy"

load MODULE_FOLDER + "GLOEseq/merge_fastq.vars.groovy"
load MODULE_FOLDER + "GLOEseq/merge_fastq.module.groovy"

load MODULE_FOLDER + "GLOEseq/barcode5ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/barcode5ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/barcode3ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/barcode3ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/AddUmiToFastq.vars.groovy"
load MODULE_FOLDER + "GLOEseq/AddUmiToFastq.module.groovy"

load MODULE_FOLDER + "GLOEseq/bowtie2.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bowtie2.module.groovy"

load MODULE_FOLDER + "GLOEseq/bowtie2_pe.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bowtie2_pe.module.groovy"

load MODULE_FOLDER + "GLOEseq/umidedup.vars.groovy"
load MODULE_FOLDER + "GLOEseq/umidedup.module.groovy"

load MODULE_FOLDER + "GLOEseq/bamindexer.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bamindexer.module.groovy"

load MODULE_FOLDER + "GLOEseq/singlereads5ends.module.groovy"
load MODULE_FOLDER + "GLOEseq/singlereads5ends.vars.groovy"

load MODULE_FOLDER + "GLOEseq/singlereads3ends.module.groovy"
load MODULE_FOLDER + "GLOEseq/singlereads3ends.vars.groovy"

load MODULE_FOLDER + "GLOEseq/bamqc.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bamqc.module.groovy"

load MODULE_FOLDER + "GLOEseq/markdups.vars.groovy"
load MODULE_FOLDER + "GLOEseq/markdups.module.groovy"

load MODULE_FOLDER + "GLOEseq/bam2bedD.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bam2bedD.module.groovy"

load MODULE_FOLDER + "GLOEseq/bam2bedI.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bam2bedI.module.groovy"

load MODULE_FOLDER + "GLOEseq/bam2bedD5ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bam2bedD5ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/bam2bedI5ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bam2bedI5ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/bedcoverage.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bedcoverage.module.groovy"

load MODULE_FOLDER + "GLOEseq/bedcoverage3ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bedcoverage3ends.module.groovy"

load MODULE_FOLDER + "GLOEseq/bedcoverage5ends.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bedcoverage5ends.module.groovy"

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

load MODULE_FOLDER + "GLOEseq/shinyreports.vars.groovy"
load MODULE_FOLDER + "GLOEseq/shinyreports.module.groovy"

//MAIN PIPELINE TASK - SINGLE END READS with Trimmomatic (INDIRECT mode - default)
run {
	"%.fastq.gz" * [ FastQC.using(subdir="raw"), Trimmomatic + [ FastQC.using(subdir:"trimmed"), bowtie2 + BAMindexer + BamQC + bam2bedI + [ bedcoverage, bed2bw + rfd ] ] ] + macs2 +  [ breaks_annotation ] + collectBpipeLogs + MultiQC
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
//               [ FastQC.using(subdir:"trimmed"), bowtie2_pe + BAMindexer + MarkDups + BAMindexer + BamQC + SingleReads3ends + BAMindexer + bam2bedI + [ bedcoverage, bed2bw + rfd ] ]
//                            ] + collectBpipeLogs + MultiQC
//     }

//MAIN PIPELINE TASK - PAIRED END READS With Cutadapt (DIRECT mode) 
//run {
//
//         "%.fastq.gz" * [
//              FastQC.using(subdir="raw")
//                         ] +
//          "%.R*.fastq.gz" * [
//               Cutadapt_pe + [ FastQC.using(subdir:"trimmed"), bowtie2_pe + BAMindexer + BamQC ]
//                         ] +
//           "%.bam" * [
//               SingleReads3ends + [ BAMindexer + bam2bedD + [ bedcoverage3ends + bedcoverage, bed2bw ] + macs2 ]
//               ] + [ breaks_detected ] + collectBpipeLogs + MultiQC
//    }

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


////MAIN PIPELINE TASK - GLOE-seq v2 (DIRECT mode)
//run {
//        "%.R*.fastq.gz" * [ 
//                         FastqScreen, Merge_fastq 
//                          ] +
//	    "%.R*.m.fastq.gz" * 
//	   [ FastQC.using(subdir="raw"),
//		Barcode5ends + [ Cutadapt_5ends + [ FastQC.using(subdir:"5ends"), bowtie2_pe + BAMindexer + BamQC + SingleReads5ends + BAMindexer + bam2bedD5ends + [ bedcoverage5ends + bedcoverage, bed2bw ] ] ], 
//		Barcode3ends + [ Cutadapt_3ends + [ FastQC.using(subdir:"3ends"), bowtie2_pe + BAMindexer + BamQC + SingleReads3ends + BAMindexer + bam2bedD + [ bedcoverage3ends + bedcoverage, bed2bw ] ] ] ] + 
//	   collectBpipeLogs + MultiQC + shinyReports
//}
