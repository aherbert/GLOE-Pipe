MODULE_FOLDER="./GLOEPipe/modules/"

load MODULE_FOLDER + "GLOEseq/essential.vars.groovy"
load MODULE_FOLDER + "GLOEseq/tool.locations.groovy"
load MODULE_FOLDER + "GLOEseq/tool.versions.groovy"

load MODULE_FOLDER + "GLOEseq/trimmomatic.vars.groovy"
load MODULE_FOLDER + "GLOEseq/trimmomatic.module.groovy"

load MODULE_FOLDER + "GLOEseq/fastqc.vars.groovy"
load MODULE_FOLDER + "GLOEseq/fastqc.module.groovy"

load MODULE_FOLDER + "GLOEseq/bowtie2.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bowtie2.module.groovy"

load MODULE_FOLDER + "GLOEseq/bamindexer.vars.groovy"
load MODULE_FOLDER + "GLOEseq/bamindexer.module.groovy"

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

load MODULE_FOLDER + "miscellaneous/collectbpipes.module.2.groovy"


//MAIN PIPELINE TASK (INDIRECT mode - default)
run {
            "%.fastq.gz" *
           [ FastQC, Trimmomatic + [ FastQC, bowtie2 + BAMindexer + bam2bedI + [ bedcoverage, bed2bw + rfd, macs2 ] ] ]  +
           [ breaks_annotation, breaks_detected ] + collectBpipeLogs
}


//MAIN PIPELINE TASK (DIRECT mode - optional) 
//run {
//	    "%.fastq.gz" * 
//	   [ FastQC, Trimmomatic + [ FastQC, bowtie2 + BAMindexer + bam2bedD + [ bedcoverage, bed2bw + rfd, macs2 ] ] ]  + 
//	   [ breaks_annotation, breaks_detected ] + collectBpipeLogs 
//}