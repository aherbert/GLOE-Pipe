##############################################################################
##
## What: Breaks_Annotation.R
## Who: Giuseppe Petrosino
## When: 
##
## Script to annotate GLOE-seq breaks using result from MACS2
##
## Args: 
## -----
## breakData=                                         # path to the bed file result from MACS2
## transcriptType=Bioconductor                        # define the transcript annotation type for the analysis
## transcriptDb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene # either GTF file or transcript annotation database as transcript annotation database
## orgDb=org.Sc.sgd.db                                # optionally genome wide annotation 
## regionTSS=200                                      # TSS region parameter from -200 to +200
## out=                                               # output directory
############################################################################
options(stringsAsFactors=FALSE)
library(ChIPseeker)
library(openxlsx)
library(Cairo)
library(regioneR)

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
 
   if(length(i <- grep(string,args,fixed=T)) == 1)
   return(do.call(convert,list(gsub(string,"",args[i]))))
   if(!is.null(default)) default else do.call(convert,list(NA))
  
}

args <- commandArgs(T)
breakData <- parseArgs(args,"breakData=") # .bed result from MACS2
transcriptType <- parseArgs(args, "transcriptType=", "Bioconductor") # transcript annotation type
transcriptDb <- parseArgs(args, "transcriptDb=", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene") # transcript annotation database
orgDb <- parseArgs(args, "orgDb=", "org.Sc.sgd.db") # genome wide annotation
regionTSS <- parseArgs(args, "regionTSS=", 200, "as.numeric") # TSS region parameter 
out <- parseArgs(args,"out=", "Breaks_Annotation") # output directory
index <- parseArgs(args, "index=")
covMax <- parseArgs(args, "covMax=", 100000, "as.numeric") # Max elements for the coverage plot

runstr <- paste0("Call with: Rscript Breaks_Annotation.R [breakData=",breakData,"] [transcriptType=",transcriptType,"] [transcriptDb=",transcriptDb,"] [orgDb=",orgDb,"] [regionTSS=",regionTSS,"] [index=",index,"] [out=",out,"]")
cat(runstr)
if (!is.numeric(regionTSS)) stop("regionTSS not numeric. Run with:\n",runstr)

breakFiles <-list.files(breakData,pattern=".bed", full.names = TRUE) # list of the full path of the .bed file 
breaks <- lapply(breakFiles, toGRanges) # read all the bed files using 'toGRanges' function
breaks <- breaks[lapply(breaks,length)>0] # select empty bed files
chr.size <- read.table(file = index, sep="\t")

if(transcriptType!="Bioconductor"){ # check the input format for the transcript annotation
   txdb <- makeTxDbFromGFF(transcriptDb, format="gtf") # if the input format is gtf file, then this file will be used to create a TxDb object
} else {
   library(transcriptDb, character.only = TRUE) # if the input format is bioconductor, then the transcript annoation library will be used 
   txdb <- eval(parse(text=transcriptDb))
  
}

if(orgDb!=""){ # check if genome wide annotation should be used
  breakAnno <- lapply(breaks, annotatePeak, TxDb=txdb, tssRegion=c(-regionTSS,regionTSS), annoDb= orgDb, sameStrand=TRUE)
} else {
  breakAnno <- lapply(breaks, annotatePeak, TxDb=txdb, tssRegion=c(-regionTSS,regionTSS), sameStrand=TRUE)  
}

breakFiles  <- breakFiles[which(lengths(breaks)>0)] # remove the filenames of empty bed files
filename <- strsplit(basename(breakFiles), "_macs2_summits.bed") # take the filenames and put it as names for the plots
names(breakAnno) <- filename

# create barplot showing the feature distribution
CairoPNG(file=paste0(out, "/GLOEseq_Feature_Distribution_Barplot.png"), width = 700, height = 500)
plot(plotAnnoBar(breakAnno))
dev.off()

# create barplot showing the feature distribution related to TSS
CairoPNG(file=paste0(out, "/GLOEseq_Feature_Distribution_Related_to_TSS_Barplot.png"), width = 700, height = 500)
plot(plotDistToTSS(breakAnno))
dev.off()

# create upsetplot 
for(i in 1:length(breakAnno)){
  CairoPNG(file=paste0(out, "/", filename[[i]], "_GLOEseq_UpSetplot.png"), width = 700, height = 500)
  plot(upsetplot(breakAnno[[i]]))
  dev.off()
}

# create Breaks coverage plot
# Note: This plot takes too long (i.e. days) when the number of annotations is large.
# The number to plot is limited to a random subset.
for(i in 1:length(breakAnno)){
  CairoPNG(file=paste0(out, "/", filename[[i]], "_GLOEseq_Breaks_Coverageplot.png"), width = 700, height = 500)
  b <- breaks[[i]]
  n <- length(b)
  if (n > covMax) {
    b <- b[sample(n, covMax)]
  }
  plot(covplot(b,weightCol="score", title= "GLOEseq breaks over Chromosomes", ylab=("-log10(qvalue))")))
  dev.off()
}

# create % of breaks per chromosome plot
# for(i in 1:length(breakAnno)){
# CairoPNG(file=paste0(out, "/", filename[[i]], "_GLOEseq_Breaks_per_Chromosomeplot.png"), width = 700, height = 500)
# chr<-as.data.frame(seqnames(breaks[[i]]))
# freq<-as.data.frame(table(chr$value))	
  
  #Merge two tables
# tab <- merge(chr.size, freq, by.x="V1", by.y="Var1", all = TRUE)
# tab[is.na(tab)] <- 0

  #Calculate percentage of breaks per chromosome 
# tab$ratio=(100/tab$V2)*tab$Freq
  
#  x1  = factor(tab$V1, levels=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM"))
#  par(mar=c(3.5, 5, 2,1))
#  plot(tab$ratio ~ x1, type="l", las=1, ylab="Percentage of breaks per chromosome", cex.axis=0.7, cex.lab=0.8)

#  dev.off()
#}

# create xls output contains the breaks annotation
outputData <- lapply(breakAnno, as.data.frame)
# Note: When the number of annotation is large saving as xlsx is not possible
# openxlsx library fails when joining the shared strings object that reuses strings across sheets
# writexl library fails as the output is above 2^20 lines
# So we save as CSV
#write.xlsx(outputData, file=paste0(out, "/Breaks_Annotation.xlsx"))
for (i in 1:length(outputData)) {
  write.table(outputData[[i]], paste0(out, "/", filename[[i]], ".csv"), quote=FALSE, sep="\t", row.names=FALSE)
}

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GLOEseq_Breaks_Annotation_session_info.txt", sep=""))
save(breakFiles,breaks,breakAnno,filename,file=paste0(out,"/Breaks_Annotation.RData"))
