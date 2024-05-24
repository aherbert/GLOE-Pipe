#############################################################################
##
## What: Breaks_signal_3ends.R
## Who: Giuseppe Petrosino
## When: 2024-05-03
##
## Script to count the number of reads for the breaks on the 3' ends,
## include qPCR percentage cut efficiency and calculate REs signal (alpha)
##                    
############################################################################
options(stringsAsFactors=FALSE)
library(GenomicRanges)

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
 
   if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)

file_a <- parseArgs(args,"ExpectedBreaks3End=") 
file_b <- parseArgs(args,"Bed=") 
qPCR <- parseArgs(args,"qPCR=")
out <- parseArgs(args,"out=")

runstr <- paste0("Call with: Rscript Breaks_signal_3ends.R [ExpectedBreaks3End=",file_a,"] [Bed=",file_b,"] [qPCR=",qPCR,"] [out=",out,"]")
cat(runstr)

qPCR <- read.table(qPCR, sep="\t", header=T)

x <- strsplit(basename(file_b), "_")[[1]][1:5]
x <- paste(x[1], x[2], x[3] , x[4], x[5],  sep="_")
qPCR_value <- qPCR[grep(x, qPCR$Samples),]$qPCR/100

# Read BED files
gr_a <- read.table(file_a, header = FALSE, col.names = c("seqnames", "start", "end", "name", "score", "strand"), sep = "\t")
gr_b <- read.table(file_b, header = FALSE, col.names = c("seqnames", "start", "end", "name", "score", "strand"), sep = "\t")

# Create GRanges objects
ranges_a <- GRanges(seqnames = gr_a$seqnames, ranges = IRanges(start = gr_a$start, end = gr_a$end), strand = gr_a$strand)
ranges_b <- GRanges(seqnames = gr_b$seqnames, ranges = IRanges(start = gr_b$start, end = gr_b$end), strand = gr_b$strand)

# Count overalaps 
ranges_a$count <- countOverlaps(ranges_a, ranges_b, ignore.strand = FALSE, minoverlap=0, type="start")

# Calculate Alpha for normalized the tracks
res <- quantile(ranges_a$count, probs = 0.015)
breaks_norm <- ranges_a[ranges_a$count > res[1],]
alpha <- (sum(ranges_a$count)/length(breaks_norm))/qPCR_value
#alpha <- (sum(ranges_a$count)/length(ranges_a[ranges_a$count > 1,]))/qPCR_value

# Total # reads (no chrM)
total_reads_noChrM <- dim(gr_b[gr_b$seqnames != "chrM",])[1]
# % reads at quantification sites
perc_reads_quant_sites <- (100/total_reads_noChrM)*sum(ranges_a$count)
# Breaks per genome
breaks_genome <- (total_reads_noChrM - sum(ranges_a$count)) / alpha
# reads per quant site / Total # reads at quant sites per sample * 100
ranges_a$norm_counts <- (ranges_a$count/sum(ranges_a$count))*100

y <- gsub(".bed", "", basename(file_b))

# Output files
output_file1 <- paste0(out, "/", y, ".alpha.txt")
output_file2 <- paste0(out, "/", y, ".breaks_norm.bed")
output_file3 <- paste0(out, "/", y, ".total_reads_noChrM.txt")
output_file4 <- paste0(out, "/", y, ".perc_reads_quant_sites.txt")
output_file5 <- paste0(out, "/", y, ".breaks_genome.txt")
output_file6 <- paste0(out, "/", y, ".reads_quant_sites.txt")

write.table(as.data.frame(alpha), file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(breaks_norm), file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(total_reads_noChrM), file = output_file3, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(perc_reads_quant_sites), file = output_file4, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(breaks_genome), file = output_file5, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(ranges_a), file = output_file6, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)