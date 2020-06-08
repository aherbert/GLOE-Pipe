#############################################################################
##
## What: Breaks_detected.R
## Who: Giuseppe Petrosino
## When: 2018-03-14
##
## Script to select breaks detected overlapping with those expected
##
## Args: 
## -----
## breakData=                       
############################################################################
options(stringsAsFactors=FALSE)
library(ChIPseeker)
library(rtracklayer)

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
 
   if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}


args <- commandArgs(T)
breakData <- parseArgs(args,"breakData=") # .xls result from MACS2

eb1 <- parseArgs(args, "ExpectedBreaks1=")
eb2 <- parseArgs(args, "ExpectedBreaks2=")
eb3 <- parseArgs(args, "ExpectedBreaks3=")

out <- parseArgs(args,"out=", "macs2")

runstr <- paste0("Call with: Rscript Breaks_detected.R [breakData=",breakData,"] [ExpectedBreaks1=",eb1,"] [ExpectedBreaks1=",eb2,"] [ExpectedBreaks1=",eb3,"] [out=",out,"]")
cat(runstr)

breakFiles <-list.files(breakData,pattern=".bed", full.names = TRUE) # list of the full path of the .xls file 
breaks <- lapply(breakFiles, readPeakFile) # read all the xls files using 'readPeakFile' function
breaks.eb1 <- import(eb1)

# check breaks overlapping with those expected
true.breaklst.eb1 <- lapply(breaks, function(x) {
	m <- x[x %over% breaks.eb1]
})

# create a summary table
filename <- strsplit(basename(breakFiles), "_macs2_summits.bed") # take the filenames 
filename.eb1 <- strsplit(basename(eb1), ".bed")

breaks.number <- sapply(breaks, function(x) { mm <-length(x) })

true.breaklst.eb1.number <- sapply(true.breaklst.eb1, function(x) { mm <-length(x) })
true.breaklst.eb1.percentage <- (100/breaks.number)*true.breaklst.eb1.number
eb1.percentage <- (100/length(breaks.re1))*true.breaklst.eb1.number

summary.table.eb1 <- cbind(filename,breaks.number,true.breaklst.eb1.number,true.breaklst.eb1.percentage,
                           length(breaks.eb1),true.breaklst.eb1.number,eb1.percentage)
colnames(summary.table.eb1) <- c("Samples", "Breaks detected", paste(filename.eb1,"Breaks detected", sep=" "),
                                 "Breaks detected (%)", paste(filename.eb1,"Breaks", sep=" "), 
                                 paste(filename.eb1,"Breaks detected", sep=" "),
                                 paste(filename.eb1,"Breaks","(%)", sep=" ") )

names(true.breaklst.eb1) <- paste0(filename, "_", filename.eb1)

outputData.eb1 <- lapply(true.breaklst.eb1, as.data.frame)
sapply(names(outputData.eb1), 
 function (x) write.table(outputData.eb1[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

write.csv(summary.table.eb1, file=paste0(out, "/", "EB1_breaks_detected_table.csv"), row.names=F)


if(file.exists(eb2)) {

	breaks.eb2 <- import(eb2)

	# check breaks overlapping with those expected
	true.breaklst.eb2 <- lapply(breaks, function(x) {
		m <- x[x %over% breaks.eb2]
	})

	filename.eb2 <- strsplit(basename(eb2), ".bed")

	true.breaklst.eb2.number <- sapply(true.breaklst.eb2, function(x) { mm <-length(x) })
	true.breaklst.eb2.percentage <- (100/breaks.number)*true.breaklst.eb2.number
        eb2.percentage <- (100/length(breaks.eb2))*true.breaklst.eb2.number

	summary.table.eb2 <- cbind(filename,breaks.number,true.breaklst.eb2.number,true.breaklst.eb2.percentage,
                                   length(breaks.eb2),true.breaklst.eb2.number,eb2.percentage)
        colnames(summary.table.eb2) <- c("Samples", "Breaks detected", paste(filename.eb2,"Breaks detected", sep=" "),
                                         "Breaks detected (%)", paste(filename.eb2,"Breaks", sep=" "), 
                                         paste(filename.eb2,"Breaks detected", sep=" "), 
                                         paste(filename.eb2,"Breaks","(%)", sep=" ") )
  
	names(true.breaklst.eb2) <- paste0(filename, "_", filename.eb2)
	
	outputData.eb2 <- lapply(true.breaklst.eb2, as.data.frame)
	sapply(names(outputData.eb2), 
	 function (x) write.table(outputData.eb2[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

	write.csv(summary.table.eb2, file=paste0(out, "/", "EB2_breaks_detected_table.csv"), row.names=F)

        save(eb2,outputData.eb2, file=paste0(out,"/Breaks_detected_eb2.RData"))
}


if(file.exists(eb3)) {

        breaks.eb3 <- import(eb3)

        # check breaks overlapping with those expected
        true.breaklst.eb3 <- lapply(breaks, function(x) {
                m <- x[x %over% breaks.eb3]
        })

        filename.eb3 <- strsplit(basename(eb3), ".bed")

        true.breaklst.eb3.number <- sapply(true.breaklst.eb3, function(x) { mm <-length(x) })
        true.breaklst.eb3.percentage <- (100/breaks.number)*true.breaklst.eb3.number
        eb3.percentage <- (100/length(breaks.eb3))*true.breaklst.eb3.number

        summary.table.eb3 <- cbind(filename,breaks.number,true.breaklst.eb3.number,true.breaklst.eb3.percentage,
                                   length(breaks.eb3),true.breaklst.eb3.number,eb3.percentage)
        colnames(summary.table.eb3) <- c("Samples", "Breaks detected", paste(filename.eb3,"Breaks detected", sep=" "),
                                         "Breaks detected (%)", paste(filename.eb3,"Breaks", sep=" "),
                                         paste(filename.eb3,"Breaks detected", sep=" "),
                                         paste(filename.eb3,"Breaks","(%)", sep=" ") )

        names(true.breaklst.eb3) <- paste0(filename, "_", filename.eb3)

        outputData.eb3 <- lapply(true.breaklst.eb3, as.data.frame)
        sapply(names(outputData.eb3),
         function (x) write.table(outputData.eb3[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

        write.csv(summary.table.eb3, file=paste0(out, "/", "EB3_breaks_detected_table.csv"), row.names=F)
         
        save(eb3,outputData.eb3, file=paste0(out,"/Breaks_detected_eb3.RData"))
}



# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GLOEseq_Breaks_detected_session_info.txt", sep=""))
save(breakFiles,breaks,eb1,filename,outputData.eb1, file=paste0(out,"/Breaks_detected_eb1.RData"))
