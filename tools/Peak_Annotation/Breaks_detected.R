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

re1 <- parseArgs(args, "RestrictionEnzyme1=")
re2 <- parseArgs(args, "RestrictionEnzyme2=")
re3 <- parseArgs(args, "RestrictionEnzyme3=")

out <- parseArgs(args,"out=", "macs2")

runstr <- paste0("Call with: Rscript Breaks_detected.R [breakData=",breakData,"] [RestrictionEnzyme1=",re1,"] [RestrictionEnzyme2=",re2,"] [RestrictionEnzyme3=",re3,"] [out=",out,"]")
cat(runstr)

breakFiles <-list.files(breakData,pattern=".bed", full.names = TRUE) # list of the full path of the .xls file 
breaks <- lapply(breakFiles, readPeakFile) # read all the xls files using 'readPeakFile' function
breaks.re1 <- import(re1)

# check breaks overlapping with those expected
true.breaklst.re1 <- lapply(breaks, function(x) {
	m <- x[x %over% breaks.re1]
})

# create a summary table
filename <- strsplit(basename(breakFiles), "_macs2_summits.bed") # take the filenames 
filename.re1 <- strsplit(basename(re1), ".bed")

breaks.number <- sapply(breaks, function(x) { mm <-length(x) })

true.breaklst.re1.number <- sapply(true.breaklst.re1, function(x) { mm <-length(x) })
true.breaklst.re1.percentage <- (100/breaks.number)*true.breaklst.re1.number
re1.percentage <- (100/length(breaks.re1))*true.breaklst.re1.number

summary.table.re1 <- cbind(filename,breaks.number,true.breaklst.re1.number,true.breaklst.re1.percentage,
                           length(breaks.re1),true.breaklst.re1.number,re1.percentage)
colnames(summary.table.re1) <- c("Samples", "Breaks detected", paste(filename.re1,"Breaks detected", sep=" "),
                                 "Breaks detected (%)", paste(filename.re1,"Breaks", sep=" "), 
                                 paste(filename.re1,"Breaks detected", sep=" "),
                                 paste(filename.re1,"Breaks","(%)", sep=" ") )

names(true.breaklst.re1) <- paste0(filename, "_", filename.re1)

outputData.re1 <- lapply(true.breaklst.re1, as.data.frame)
sapply(names(outputData.re1), 
 function (x) write.table(outputData.re1[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

write.csv(summary.table.re1, file=paste0(out, "/", "RE1_breaks_detected_table.csv"), row.names=F)


if(file.exists(re2)) {

	breaks.re2 <- import(re2)

	# check breaks overlapping with those expected
	true.breaklst.re2 <- lapply(breaks, function(x) {
		m <- x[x %over% breaks.re2]
	})

	filename.re2 <- strsplit(basename(re2), ".bed")

	true.breaklst.re2.number <- sapply(true.breaklst.re2, function(x) { mm <-length(x) })
	true.breaklst.re2.percentage <- (100/breaks.number)*true.breaklst.re2.number
        re2.percentage <- (100/length(breaks.re2))*true.breaklst.re2.number

	summary.table.re2 <- cbind(filename,breaks.number,true.breaklst.re2.number,true.breaklst.re2.percentage,
                                   length(breaks.re2),true.breaklst.re2.number,re2.percentage)
        colnames(summary.table.re2) <- c("Samples", "Breaks detected", paste(filename.re2,"Breaks detected", sep=" "),
                                         "Breaks detected (%)", paste(filename.re2,"Breaks", sep=" "), 
                                         paste(filename.re2,"Breaks detected", sep=" "), 
                                         paste(filename.re2,"Breaks","(%)", sep=" ") )
  
	names(true.breaklst.re2) <- paste0(filename, "_", filename.re2)
	
	outputData.re2 <- lapply(true.breaklst.re2, as.data.frame)
	sapply(names(outputData.re2), 
	 function (x) write.table(outputData.re2[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

	write.csv(summary.table.re2, file=paste0(out, "/", "RE2_breaks_detected_table.csv"), row.names=F)

        save(re2,outputData.re2, file=paste0(out,"/Breaks_detected_re2.RData"))
}


if(file.exists(re3)) {

        breaks.re3 <- import(re3)

        # check breaks overlapping with those expected
        true.breaklst.re3 <- lapply(breaks, function(x) {
                m <- x[x %over% breaks.re3]
        })

        filename.re3 <- strsplit(basename(re3), ".bed")

        true.breaklst.re3.number <- sapply(true.breaklst.re3, function(x) { mm <-length(x) })
        true.breaklst.re3.percentage <- (100/breaks.number)*true.breaklst.re3.number
        re3.percentage <- (100/length(breaks.re3))*true.breaklst.re3.number

        summary.table.re3 <- cbind(filename,breaks.number,true.breaklst.re3.number,true.breaklst.re3.percentage,
                                   length(breaks.re3),true.breaklst.re3.number,re3.percentage)
        colnames(summary.table.re3) <- c("Samples", "Breaks detected", paste(filename.re3,"Breaks detected", sep=" "),
                                         "Breaks detected (%)", paste(filename.re3,"Breaks", sep=" "),
                                         paste(filename.re3,"Breaks detected", sep=" "),
                                         paste(filename.re3,"Breaks","(%)", sep=" ") )

        names(true.breaklst.re3) <- paste0(filename, "_", filename.re3)

        outputData.re3 <- lapply(true.breaklst.re3, as.data.frame)
        sapply(names(outputData.re3),
         function (x) write.table(outputData.re3[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

        write.csv(summary.table.re3, file=paste0(out, "/", "RE3_breaks_detected_table.csv"), row.names=F)
         
        save(re3,outputData.re3, file=paste0(out,"/Breaks_detected_re3.RData"))
}



# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/GLOEseq_Breaks_detected_session_info.txt", sep=""))
save(breakFiles,breaks,re1,filename,outputData.re1, file=paste0(out,"/Breaks_detected_re1.RData"))
