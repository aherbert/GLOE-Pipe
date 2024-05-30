##################################
##
## helper functions to create the plots for the Shiny report
##
##################################
library("knitr")        # for markdown output
library("ChIPpeakAnno")	#for peak venn diagrams
library("RColorBrewer")
library("ggplot2")
library("gridExtra")

##
## loadGlobalVars: read configuration from bpipe vars
##
loadGlobalVars <- function(f="shinyReports.txt") {

    # read in the conf file
    conf <- readLines(f)
    conf <- conf[grep("^SHINYREPS_", conf)]
    
    # create the vars
    sapply(conf, function(x) {
        x <- unlist(strsplit(x, "=", fixed=T))
        assign(x[1], x[2], envir=.GlobalEnv)
    })
    
    invisible(0)
}

##
## Some generic functions
##
# shorten: if a text string is longer than certain length, shorten by showing the first and last characters
shorten <- function(x, max.len=40, ini=20, end=15) {
    l <- nchar(x)
    if(l > max.len) paste(substr(x, 1, ini), substr(x, (l-end), l), sep="...") else x
}

##
## GLOEhelper.init: some time consuming tasks that can be done in advance
##
GLOEhelper.init <- function(task) {
    
    # read targets.txt
    readTargets <- function() {
        TARGETS <- paste0(SHINYREPS_PROJECT, "/", SHINYREPS_TARGETS)
        if(!file.exists(TARGETS)) {
            return("Targets file not available")
        }
    
        return(read.delim(TARGETS))
    }
    
    # read peaks from MACS2 output
    readPeaks <- function() {
    
        # check which macs output exist
        if(!file.exists(SHINYREPS_MACS2)) {
            return("MACS2 results not available")
        }
        comparisons <- paste0(targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls")
        exist <- sapply(paste0(SHINYREPS_MACS2, "/", comparisons), file.exists)
        targets <- targets[exist, ]

        # and return the tables
        peaks <- lapply(paste0(SHINYREPS_MACS2, "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls"), function(x) {
            x <- tryCatch(read.delim(x, comment.char="#"), error=function(e) as.data.frame(matrix(ncol=10)))
            colnames(x) <- c("chr", "start", "end", "length", "summit", "tags", "-log10 pval", "fold enrichment", "-log10 FDR", "name")
            x[order(x$chr, x$start, x$end), c(-7, -10)]
        })
        names(peaks) <- paste0(targets$IPname, " vs. ", targets$INPUTname)

        return(peaks)
    }
    
    # dispatch tasks
    switch(task, 
           readTargets=readTargets(), 
           readPeaks=readPeaks())
}

##
## GLOEhelper.ComparisonsFromTargets: get te comparisons performed by MACS2 from the targets file
##
GLOEhelper.ComparisonsFromTargets <- function() {
    
    # check for targets.txt and macs2 results
    TARGETS <- paste0(SHINYREPS_PROJECT, "/", SHINYREPS_TARGETS)
    if(!file.exists(TARGETS)) {
        return("Targets file not available")
    }

    if(!file.exists(SHINYREPS_MACS2)) {
        return("MACS2 results not available")
    }
    
    # get the comparisons and clean the names
    x <- read.delim(TARGETS)
    comparisons <- paste0(x$IPname, ".vs.", x$INPUTname, "_macs2_peaks.xls")
    exist <- sapply(paste0(SHINYREPS_MACS2, "/", comparisons), file.exists)
    comparisons <- gsub(".vs.", " vs. ", comparisons)
    comparisons <- gsub("_macs2_peaks.xls", "", comparisons)
    
    return(comparisons[exist])
}

##
## GLOEhelper.Peaks: show the peaks called by MACS2
##
GLOEhelper.Peaks <- function(i=1) {
    ord  <- order(peaks[[i]][, "-log10 FDR"], 
                  peaks[[i]][, "fold enrichment"], 
                  decreasing=TRUE)
    peaks[[i]][ord, ]
}

##
## GLOEhelper.VennDiagram: shows a venn diagram per group of peaks called
##
GLOEhelper.VennDiagram <- function(){
	groups <- unique(targets$group)
	#create granges from the peaks
	peak.ranges <- lapply(peaks, function(x){
		  x <- GRanges(seqnames=x$chr,
			  IRanges(x$start,
				  end=x$end),
			  strand="*"
			   )
		  
		  })
	peak.groups <- targets$group
	for(group in groups){
		cat(paste0("#### ", group), fill=T)
		cat("\n", fill=T)
		peak <- peak.ranges[peak.groups==group]
		peaks.ov <- findOverlapsOfPeaks(peak)
		makeVennDiagram(peaks.ov,
				margin=0.1,
				cat.fontface=rep("bold", length(peak)),
				fill=brewer.pal(length(peak), "Accent")[1:length(peak)]
				)
		cat("\n", fill=T)
	}

}

##
## GLOEhelper.BOWTIE: parse bowtie log files and create a md table
##
GLOEhelper.Bowtie <- function() {
    
    # log file
	LOG <- SHINYREPS_BOWTIE_LOG   
    if(!file.exists(LOG)) {
        return("Bowtie statistics not available")
    }
    
    # look for the lines containing the strings
    # and get the values associated with this strings
    x <- sapply(list.files(LOG), function(f) {
        
        x <- file(paste0(LOG, "/", f))
        l <- readLines(x)
        close(x)

        names = unlist(
                  lapply(1:length(l), 
                  function (x) sub("[^[:alpha:]]+", "", l[x] ))
                 )


        stats = lapply(1:length(l), 
                     function(x){ 
                       a = sub("^\\s+", "", l[x]); 
                       b = sub("\\).*","\\)", a) ; 
                       c = unlist(strsplit(b, " "))[1:2]; 
                       d = c(c[1] , sub("[^[:digit:]]+", "", c[2])); 
                       e = gsub(")", "", d)
                      }
                    )

   
        # and add the duplicates information
        f <- gsub(".bam.log", ".duprm.bam.log", f)
        dups <- if(file.exists(paste0(SHINYREPS_MARKDUPS_LOG, "/", f))) {
            x <- file(paste0(SHINYREPS_MARKDUPS_LOG, "/", f))
            l <- readLines(x)
            close(x)
            gsub(".+Marking (\\d+) records as duplicates.+", "\\1", l[grep("Marking \\d+ records as duplicates", l)])
        } else {
            "not available"
        }
        

        names(stats) = names
        tmp1 = stack(stats)
        tmp2 = subset(tmp1, ind != "")
        tmp2$indm = paste(tmp2$ind, rep(c(1:2), length(tmp2$ind)/2, sep = ""))
        tmp3 = subset(tmp2,values != "")

        stats.return = as.character(tmp3$values)
        names(stats.return) = tmp3$indm

        
        c(stats.return, dups)
    })

    # set row and column names, and output the md table
    colnames(x) <- gsub(paste0("^", SHINYREPS_PREFIX), "", colnames(x))
    colnames(x) <- gsub(".bam.log$", "", colnames(x))
    df <- data.frame(sample_names=sapply(colnames(x), shorten), 
                     input_reads=format( as.numeric(x[1, ]), big.mark=","), 
                     times0=paste( format( as.numeric(x[4, ]), big.mark=","), " (", x[5, ],")", sep=""), 
                     time1=paste( format( as.numeric(x[6, ]), big.mark=","), " (", x[7, ], ")", sep=""), 
                     timesg1=paste0(format(as.numeric(x[8, ]), big.mark=","), " (", x[9, ], ")", sep="")
                     #duplicates=paste0(format(as.numeric(x[8, ]), big.mark=", "), " (", round(100 * as.numeric(x[8, ]) / as.numeric(x[2, ]), 2), "%)")
                     )
    kable(df, align=c("l", "r", "r", "r", "r", "r"), output=F, format="markdown", row.names=FALSE,
          col.names=c("sample names", "all reads", "unmapped (%)", "uniquely mapped (%)", "multiple mapped (%)"))
}

##
## GLOEhelper.Fastqc: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
GLOEhelper.Fastqc <- function(web=TRUE) {
    
    # logs folder
    if(!file.exists(SHINYREPS_FASTQC)) {
        return("Fastqc statistics not available")
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/fastqc" else SHINYREPS_FASTQC
    
    # construct the image url from the folder ents (skip current dir .)
    samples <- list.dirs(SHINYREPS_FASTQC, recursive=F)
    df <- sapply(samples, function(f) {
        c(paste0("![fastq dup img](", QC, "/", basename(f), "/Images/duplication_levels.png)"), 
          paste0("![fastq qual img](", QC, "/", basename(f), "/Images/per_base_quality.png)"),
          paste0("![fastq adap img](", QC, "/", basename(f), "/Images/adapter_content.png)"), 
          paste0("![fastq sequ img](", QC, "/", basename(f), "/Images/per_base_sequence_content.png)"))
    })

    # set row and column names, and output the md table
    df <- as.data.frame(t(df))
    x <- gsub(paste0("^", SHINYREPS_PREFIX), "", basename(samples))
    x <- gsub("_fastqc$", "", x)
    rownames(df) <- sapply(x, shorten)
    colnames(df) <- c("Duplication levels", "Read qualities", "Adapter content", "Sequence bias")
    kable(df, output=F, align="c", format="markdown")
}

##
## GLOEhelper.IPstrength: go through IPstrength output dir and create a md table with
##     the plots
##
GLOEhelper.IPstrength<- function(web=TRUE) {
    
    # logs folder
    if(!file.exists(SHINYREPS_IPSTRENGTH)) {
        return("IPstrength statistics not available")
    }
    
    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 4L    # default to 4 columns
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/ipstrength" else SHINYREPS_IPSTRENGTH
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_IPSTRENGTH, pattern=".png$")
    COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
    df <- sapply(samples, function(f) {
        paste0("![IPstrength img](", QC, "/", basename(f), ")")
    })
    
    # put sample names and output an md table of COLUMNS columns
    while(length(df) %% COLUMNS != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        sapply(gsub("_ipstrength.png)$", "", x), shorten)
    })
    df      <- matrix(df     , ncol=COLUMNS, byrow=T)
    samples <- matrix(samples, ncol=COLUMNS, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=COLUMNS, byrow=T)
    colnames(df.names) <- rep(" ", COLUMNS)
    
    kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## GLOEhelper.breakperChromosome: go through Peak_Annotation output dir and create a md table with
##      the Breaks per Chromosome plots
##
GLOEhelper.breakperChromosome <- function(web=TRUE) {
  # check if peak annotation results are available
  if(!file.exists(SHINYREPS_PEAK_ANNOTATION)){
    return("Breaks per Chromosome results not available")
  }

  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 2L    # default to 3 columns
  }

  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(SHINYREPS_PEAK_ANNOTATION, pattern="_GLOEseq_Breaks_per_Chromosomeplot.png$")
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![Peak_Annotation img](", SHINYREPS_PEAK_ANNOTATION, "/", basename(f), ")")
  })

  # put sample names and output an md table of COLUMN columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_GLOEseq_Breaks_per_Chromosomeplot.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)

  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }),
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)

  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}



##
## GLOEhelper.peakAnnotation: go through Peak_Annotation output dir and create a md table with 
##      the coverage plots
##
GLOEhelper.peakAnnotationCoverage <- function(web=TRUE) {
  # check if peak annotation results are available
  if(!file.exists(SHINYREPS_PEAK_ANNOTATION)){
    return("Peak annotation results not available")  
  }
  
  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 2L    # default to 3 columns
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(SHINYREPS_PEAK_ANNOTATION, pattern="_GLOEseq_Breaks_Coverageplot.png$")
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![Peak_Annotation img](", SHINYREPS_PEAK_ANNOTATION, "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMN columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_GLOEseq_Breaks_Coverageplot.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## GLOEhelper.peakAnnotationUpSet: go through Peak_Annotation output dir and create a md table with 
##      the UpSet plots
##
GLOEhelper.peakAnnotationUpSet <- function(web=TRUE) {
  # check if peak annotation results are available
  if(!file.exists(SHINYREPS_PEAK_ANNOTATION)){
    return("Peak annotation results not available")  
  }
  
  if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
    SHINYREPS_PLOTS_COLUMN <- 2L    # default to 3 columns
  }
  
  # construct the image url from the folder contents (skip current dir .)
  samples <- list.files(SHINYREPS_PEAK_ANNOTATION, pattern="_GLOEseq_UpSetplot.png$")
  COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
  df <- sapply(samples, function(f) {
    paste0("![Peak_Annotation img](", SHINYREPS_PEAK_ANNOTATION, "/", basename(f), ")")
  })
  
  # put sample names and output an md table of COLUMN columns
  while(length(df) %% COLUMNS != 0) df <- c(df, "")
  samples <- sapply(df, function(x) {
    x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
    sapply(gsub("_GLOEseq_UpSetplot.png)$", "", x), shorten)
  })
  df      <- matrix(df     , ncol=COLUMNS, byrow=T)
  samples <- matrix(samples, ncol=COLUMNS, byrow=T)
  
  # add a row with the sample names
  df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                     ncol=COLUMNS, byrow=T)
  colnames(df.names) <- rep(" ", COLUMNS)
  
  kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## GLOEhelper.PhantomPeak: go through PhantomPeak output dir and create a md table with
##     the plots
##
GLOEhelper.PhantomPeak <- function(web=TRUE) {
    
    # logs folder
    if(!file.exists(SHINYREPS_PHANTOMPEAK)) {
        return("PhantomPeak statistics not available")
    }
    
    if(!is.integer(SHINYREPS_PLOTS_COLUMN) | SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 2L    # default to 4 columns
    }
    
    # construct the folder name, which is different for web and noweb
    QC <- if(web) "/phantompeak" else SHINYREPS_PHANTOMPEAK
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_PHANTOMPEAK, pattern=".png$")
    COLUMNS <- min(length(samples), SHINYREPS_PLOTS_COLUMN)
    df <- sapply(samples, function(f) {
        paste0("![PhantomPeak img](", QC, "/", basename(f), ")")
    })
    
    # put sample names and output an md table of COLUMN columns
    while(length(df) %% COLUMNS != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        sapply(gsub("_phantompeak.png)$", "", x), shorten)
    })
    df      <- matrix(df     , ncol=COLUMNS, byrow=T)
    samples <- matrix(samples, ncol=COLUMNS, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=COLUMNS, byrow=T)
    colnames(df.names) <- rep(" ", COLUMNS)
    
    kable(as.data.frame(df.names), output=F, align="c", format="markdown")
}

##
## GLOEhelper.PBC: go through PBC output dir and create a md table with
##     the PBC stats
##
GLOEhelper.PBC <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_PBC)) {
        return("PCR bottleneck coefficient statistics not available")
    }
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_PBC, pattern="*.csv")
    df <- sapply(samples, function(f) {
        read.csv(paste0(SHINYREPS_PBC, "/", f))$PBC
    })
    
    # output md table
    df <- as.data.frame(df)
    colnames(df) <- "PBC"
    rownames(df) <- gsub("_PBC.csv", "", rownames(df))
    kable(as.data.frame(df), output=F, format="markdown")
}

##
## GLOEhelper.Bustard: call the perl XML interpreter and get the MD output
##
GLOEhelper.Bustard <- function() {
    f  <- SHINYREPS_BUSTARD
    
    if(!file.exists(f)) {
        return("Bustard statistics not available")
    }
    
    # call the perl XSL inetrpreter
    cmd <- paste(" bustard.pl", f)
    try(ret <- system2("perl", cmd, stdout=TRUE, stderr=FALSE))
    
    # check RC
    if(!is.null(attributes(ret))) {
        return(paste("Error parsing bustard statistics. RC:", attributes(ret)$status, "in command: perl", cmd))
    }
    
    ret     # ret contains already MD code
}

##
## extract tool versions
##
## report version of used tools
Toolhelper.ToolVersions <- function() {
    ver <- read.table(file=SHINYREPS_TOOL_VERSIONS,sep="=")
        ver$V1 <- strsplit(as.character(ver$V1),"_VERSION")
        colnames(ver) <- c("Tool name","Version")

            kable(as.data.frame(ver),output=F)
}

##
## GLOEhelper.BREAKSDETECTED: go through MACS2 output dir and show Breaks detected 
##                            for three restriction enzymes
##     Filter module file
##
GLOEhelper.BREAKSDETECTED.RE1 <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_BREAKS_DETECTED_RE1)) {
        return("Breaks overlapping Restriction enzyme breaks not available")
    }
    
    # construct the image url from the folder contents (skip current dir .)
    samples_RE1 <- read.csv(SHINYREPS_BREAKS_DETECTED_RE1, row.names=1, check.names=FALSE)
    kable(samples_RE1, row.names=T, align="r")
}

GLOEhelper.BREAKSDETECTED.RE2 <- function() {

    # logs folder
    if(!file.exists(SHINYREPS_BREAKS_DETECTED_RE2)) {
        return("Breaks overlapping Restriction enzyme breaks not available")
    }
    
    # construct the image url from the folder contents (skip current dir .)
    samples_RE2 <- read.csv(SHINYREPS_BREAKS_DETECTED_RE2, row.names=1, check.names=FALSE)
    kable(samples_RE2, row.names=T, align="r")
}

GLOEhelper.BREAKSDETECTED.RE3 <- function() {
 
    # logs folder
    if(!file.exists(SHINYREPS_BREAKS_DETECTED_RE3)) {
        return("Breaks overlapping Restriction enzyme breaks not available")     
    }   
         
    # construct the image url from the folder contents (skip current dir .)
    samples_RE3 <- read.csv(SHINYREPS_BREAKS_DETECTED_RE3, row.names=1, check.names=FALSE)
    kable(samples_RE3, row.names=T, align="r")

}


##
## smallRNAhelper.fastqscreen: go through Fastqc output dir and create a md table with the duplication & read quals & sequence bias plots
##
GLOEhelper.fastqscreen <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_FASTQSCREEN_OUT)) {
        return("FastQScreen statistics not available")
    }

    # construct the folder name, which is different for web and noweb
    QC <- SHINYREPS_FASTQSCREEN_OUT

    SHINYREPS_PLOTS_COLUMN <- tryCatch(as.integer(SHINYREPS_PLOTS_COLUMN),error=function(e){4})
    if(SHINYREPS_PLOTS_COLUMN < 2) {
        SHINYREPS_PLOTS_COLUMN <- 1L    # default to 4 columns
    }
    
    # construct the image url from the folder contents (skip current dir .)
    samples <- list.files(SHINYREPS_FASTQSCREEN_OUT, pattern=".png$", recursive=T, full.names=T)
    df <- sapply(samples, function(f) {
        paste0("![fastqscreen img](", f, ")")
    })
    
    # put sample names and output an md table of SHINYREPS_PLOTS_COLUMN columns
    while(length(df) %% SHINYREPS_PLOTS_COLUMN != 0) df <- c(df, "")
    samples <- sapply(df, function(x) {
        x <- sapply(x, function(x) gsub(paste0("^", SHINYREPS_PREFIX), "", basename(x)))
        gsub(".cutadapt|.highQ|.deduped|.trimmed|_screen.png)","", x)
    })

    df      <- matrix(df     , ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    samples <- matrix(samples, ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    
    # add a row with the sample names
    df.names <- matrix(sapply(1:nrow(df), function(i) { c(df[i, ], samples[i, ]) }), 
                       ncol=SHINYREPS_PLOTS_COLUMN, byrow=T)
    colnames(df.names) <- rep(" ", SHINYREPS_PLOTS_COLUMN)
    
    kable(as.data.frame(df.names), align="c", output=F, format="markdown")
}




GLOEhelper.ALPHA.3end <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_ALPHA)) {
        return("Alpha values are not available")
    }
    
    # Alpha values
    samples_3ends_alpha <- list.files(paste0(SHINYREPS_ALPHA, "/3end"), pattern=".alpha.txt$", recursive=T, full.names=T)
    alpha <- lapply(samples_3ends_alpha, read.delim, head=F)
    alpha <- do.call(rbind, alpha)

    # Total reads (without chrM)
    samples_3ends_total_reads_noChrM <- list.files(paste0(SHINYREPS_ALPHA, "/3end"), pattern=".total_reads_noChrM.txt$", recursive=T, full.names=T)
    total_reads_noChrM <- lapply(samples_3ends_total_reads_noChrM, read.delim, head=F)
    total_reads_noChrM <- do.call(rbind, total_reads_noChrM)

    # Perc nuclear reads
    samples_3ends_perc_nuclear_reads <- list.files(paste0(SHINYREPS_ALPHA, "/3end"), pattern=".perc_nuclear_reads.txt$", recursive=T, full.names=T)
    perc_nuclear_reads <- lapply(samples_3ends_perc_nuclear_reads, read.delim, head=F)
    perc_nuclear_reads <- do.call(rbind, perc_nuclear_reads)

    # Breaks genome
    samples_3ends_breaks_genome <- list.files(paste0(SHINYREPS_ALPHA, "/3end"), pattern=".breaks_genome.txt$", recursive=T, full.names=T)
    breaks_genome <- lapply(samples_3ends_breaks_genome, read.delim, head=F)
    breaks_genome <- do.call(rbind, breaks_genome)

    # Perc reads quant sites
    samples_3ends_perc_reads_quant_sites <- list.files(paste0(SHINYREPS_ALPHA, "/3end"), pattern=".perc_reads_quant_sites.txt$", recursive=T, full.names=T)
    perc_reads_quant_sites <- lapply(samples_3ends_perc_reads_quant_sites, read.delim, head=F)
    perc_reads_quant_sites <- do.call(rbind, perc_reads_quant_sites)


    samples <- gsub(".R1.m.end3.trimmed.sr.alpha.txt", "", basename(samples_3ends_alpha))
    samples <- gsub(SHINYREPS_PREFIX, "", samples)

    df <- cbind(samples, "3end", total_reads_noChrM, perc_nuclear_reads, perc_reads_quant_sites, alpha, breaks_genome)
    colnames(df) <- c("sample_name", "end_type", "total_count", "perc_nuclear_reads", "perc_quant_site", "alpha", "breaks_per_genome")

    pal <- brewer.pal(4, "BuGn")

    total_count.p <- ggplot(df, aes(sample_name,total_count/1e7)) +
                    geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[1]) +
                    labs(x="", y="Total # reads (1e7)") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))

    nuclear_reads.p <- ggplot(df, aes(sample_name,perc_nuclear_reads)) +
                    geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[2]) +
                    labs(x="", y="% nuclear reads") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))

    quant_site.p <- ggplot(df, aes(sample_name,perc_quant_site)) +
                   geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[3]) +
                   labs(x="", y="% reads at quantification sites") +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))

    breaks_genome.p <- ggplot(df, aes(sample_name,breaks_per_genome, fill = end_type)) +
                      geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[4]) +
                      labs(x="", y="Breaks per genome") +
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))
       
    kable(as.data.frame(df), align="c", output=F, format="markdown")


    samples_3ends_reads_quant <- list.files(paste0(SHINYREPS_ALPHA, "/3end"), pattern="sr.reads_quant_sites.txt$", recursive=T, full.names=T)
    reads_quant <- lapply(samples_3ends_reads_quant, read.delim, head=F)


    samples <- gsub(".R1.m.end3.trimmed.sr.reads_quant_sites.txt", "", basename(samples_3ends_reads_quant))
    samples <- gsub(SHINYREPS_PREFIX, "", samples)
    
    for (x in 1:length(reads_quant)) {
    reads_quant[[x]]$samples <- samples[x]
    }

    reads_quant <- do.call(rbind, reads_quant)
      


    reads_quant$sites <- paste0(reads_quant$V1,",",reads_quant$V2)

    pal <- rep(brewer.pal(10, "Spectral"), 8)

    reads_quant.plot <- ggplot(reads_quant) +
                               geom_point(aes(x = sites, y = V7, color = sites), show.legend = FALSE) +
                               scale_color_manual(values=pal) +
                               theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                                     axis.title = element_text(size = 7)) +
                               labs(title = "3' end", y = "# reads per quant site / Total # reads at quant sites per samples * 100")

    return(list(table.3end = kable(as.data.frame(df), align="c", output=F, format="markdown" , row.names=FALSE), total_count.p = total_count.p, nuclear_reads.p = nuclear_reads.p, quant_site.p = quant_site.p, breaks_genome.p = breaks_genome.p, reads_quant.plot = reads_quant.plot ))
     
}



GLOEhelper.ALPHA.5end <- function() {
    
    # logs folder
    if(!file.exists(SHINYREPS_ALPHA)) {
        return("Alpha values are not available")
    }
    
    # Alpha values
    samples_5ends_alpha <- list.files(paste0(SHINYREPS_ALPHA, "/5end"), pattern=".alpha.txt$", recursive=T, full.names=T)
    alpha <- lapply(samples_5ends_alpha, read.delim, head=F)
    alpha <- do.call(rbind, alpha)

    # Total reads (without chrM)
    samples_5ends_total_reads_noChrM <- list.files(paste0(SHINYREPS_ALPHA, "/5end"), pattern=".total_reads_noChrM.txt$", recursive=T, full.names=T)
    total_reads_noChrM <- lapply(samples_5ends_total_reads_noChrM, read.delim, head=F)
    total_reads_noChrM <- do.call(rbind, total_reads_noChrM)


    # Perc nuclear reads
    samples_5ends_perc_nuclear_reads <- list.files(paste0(SHINYREPS_ALPHA, "/5end"), pattern=".perc_nuclear_reads.txt$", recursive=T, full.names=T)
    perc_nuclear_reads <- lapply(samples_5ends_perc_nuclear_reads, read.delim, head=F)
    perc_nuclear_reads <- do.call(rbind, perc_nuclear_reads)


    # Breaks genome
    samples_5ends_breaks_genome <- list.files(paste0(SHINYREPS_ALPHA, "/5end"), pattern=".breaks_genome.txt$", recursive=T, full.names=T)
    breaks_genome <- lapply(samples_5ends_breaks_genome, read.delim, head=F)
    breaks_genome <- do.call(rbind, breaks_genome)

    # Perc reads quant sites
    samples_5ends_perc_reads_quant_sites <- list.files(paste0(SHINYREPS_ALPHA, "/5end"), pattern=".perc_reads_quant_sites.txt$", recursive=T, full.names=T)
    perc_reads_quant_sites <- lapply(samples_5ends_perc_reads_quant_sites, read.delim, head=F)
    perc_reads_quant_sites <- do.call(rbind, perc_reads_quant_sites)


    samples <- gsub(".R1.m.end5.trimmed.sr.alpha.txt", "", basename(samples_5ends_alpha))
    samples <- gsub(SHINYREPS_PREFIX, "", samples)

    df <- cbind(samples, "5end", total_reads_noChrM, perc_nuclear_reads, perc_reads_quant_sites, alpha, breaks_genome)
    colnames(df) <- c("sample_name", "end_type", "total_count", "perc_nuclear_reads", "perc_quant_site", "alpha", "breaks_per_genome")

    pal <- brewer.pal(4, "BuGn")

    total_count.p <- ggplot(df, aes(sample_name,total_count/1e7)) +
                    geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[1]) +
                    labs(x="", y="Total # reads (1e7)") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))

    nuclear_reads.p <- ggplot(df, aes(sample_name,perc_nuclear_reads)) +
                    geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[2]) +
                    labs(x="", y="% nuclear reads") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))


    quant_site.p <- ggplot(df, aes(sample_name,perc_quant_site)) +
                   geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[3]) +
                   labs(x="", y="% reads at quantification sites") +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))

    breaks_genome.p <- ggplot(df, aes(sample_name,breaks_per_genome)) +
                      geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, fill = pal[4]) +
                      labs(x="", y="Breaks per genome") +
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                         axis.title = element_text(size = 7))
       
    kable(as.data.frame(df), align="c", output=F, format="markdown")


    samples_5ends_reads_quant <- list.files(paste0(SHINYREPS_ALPHA, "/5end"), pattern="sr.reads_quant_sites.txt$", recursive=T, full.names=T)
    reads_quant <- lapply(samples_5ends_reads_quant, read.delim, head=F)


    samples <- gsub(".R1.m.end5.trimmed.sr.reads_quant_sites.txt", "", basename(samples_5ends_reads_quant))
    samples <- gsub(SHINYREPS_PREFIX, "", samples)
    
    for (x in 1:length(reads_quant)) {
    reads_quant[[x]]$samples <- samples[x]
    }

    reads_quant <- do.call(rbind, reads_quant)
      


    reads_quant$sites <- paste0(reads_quant$V1,",",reads_quant$V2)

    pal <- rep(brewer.pal(10, "Spectral"), 8)

    reads_quant.plot <- ggplot(reads_quant) +
                               geom_point(aes(x = sites, y = V7, color = sites), show.legend = FALSE) +
                               scale_color_manual(values=pal) +
                               theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
                                     axis.title = element_text(size = 7)) +
                               labs(title = "5' end", y = "# reads per quant site / Total # reads at quant sites per samples * 100")

    return(list(table.5end = kable(as.data.frame(df), align="c", output=F, format="markdown" , row.names=FALSE), total_count.p = total_count.p, nuclear_reads.p = nuclear_reads.p, quant_site.p = quant_site.p, breaks_genome.p = breaks_genome.p, reads_quant.plot = reads_quant.plot ))
     
}