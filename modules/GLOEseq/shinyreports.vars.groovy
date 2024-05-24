//vars for task shinyReports from catalog ChIPseq, version 1
SHINYREPS_PROJECT=PROJECT	//project directory
SHINYREPS_LOG=LOGS			//where the logs lie
SHINYREPS_QC=QC				//where the QC lie
SHINYREPS_RES=RESULTS		//where the results lie
SHINYREPS_PREFIX="sfb_[A-Za-z]+_\\\\d+_\\\\d+_"  //standard sample prefix
SHINYREPS_TARGETS="targets.txt"
SHINYREPS_BOWTIE_LOG=LOGS + "/bowtie2_pe"	//where the Bowtie logs lie
SHINYREPS_BAMINDEX_LOG=LOGS + "/BAMindexer"	//where the Samtools/BamIndexer logs lie
SHINYREPS_MARKDUPS_LOG=LOGS + "/MarkDups"	//where the MarkDups logs lie
SHINYREPS_FASTQC=FASTQC_OUTDIR + "/raw"		//where the Fastqc logs lie
SHINYREPS_FASTQC_LOG=LOGS + "/FastQC"  //where the Fastqc logs lie
SHINYREPS_FASTQSCREEN_OUT=FASTQSCREEN_OUTDIR		        //where the FastqScreen logs lie
SHINYREPS_FASTQSCREEN_LOG=LOGS + "/FastqScreen"  //where the FastqScreen logs lie
SHINYREPS_MACS2=RESULTS + "/macs2"	//where the MACS2 results lie
SHINYREPS_MACS2_LOG=LOGS + "/macs2" //where the macs2 logs lie
SHINYREPS_PLOTS_COLUMN=2            //number of columns to splits the plot grids (ipstrength, phantompeaks...). Min=2
SHINYREPS_PEAK_ANNOTATION=RESULTS + "/Break_Annotation" // where the peak annotation results lie
SHINYREPS_TOOL_VERSIONS= PROJECT + "/GLOEPipe/modules/GLOEseq/tool.versions.groovy" //where the tool versions listed
SHINYREPS_BREAKS_DETECTED_RE1=RESULTS + "/macs2/EB1_breaks_detected_table.csv" 
SHINYREPS_BREAKS_DETECTED_RE2=RESULTS + "/macs2/EB2_breaks_detected_table.csv" 
SHINYREPS_BREAKS_DETECTED_RE3=RESULTS + "/macs2/EB3_breaks_detected_table.csv"
SHINYREPS_ALPHA=RESULTS + "/alpha"