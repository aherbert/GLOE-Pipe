//rule for task shinyReports from catalog GLOEseq, version 1
//desc: Creates the source code to compile the shiny and markdown reports
shinyReports = {
    doc title: "shinyReports",
        desc:  "creates the source code to compile the shiny and markdown reports",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Giuseppe Petrosino"

    output.dir = REPORTS
    
    produce("shinyReports.txt") {
        exec """
            cp ${MODULE_FOLDER}/../tools/reports/shiny_gloeseq_reporting_tool/server.R ${REPORTS}                &&
            cp ${MODULE_FOLDER}/../tools/reports/shiny_gloeseq_reporting_tool/ui.R ${REPORTS}                    &&
            cp ${MODULE_FOLDER}/../tools/reports/shiny_gloeseq_reporting_tool/GLOE.shinyrep.helpers.R ${REPORTS} &&
            cp ${MODULE_FOLDER}/../tools/reports/shiny_gloeseq_reporting_tool/styles.css ${REPORTS}               &&
            
            if [ -e "${REPORTS}/GLOEreport.Rmd" ]; then
                echo 'GLOEreport.Rmd already exists. Older copy will be kept and not overwritten';
            else
                cp ${MODULE_FOLDER}/../tools/reports/shiny_gloeseq_reporting_tool/GLOEreport.Rmd ${REPORTS};
            fi &&
    
            PROJECT=\$(basename ${SHINYREPS_PROJECT})                              &&
            sed -i "2,2s/SHINYREPS_PROJECT/\${PROJECT}/" ${REPORTS}/GLOEreport.Rmd &&
    
            echo "SHINYREPS_PROJECT=${SHINYREPS_PROJECT}" >  $output &&
            echo "SHINYREPS_LOG=${SHINYREPS_LOG}"         >> $output &&
            echo "SHINYREPS_QC=${SHINYREPS_QC}"           >> $output &&
            echo "SHINYREPS_RES=${SHINYREPS_RES}"         >> $output &&
            echo "SHINYREPS_TARGETS=${SHINYREPS_TARGETS}" >> $output &&
            echo "SHINYREPS_BOWTIE_LOG=${SHINYREPS_BOWTIE_LOG}"      >> $output &&
            echo "SHINYREPS_BAMINDEX_LOG=${SHINYREPS_BAMINDEX_LOG}"  >> $output &&
            echo "SHINYREPS_MARKDUPS_LOG=${SHINYREPS_MARKDUPS_LOG}"  >> $output &&
            echo "SHINYREPS_FASTQC=${SHINYREPS_FASTQC}"   >> $output &&
            echo "SHINYREPS_FASTQC_LOG=${SHINYREPS_FASTQC_LOG}"   >> $output &&
			echo "SHINYREPS_FASTQSCREEN_OUT=${SHINYREPS_FASTQSCREEN_OUT}"     >> $output &&
			echo "SHINYREPS_FASTQSCREEN_LOG=${SHINYREPS_FASTQSCREEN_LOG}"     >> $output &&
            echo "SHINYREPS_MACS2=${SHINYREPS_MACS2}"     >> $output &&
            echo "SHINYREPS_MACS2_LOG=${SHINYREPS_MACS2_LOG}"         >> $output &&
            echo "SHINYREPS_PREFIX=${SHINYREPS_PREFIX}"   >> $output &&
            echo "SHINYREPS_PLOTS_COLUMN=${SHINYREPS_PLOTS_COLUMN}" >> $output &&
            echo "SHINYREPS_PEAK_ANNOTATION=${SHINYREPS_PEAK_ANNOTATION}" >> $output &&
            echo "SHINYREPS_BREAKS_DETECTED_RE1=${SHINYREPS_BREAKS_DETECTED_RE1}" >> $output &&            
            echo "SHINYREPS_BREAKS_DETECTED_RE2=${SHINYREPS_BREAKS_DETECTED_RE2}" >> $output &&            
            echo "SHINYREPS_BREAKS_DETECTED_RE3=${SHINYREPS_BREAKS_DETECTED_RE3}" >> $output && 
            echo "SHINYREPS_ALPHA=${SHINYREPS_ALPHA}" >> $output && 
            echo "SHINYREPS_TOOL_VERSIONS=${SHINYREPS_TOOL_VERSIONS}" >> $output
        ""","shinyReports"
    }
}