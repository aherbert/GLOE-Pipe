breaks_detected = {

    doc title: "Breaks_detected",
        desc: "Select breaks detected overlapping with those expected",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino"

    output.dir=breaks_detected_OUTDIR.replaceFirst("out=","")
    def breaks_detected_FLAGS = breaks_detected_FILES + " " +
        breaks_detected_EB1 + " " + 
        breaks_detected_EB2 + " " + 
        breaks_detected_EB3 + " " +
        breaks_detected_OUTDIR + " " +
        breaks_detected_EXTRA

    produce("Breaks_detected_eb1.RData") {
        exec """
           module load bedtools/${BEDTOOLS_VERSION} &&
           module load R/${R_VERSION} &&
            Rscript ${TOOL_BREAKS_DETECTED}/Breaks_detected.R $breaks_detected_FLAGS;
        ""","breaks_detected"
    }
	forward input
}