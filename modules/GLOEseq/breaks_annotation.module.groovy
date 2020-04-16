breaks_annotation = {

    doc title: "breaks_annotation",
        desc: "",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino"

    output.dir=breaks_annotation_OUTDIR.replaceFirst("out=","")
    def breaks_annotation_FLAGS = breaks_annotation_FILES + " " +
        breaks_annotation_TRANSCRIPT_TYPE + " " + 
        breaks_annotation_TRANSCRIPT_DB + " " + 
        breaks_annotation_ORGDB + " " + 
        breaks_annotation_REGIONTSS + " " + 
        breaks_annotation_OUTDIR + " " +
        breaks_annotation_INDEX + " " +
        breaks_annotation_EXTRA


    produce("Breaks_Annotation.RData") {
        exec """

                     if [ -e ${TMP} ]; then
                           rm -r ${TMP};
                     fi &&

            module load R/${R_VERSION} &&

            Rscript ${TOOL_BREAKS_ANNOTATION}/Breaks_Annotation.R $breaks_annotation_FLAGS;
        ""","breaks_annotation"
    }
}
