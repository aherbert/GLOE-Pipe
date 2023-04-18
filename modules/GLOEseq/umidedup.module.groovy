umidedup = {
    doc title: "deduplication based on UMIs",
        desc:  "Deduplication of mapped data using UMIs with umi_tools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.5.slurm",
        author: "Giuseppe Petrosino"

    output.dir = MAPPED


	transform(".bam") to (".umidedup.bam") {
		exec """

		module load umitools/${UMITOOLS_VERSION} &&

                if [[ "${ESSENTIAL_PAIRED}" == "yes" ]]; then
			                    umi_tools dedup -I $input --paired -S $output --output-stats=${output.prefix}.stats;
                else
			                    umi_tools dedup -I $input -S $output --output-stats=${output.prefix}.stats;
                fi;

        ""","umidedup"
    }
}
