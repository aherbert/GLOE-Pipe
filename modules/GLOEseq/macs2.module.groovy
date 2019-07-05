macs2 = {
    doc title: "MACS2",
        desc:  "MACS2 wrapper",
        constraints: "Only performs treatment versus control peakcalling",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Giuseppe Petrosino"

    output.dir = RESULTS + "/macs2"
    MACS2_FLAGS= MACS2_MFOLD  + " " + 
                 MACS2_GSIZE  + " " + 
                 MACS2_BWIDTH + " " + 
		 MACS2_FORMAT + " " +
                 MACS2_EXTRA
   
    transform(".bed") to("_macs2.done") {
        exec """
            module load macs2/${MACS2_VERSION} &&
            module load bedtools/${BEDTOOLS_VERSION} &&

            touch $output;
            if [ ! -e $MACS2_TARGETS ]; then
                echo "Targets file $MACS2_TARGETS doesn't exist" >> $output &&
                exit 0;
            fi;
            
            BED=\$(basename $input) &&
            grep \$BED $MACS2_TARGETS | while read -r TARGET; do
                IP=\$(       echo \$TARGET | cut -f1 -d" ") &&
                IPname=\$(   echo \$TARGET | cut -f2 -d" ") &&
                INPUT=\$(    echo \$TARGET | cut -f3 -d" ") &&
                INPUTname=\$(echo \$TARGET | cut -f4 -d" ") &&
                IPFOR=\${IP%%.bed}".for.bed" && 
                IPREV=\${IP%%.bed}".rev.bed" && 
                IPnameFOR=\${IPname}"_for" && 
                IPnameREV=\${IPname}"_rev" && 
                INPUTFOR=\${INPUT%%.bed}".for.bed" && 
                INPUTREV=\${INPUT%%.bed}".rev.bed" && 
                INPUTnameFOR=\${INPUTname}"_for" && 
                INPUTnameREV=\${INPUTname}"_rev";   
                
                if [ "\$BED" != "\$INPUT" ]; then
                    echo "\${IPname} vs \${INPUTname}" >> $output &&
                    macs2 callpeak -t $MACS2_INPUT/\$IPFOR -c $MACS2_INPUT/\$INPUTFOR -n \${IPnameFOR}.vs.\${INPUTnameFOR}_macs2 $MACS2_FLAGS &&
                    macs2 callpeak -t $MACS2_INPUT/\$IPREV -c $MACS2_INPUT/\$INPUTREV -n \${IPnameREV}.vs.\${INPUTnameREV}_macs2 $MACS2_FLAGS &&
                    if [ \$? -ne 0 ]; then rm $output; fi &&
                    awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" "+"}' \${IPnameFOR}.vs.\${INPUTnameFOR}_macs2_summits.bed > \${IPname}.vs.\${INPUTname}_macs2.bed &&
                    awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" "-"}' \${IPnameREV}.vs.\${INPUTnameREV}_macs2_summits.bed >> \${IPname}.vs.\${INPUTname}_macs2.bed &&
                    bedtools sort -i \${IPname}.vs.\${INPUTname}_macs2.bed > \${IPname}.vs.\${INPUTname}_macs2_summits.bed &&
                    awk '{if(NR>20)print}' \${IPnameFOR}.vs.\${INPUTnameFOR}_macs2_peaks.xls > \${IPname}.vs.\${INPUTname}_macs2_peaks.xls &&
                    awk '{if(NR>21)print}' \${IPnameREV}.vs.\${INPUTnameREV}_macs2_peaks.xls >> \${IPname}.vs.\${INPUTname}_macs2_peaks.xls &&
                    cat \${IPnameFOR}.vs.\${INPUTnameFOR}_macs2_peaks.narrowPeak \${IPnameREV}.vs.\${INPUTnameREV}_macs2_peaks.narrowPeak > \${IPname}.vs.\${INPUTname}_macs2_peaks.narrowPeak && 
                    rm \${IPnameFOR}.vs.\${INPUTnameFOR}_macs2* &&
                    rm \${IPnameREV}.vs.\${INPUTnameREV}_macs2* &&
                    rm \${IPname}.vs.\${INPUTname}_macs2.bed    &&
                    mv \${IPname}.vs.\${INPUTname}_macs2* $output.dir ;
                fi;
            done
        ""","macs2"
    }

    forward input
}