macs2 = {
    doc title: "MACS2",
        desc:  "MACS2 wrapper",
        constraints: "Only performs treatment versus control breakcalling",
        bpipe_version: "tested with bpipe 0.9.9.5.slurm",
        author: "Giuseppe Petrosino"

    output.dir = RESULTS + "/macs2"
    def MACS2_FLAGS= MACS2_MFOLD  + " " + 
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
                TREATMENT=\$(       echo \$TARGET | cut -f1 -d" ") &&
                Tname=\$(   echo \$TARGET | cut -f2 -d" ") &&
                CONTROL=\$(    echo \$TARGET | cut -f3 -d" ") &&
                Cname=\$(echo \$TARGET | cut -f4 -d" ") &&
                CompName=\$(echo \$TARGET | cut -f5 -d" ") &&
                TREATMENTFOR=\${TREATMENT%%.bed}".for.bed" && 
                TREATMENTREV=\${TREATMENT%%.bed}".rev.bed" && 
                TREATMENTnameFOR=\${Tname}"_for" && 
                TREATMENTnameREV=\${Tname}"_rev" && 
                CONTROLFOR=\${CONTROL%%.bed}".for.bed" && 
                CONTROLREV=\${CONTROL%%.bed}".rev.bed" && 
                CONTROLnameFOR=\${Cname}"_for" && 
                CONTROLnameREV=\${Cname}"_rev";   
                
                if [ "\$BED" != "\$INPUT" ]; then
                    echo "\${Tname} vs \${Cname}" >> $output &&
                    macs2 callpeak -t $MACS2_INPUT/\$TREATMENTFOR -c $MACS2_INPUT/\$CONTROLFOR -n \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2 $MACS2_FLAGS &&
                    macs2 callpeak -t $MACS2_INPUT/\$TREATMENTREV -c $MACS2_INPUT/\$CONTROLREV -n \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2 $MACS2_FLAGS &&
                    if [ \$? -ne 0 ]; then rm $output; fi &&
                    awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" "+"}' \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2_summits.bed > \${CompName}_macs2.bed &&
                    awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" "-"}' \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2_summits.bed >> \${CompName}_macs2.bed &&
                    bedtools sort -i \${CompName}_macs2.bed > \${CompName}_macs2_summits.bed &&
                    awk '{if(NR>20)print}' \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2_peaks.xls > \${CompName}_macs2_peaks.xls &&
                    awk '{if(NR>21)print}' \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2_peaks.xls >> \${CompName}_macs2_peaks.xls &&
                    cat \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2_peaks.narrowPeak \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2_peaks.narrowPeak > \${CompName}_macs2_peaks.narrowPeak && 
                    rm \${TREATMENTnameFOR}.vs.\${CONTROLnameFOR}_macs2* &&
                    rm \${TREATMENTnameREV}.vs.\${CONTROLnameREV}_macs2* &&
                    rm \${CompName}_macs2.bed    &&
                    mv \${CompName}_macs2* $output.dir ;
                fi;
            done
        ""","macs2"
    }

    forward input
}
