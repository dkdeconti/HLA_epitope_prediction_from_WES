for SAMPLENAME in $(cut -f 1 $1 | sort | uniq); do
    grep $SAMPLENAME $1 | cut -f 2 | sort | awk 'NR % 2 == 1' > first.txt;
    grep $SAMPLENAME $1 | cut -f 2 | sort | awk 'NR % 2 == 0' > second.txt;
    paste first.txt second.txt > ${SAMPLENAME}.tsv;
    echo $SAMPLENAME > tsv_list.txt;
    rm first.txt second.txt;
done;