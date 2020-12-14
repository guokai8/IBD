for i in *R1_clean*
do
    kraken2 --db /extDataS/kai/microbiome/kraken2-microbial/ --threads 40 --paired $i ${i%%_R1_clean.fq.gz}_R2_clean.fq.gz --report ${i%%_R1_clean.fq.gz}.report.txt >${i%%_R1_clean.fq.gz}.kraken.txt --report-zero-counts --use-names --classified-out ${i%%_R1_clean.fq.gz}#.seq
done

