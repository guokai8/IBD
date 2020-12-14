for i in *R1.fq.gz; 
do 
     n=${i%%_R1.fq.gz}; 
     java -jar /data/biotools/Trimmomatic/trimmomatic-0.36.jar PE -threads 40 -phred33 $i ${n}_R2.fq.gz ${n}_R1_clean.fq.gz ${n}_R1_unpair.fq.gz ${n}_R2_clean.fq.gz ${n}_R2_unpair.fq.gz\
         ILLUMINACLIP:/data/biotools/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 HEADCROP:10;
done 
