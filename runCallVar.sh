genome="/scratch/bucket/Rattus_norvegicus_nor6.fa"
dbsnp="/scratch/bucket/00-All.vcf.gz"
intervals="/scratch/bucket/Rnor_6.0.102.interval_list"
tmp="/data/tmp"
volume="/scratch"

# run HaplotypeCaller
for srafile in $(cat $1)
do
        sample=$(basename $srafile)

        docker run -v $volume:$volume -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller \
        -R $genome -I /data/output/$sample.sorted.dup.recal.bam \
        -L $intervals \
        --native-pair-hmm-threads 8 \
        -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        -O /data/output/$sample.vcf.gz -ERC GVCF

done
