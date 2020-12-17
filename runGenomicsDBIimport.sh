intervals="/scratch/bucket/Rnor_6.0.102.interval_list"
volume="/scratch"

# run GenomicsDBImport
docker run  -v $volume:$volume -v $(pwd):/data/ broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=5" GenomicsDBImport \
        --genomicsdb-workspace-path /data/genomicsdb \
        -V /data/output/SRR1594108.2.vcf.gz \
        -V /data/output/SRR1594109.2.vcf.gz \
        -V /data/output/SRR1594110.2.vcf.gz \
        -V /data/output/SRR1594111.2.vcf.gz \
        -V /data/output/SRR1594112.2.vcf.gz \
        -V /data/output/SRR1594113.2.vcf.gz \
        -V /data/output/SRR1594114.2.vcf.gz \
        -V /data/output/SRR1594115.2.vcf.gz \
        -V /data/output/SRR1594116.2.vcf.gz \
        -V /data/output/SRR1594117.2.vcf.gz \
        -L $intervals
