# scratch
genome="/scratch/bucket/Rattus_norvegicus_nor6.fa"
dbsnp="/scratch/bucket/00-All.vcf.gz"
intervals="/scratch/bucket/Rattus_norvegicus_nor6.interval_list"

tmp="/data/tmp"

sample=$1
LB="WES"
PL="illumina"
PU="HiSeq"

# listar amostras
listR1=$(ls -1 input/$sample/*_1*.fastq)

# run bwa-mem 
for i in $listR1
do
	R1=$i
	R2=$(echo $i | sed -e "s/_1/_2/g")

	docker run --user "$(id -u):$(id -g)" -v /scratch/:/scratch -v $(pwd):/data/ comics/bwa bwa mem -t 5 -M -R '@RG\tID:'$sample'\tLB:'$LB'\tSM:'$sample'\tPL:'$PL'\tPU:'$PU'' $genome /data/$R1 /data/$R2 | samtools view -F4 -Sbu -@2 - | samtools sort -m4G -@2 -o output/$sample.sorted.bam
done

# run MarkDuplicates
docker run -v /tmp:/tmp -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk --java-options "-Djava.io.tmpdir=${tmp}  -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" MarkDuplicates \
 --TMP_DIR $tmp \
 -I /data/output/$sample.sorted.bam -O /data/output/$sample.sorted.dup.bam \
 -M /data/output/$sample.sorted.dup_metrics \
 --VALIDATION_STRINGENCY SILENT \
 --CREATE_INDEX true \

# run BaseRecalibrator
docker run --user "$(id -u):$(id -g)" -v /scratch/:/scratch -v $(pwd):/data/ broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator -L $intervals -R $genome -I /data/output/$sample.sorted.dup.bam --known-sites $dbsnp -O /data/output/$sample.sorted.dup.recal.data.table

# run ApplyBQSR
docker run --user "$(id -u):$(id -g)" -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSR -R $genome -I /data/output/$sample.sorted.dup.bam -bqsr /data/output/$sample.sorted.dup.recal.data.table -L $intervals --create-output-bam-index true --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O /data/output/$sample.sorted.dup.recal.bam

# remove tmp files
rm -rf input/$sample*
rm -f output/$sample.sorted.dup.recal.data.table
rm -f output/$sample.sorted.dup.bam 
rm -f output/$sample.sorted.dup.bai
rm -f output/$sample.sorted.bam
