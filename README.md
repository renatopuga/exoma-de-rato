# Exoma de Rattus norvegicus

Variantes genéticas de amostras de Rattus norvegicus utilizando GATK4



### Requisitos

* docker
* wget
* sratools

## Workflow

[TOC]

### Projeto e Amostras Utilizadas (SRA)

Rattus norvegicus strain: Selectively bred alcohol-preferring (P) and nonpreferring (NP) rats (Norway rat). ExomeSeq of selectively bred alcohol-preferring (P) and nonpreferring (NP) rats.

* https://www.ncbi.nlm.nih.gov/bioproject/PRJNA262169



| Accession    | PRJNA262169                                                  |
| ------------ | ------------------------------------------------------------ |
| Data Type    | Exome                                                        |
| Scope        | Multiisolate                                                 |
| Organism     | [Rattus norvegicus](https://www.ncbi.nlm.nih.gov/taxonomy/10116)[Taxonomy ID: 10116]Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Rattus; Rattus norvegicus |
| Publications | [Zhou Z *et al.,*](https://www.ncbi.nlm.nih.gov/pubmed/24082084) "Loss of metabotropic glutamate receptor 2 escalates alcohol consumption.", *Proc Natl Acad Sci U S A*, 2013 Oct 15;110(42):16963-8 |
| Submission   | Registration date: 26-Sep-2014 **NIH/NIAAA**                 |
| Relevance    | Medical                                                      |



### Classificação

| ID           | Group                                   | Name |
| ------------ | --------------------------------------- | ---- |
| SRR1594108.2 | Selectively bred alcohol-preferring (P) | P4   |
| SRR1594109.2 | Selectively bred alcohol-preferring (P) | P5   |
| SRR1594110.2 | Selectively bred alcohol-preferring (P) | P7   |
| SRR1594111.2 | Selectively bred alcohol-preferring (P) | P8   |
| SRR1594112.2 | nonpreferring (NP)                      | NP2  |
| SRR1594113.2 | nonpreferring (NP)                      | NP4  |
| SRR1594114.2 | nonpreferring (NP)                      | NP5  |
| SRR1594115.2 | nonpreferring (NP)                      | NP6  |
| SRR1594116.2 | nonpreferring (NP)                      | NP7  |
| SRR1594117.2 | nonpreferring (NP)                      | NP8  |



# Genome

### Rattus norvegicus (rn6)

* Genes: Rat kidney disease: https://rgd.mcw.edu/rgdweb/portal/home.jsp?p=9
* dbsnp Rat: 
  * https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/archive/rat_10116/VCF/00-All.vcf.gz
  * https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/archive/rat_10116/VCF/00-All.vcf.gz.tbi
* Genome Rat UCSC (rn6): https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz



### Genoma rn6 `Index`

```bash
cd /scratch/bucket
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz
gunzip rn6.fa.gz
docker run -v /scratch/:/scratch -v $(pwd):/data/ comics/bwa bwa index /data/rn6.fa 
```



### Lista das Amostras `srafile.txt`

```bash
SRR1594108/SRR1594108.2
SRR1594109/SRR1594109.2
SRR1594110/SRR1594110.2
SRR1594111/SRR1594111.2
SRR1594112/SRR1594112.2
SRR1594113/SRR1594113.2
SRR1594114/SRR1594114.2
SRR1594115/SRR1594115.2
SRR1594116/SRR1594116.2
SRR1594117/SRR1594117.2
```



### Code `getRun.sh`

```bash
# caminho do comando fastq-dump
fastqdump="/scratch/app/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump"

# criar diretorios de input e output
mkdir -p input output

# copiar, extrair, 
for srafile in $(cat $1)
do
    wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/$srafile
    $fastqdump --split-files $(basename $srafile)
    rm -rf $(basename $srafile)
    mkdir input/$(basename $srafile)
    mv $(basename $srafile)*.fastq input/$(basename $srafile)
    sh pipe.sh $(basename $srafile)
done
```



### Executar `getRun.sh`

```bash
sh getRun.sh srafile.txt
```



### Pipeline `pipe.sh`

```bash
# scratch
genome="/scratch/bucket/rn6/rn6.fa"
dbsnp="/scratch/bucket/rn6/00-All.vcf.gz"
intervals="/scratch/bucket/rn6/rn6.interval_list"

tmp="/data/tmp"

sample=$1
LB="WES"
PL="illumina"
PU="HiSeq"

# listar o R1 de cada lanes
listR1=$(ls -1 input/$sample/*_1*.fastq)
# executar o bwa para cada lane
for i in $listR1
do
	R1=$i
	R2=$(echo $i | sed -e "s/_1/_2/g")

	docker run --user "$(id -u):$(id -g)" -v /scratch/:/scratch -v $(pwd):/data/ comics/bwa bwa mem -t 5 -M -R '@RG\tID:'$sample'\tLB:'$LB'\tSM:'$sample'\tPL:'$PL'\tPU:'$PU'' $genome /data/$R1 /data/$R2 | samtools view -F4 -Sbu -@2 - | samtools sort -m4G -@2 -o output/$sample.sorted.bam
done

docker run -v /tmp:/tmp -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk --java-options "-Djava.io.tmpdir=${tmp}  -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" MarkDuplicates \
 --TMP_DIR $tmp \
 -I /data/output/$sample.sorted.bam -O /data/output/$sample.sorted.dup.bam \
 -M /data/output/$sample.sorted.dup_metrics \
 --VALIDATION_STRINGENCY SILENT \
 --CREATE_INDEX true \

	#Rodando o conteiner de recalibrador de base para aumentar a qualidade liberando uma tabela com todos os pontos considerados
	docker run --user "$(id -u):$(id -g)" -v /scratch/:/scratch -v $(pwd):/data/ broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator -L $intervals -R $genome -I /data/output/$sample.sorted.dup.bam --known-sites $dbsnp -O /data/output/$sample.sorted.dup.recal.data.table

	docker run --user "$(id -u):$(id -g)" -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSR -R $genome -I /data/output/$sample.sorted.dup.bam -bqsr /data/output/$sample.sorted.dup.recal.data.table -L $intervals --create-output-bam-index true --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O /data/output/$sample.sorted.dup.recal.bam


rm -rf input/$sample*
rm -f output/$sample.sorted.dup.recal.data.table
rm -f output/$sample.sorted.dup.bam 
rm -f output/$sample.sorted.dup.bai
rm -f output/$sample.sorted.bam
```

