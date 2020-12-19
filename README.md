# Exoma de Rattus norvegicus

Variantes genéticas de amostras de Rattus norvegicus utilizando GATK4

# Recurso Computacional
 * Ubuntu 20
 * 6 cpus (core i5)
 * 32 Gb Memória RAM
 * 256 GB SSD storage

# Requisitos

* docker

* wget

  * ```bash
    # ubuntu
    sudo apt-get install -y wget
    
    # mac
    brew install wget
    ```

* samtools

  * ```bash
    # ubuntu
    sudo apt-get install -y samtools
    
    # mac
    brew install samtools
    ```

* sratools (SRA Toolkit provides **64-bit** binary)

  | OS           | are available here                                           |
  | ------------ | ------------------------------------------------------------ |
  | **Windows**  | [sratoolkit.current-win64.zip] (http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-win64.zip) |
  | **Ubuntu**   | [sratoolkit.current-ubuntu64.tar.gz](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz) |
  | **CentOS**   | [sratoolkit.current-centos_linux64.tar.gz](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz) |
  | **Mac OS X** | [sratoolkit.current-mac64.tar.gz](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz) |


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



### Genoma nor6 `Index`, `dict`, `faidx`, `interval_list`e `BED`

```bash
# entrar no dir de referencias
cd /scratch/bucket

# download
wget -c ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz

# descompactar
gunzip Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
# renomear
mv Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rattus_norvegicus_nor6.fa 

# faidx
samtools fadix Rattus_norvegicus_nor6.fa 

# BWA index
docker run -v /scratch/:/scratch -v $(pwd):/data/ comics/bwa bwa index /data/Rattus_norvegicus_nor6.fa 

# GATK dict
docker run  -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk CreateSequenceDictionary \
	-R /data/Rattus_norvegicus_nor6.fa \
	-O /data/Rattus_norvegicus_nor6.dict
	
# GATK intereval_list (FASTA)
docker run  -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk ScatterIntervalsByNs \
      -R /data/Rattus_norvegicus_nor6.fa \
      -O /data/Rattus_norvegicus_nor6.interval_list \
      -OT ACGT

# GATK intereval_list (BED)
wget -c ftp://ftp.ensembl.org/pub/release-102/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gff3.gz

# separar regions e colunas interesse
# ensembl, insdc e gene
zgrep -v "\#" Rattus_norvegicus.Rnor_6.0.102.gff3.gz  | grep -w "ensembl\|insdc"| grep -w gene | grep  ^[0-9M] | cut -f1,4-7 | sort -Vu -k1,2 > Rnor_6.0.102.bed

docker run  -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk BedToIntervalList \
        -I /data/Rnor_6.0.102.bed  \
        -O /data/Rnor_6.0.102.interval_list \
        -SD /data/Rattus_norvegicus_nor6.dict
```



### Pipeline `mapRun.sh`

```bash
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
    sh mapRun.sh $(basename $srafile)
done
```



### Executar `getRun.sh`

```bash
sh getRun.sh srafile.txt
```



### Chamar Variantes `runCallVar.sh`

```bash
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
```

### Executar `runCallVar.sh` (~453m / ~7,5h)

```bash
sh runCallVar.sh srafile.txt
```



### GenomicsDBImport `runGenomicsDBIimport.sh` (~272m)

```bash
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
```

### Executar `runGenomicsDBIimport.sh`

```bash
sh runGenomicsDBIimport.sh
```



### GenotypeGVCFs `runGenotypeGVCFs.sh` 

```bash
genome="/scratch/bucket/Rattus_norvegicus_nor6.fa"
volume="/scratch"

docker run -v $volume:$volume -v $(pwd):/data/ broadinstitute/gatk:4.1.4.1 gatk GenotypeGVCFs \
        -V gendb:///data/genomicsdb \
        -R $genome \
        -O /data/all.samples.vcf
```



### Annotation - VEP



```bash
# docker pull vep
docker pull ensemblorg/ensembl-vep

# diretorio vep no seu computador
sudo mkdir /vep

# permissão de completo no diretorio
sudo chmod a+rwx /vep   

docker run -t -i -v /vep:/opt/vep/.vep ensemblorg/ensembl-vep
```



### Instalação `References` e `Plugins` (1.2Gb)

> -g all
>
> **Available plugins:** AncestralAllele,Blosum62,CADD,CSN,Carol,Condel,Conservation,DisGeNET,Downstream,Draw,ExAC,ExACpLI,FATHMM,FATHMM_MKL,FlagLRG,FunMotifs,G2P,GO,GeneSplicer,Gwava,HGVSIntronOffset,LD,LOVD,LoF,LoFtool,LocalID,MPC,MTR,Mastermind,MaxEntScan,NearestExonJB,NearestGene,PON_P2,Phenotypes,PostGAP,ProteinSeqs,REVEL,ReferenceQuality,SameCodon,SingleLetterAA,SpliceAI,SpliceRegion,StructuralVariantOverlap,SubsetVCF,TSSDistance,dbNSFP,dbscSNV,gnomADc,miRNA,neXtProt,satMutMPRA

```bash
# download das referencias Rnor 6.0
docker run -t -i -v /vep:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a cfp -s rattus_norvegicus -y Rnor_6.0  -g all
```



### Executar VEP

```bash
mkdir vep_data
chmod 777 vep_data

docker run -t -i -v /vep:/opt/vep/.vep  -v $(pwd):/data ensemblorg/ensembl-vep vep -i /data/all.samples.vcf -o /data/vep_data/all.samples.vep.txt  --appris --biotype --check_existing --distance 5000 --mane --sift b --species rattus_norvegicus --symbol --transcript_version --tsl --cache  --force_overwrite --tab --pick_allele --pick --pubmed --var_synonyms --variant_class --mane

less -SN vep_data/all.samples.vep.txt
```



### VEP run statistics

Veja o report completo [VEP SUMMARY](https://github.com/renatopuga/exoma-de-rato/blob/main/all.samples.vep_summary.html)

| VEP version (API)    | 102 (102)                                                    |
| -------------------- | ------------------------------------------------------------ |
| Annotation sources   | Cache: /opt/vep/.vep/rattus_norvegicus/102_Rnor_6.0; rattus_norvegicus_core_102_6 on ensembldb.ensembl.org |
| Species              | rattus_norvegicus                                            |
| Command line options | `--appris --biotype --cache --check_existing --distance 5000 --force_overwrite --input_file /data/all.samples.vcf --mane --output_file /data/vep_data/all.samples.vep --pick --pick_allele --pubmed --sift b --species rattus_norvegicus --symbol --tab --transcript_version --tsl --var_synonyms --variant_class` |
| Start time           | 2020-12-19 13:51:14                                          |
| End time             | 2020-12-19 13:54:14                                          |
| Run time             | 180 seconds                                                  |
| Input file           | /data/all.samples.vcf                                        |
| Output file          | /data/vep_data/all.samples.vep                               |

### General statistics

| Lines of input read            | 197276                       |
| ------------------------------ | ---------------------------- |
| Variants processed             | 197276                       |
| Variants filtered out          | 0                            |
| Novel / existing variants      | 121166 (61.4) / 76110 (38.6) |
| Overlapped genes               | 15071                        |
| Overlapped transcripts         | 15375                        |
| Overlapped regulatory features | -                            |

