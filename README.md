# Exoma de Rattus norvegicus
Variantes genéticas de amostras de Rattus norvegicus utilizando GATK4

### Amostras Utilizadas

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


### srafile.txt

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

# Pipeline


### Script getRun.sh

```bash
# criar diretorios de input e output
mkdir -p input output

# copiar, extrair, 
for srafile in $(cat $1)
do
    wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/$srafile
    ~/rat/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump --split-files $(basename $srafile)
    rm -rf $(basename $srafile)
    mkdir input/$(basename $srafile)
    mv $(basename $srafile)*.fastq input/$(basename $srafile)
    sh pipe.sh $(basename $srafile)
done
```

> Rodar getRub.sh

```bash
sh getRun.sh srafile.txt
```



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

