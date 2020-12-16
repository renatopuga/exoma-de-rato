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
