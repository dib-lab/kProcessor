mkdir -p inputData
cd inputData
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
tail -n 1785 chr22_with_ERCC92.fa > genesERCC.fa
paste <(grep ">" genesERCC.fa |tr -d '>') <(grep ">" genesERCC.fa |tr -d '>') > genesERCC.fa.names
