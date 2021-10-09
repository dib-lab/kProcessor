mkdir -p inputData
cd inputData

#downlad references
mkdir -p references
cd references
wget http://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/PfHB3.April2018.fasta.gz
gzip -d PfHB3.April2018.fasta.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.5_GCA_000002765/GCF_000002765.5_GCA_000002765_genomic.fna.gz
gzip -d GCF_000002765.5_GCA_000002765_genomic.fna.gz

cd ..
mkdir -p illumina/
cd illumina/
wget  ftp.sra.ebi.ac.uk/vol1/fastq/ERR019/ERR019054/ERR019054_1.fastq.gz
wget  ftp.sra.ebi.ac.uk/vol1/fastq/ERR019/ERR019054/ERR019054_2.fastq.gz
wget  ftp.sra.ebi.ac.uk/vol1/fastq/ERR019/ERR019061/ERR019061_1.fastq.gz
wget  ftp.sra.ebi.ac.uk/vol1/fastq/ERR019/ERR019061/ERR019061_2.fastq.gz



