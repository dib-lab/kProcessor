mkdir -p inputData
cd inputData
wget https://zenodo.org/record/5906579/files/bovine_pangenome_assemblies.agc
git clone https://github.com/refresh-bio/agc
cd agc && make
cd -
agc/agc getcol -o ./ bovine_pangenome_assemblies.agc