<p align="center">
  <img src="https://i.imgur.com/YPtoUI7.png" alt="Logo"/>

</p>
<h1 align="center"> @DIB-LAB/kProcessor </h1>

## Requirement

```bash
sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils libboost-all-dev wget
```

Install Corticall
```bash
git clone https://github.com/mcveanlab/Corticall
cd Corticall
ant
java -jar dist/corticall.jar -h
```

Install KMC
```bash
git clone https://github.com/refresh-bio/KMC.git
make

cd KMC
```

snakemake >=5.23

## Build

```bash
git checkout v2
cmake -Bbuild -DBUILD_USECASES=1 
cmake --build build -j4
```

## Download the Data
```bash
cd build/usecases/kProcessor_corticall
snakemake --use-conda -s download.smk
```


## run
snakemake -j2 --use-conda


## Change Dataset
edit subsample_table.csv

## Change kSize or Config
edit config.yaml


