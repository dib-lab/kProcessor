<p align="center">
  <img src="https://i.imgur.com/YPtoUI7.png" alt="Logo"/>

</p>
<h1 align="center"> @DIB-LAB/kProcessor </h1>

## Requirement

```bash
sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils libboost-all-dev wget
```

snakemake >=5.23

## Build

```bash
cmake -Bbuild -DBUILD_USECASES=1 
cmake --build build -j4
```

## Download the Data
cd build/usecases/kDifferentialExpression
./downloadData.sh


## run
snakemake -j2 --use-conda


## Change Dataset
edit subsample_table.csv

## Change kSize or Config
edit config.yaml

S